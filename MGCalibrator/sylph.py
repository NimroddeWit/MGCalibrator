"""
sylph.py

Functions for reading Sylph results for MGCalibrator.

Two kinds of files are supported:

* ``.sylphmpa`` files written by ``sylph-tax taxprof`` (required). These contain the
  full taxonomic string of every detected genome and, for strain-level rows
  (``t__<genome>``), the coverage (column "Coverage (if strain-level)", which is the
  ``Eff_cov`` or ``True_cov`` column of Sylph). ``sylph-tax taxprof`` writes one
  ``.sylphmpa`` file per sample; the sample name is stored in the ``#SampleID`` header
  line of each file.
* ``sylph profile`` output TSV files (optional). Their ``Lambda_5-95_percentile``
  column gives a confidence interval around the coverage (``Eff_lambda``), which
  MGCalibrator uses as ``error_down`` and ``error_up``.

Sample name matching is implemented in samples.py.
"""

import gzip
import logging
import os
import re
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd

SYLPHMPA_REQUIRED_COLUMNS = ["clade_name", "Coverage (if strain-level)"]
PROFILE_REQUIRED_COLUMNS = ["Sample_file", "Genome_file", "Contig_name", "Eff_lambda", "Lambda_5-95_percentile"]

STRAIN_PREFIX = "t__"

_NUMBER = r"[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"
_CI_PATTERN = re.compile(rf"^\s*({_NUMBER})-({_NUMBER})\s*$")


# ----------------------------------------------------------------------
# .sylphmpa files (taxonomy + coverage)
# ----------------------------------------------------------------------
def _open_text(path: str):
    """Open a plain or gzip-compressed text file."""
    if path.endswith((".gz", ".gzip")):
        return gzip.open(path, "rt")
    return open(path, "r")


def _read_sylphmpa(path: str) -> pd.DataFrame:
    """
    Read the strain-level rows of one .sylphmpa file.

    Returns
    -------
    pandas.DataFrame
        Columns 'sample' (the sample name from the '#SampleID' header line),
        'reference' (full taxonomic string, exactly as written by sylph-tax),
        'depth' (coverage), 'accession' (the identifier after 't__') and 'source_file'.
    """
    sample_name = None
    n_comment_lines = 0
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            n_comment_lines += 1
            if line.startswith("#SampleID"):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) > 1 and fields[1]:
                    sample_name = fields[1]

    table = pd.read_csv(
        path,
        sep="\t",
        skiprows=n_comment_lines,
        dtype=str,
        keep_default_na=False,
        na_values=[""],
    )
    missing = [c for c in SYLPHMPA_REQUIRED_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(
            f"'{path}' is not a .sylphmpa file with coverage (missing column(s): {', '.join(missing)}). "
            "Use the output of `sylph-tax taxprof` without the --pavian option."
        )

    if sample_name is None:
        base = os.path.basename(path)
        sample_name = re.sub(r"\.sylphmpa(\.gz)?$", "", base)
        logging.warning(
            f"No '#SampleID' header line found in '{path}'; using the file name ('{sample_name}') as sample name."
        )

    clade = table["clade_name"].astype(str)
    last_rank = clade.str.rsplit("|", n=1).str[-1]
    is_strain = last_rank.str.startswith(STRAIN_PREFIX)
    coverage = pd.to_numeric(table["Coverage (if strain-level)"], errors="coerce")

    keep = is_strain & coverage.notna()
    strains = pd.DataFrame({
        "sample": sample_name,
        "reference": clade[keep],
        "depth": coverage[keep],
        "accession": last_rank[keep].str[len(STRAIN_PREFIX):],
        "source_file": os.path.basename(path),
    }).reset_index(drop=True)

    logging.info(
        f"Read {len(strains)} strain-level rows (of {len(table)} rows) for sample '{sample_name}' "
        f"from {os.path.basename(path)}"
    )
    return strains


def read_sylphmpa_files(sylphmpa_files: List[str]) -> pd.DataFrame:
    """
    Read one or more .sylphmpa files and keep only the strain-level rows.

    Files for the same sample (for example from separate Sylph runs with different
    databases) are combined; they are recognised by the sample name in the file header.

    Returns
    -------
    pandas.DataFrame
        Columns 'sample', 'reference', 'depth', 'accession' and 'source_file'.
    """
    strains = pd.concat([_read_sylphmpa(path) for path in sylphmpa_files], ignore_index=True)
    if strains.empty:
        raise ValueError(
            "No strain-level rows ('t__') with coverage found in the .sylphmpa files. "
            "Was `sylph-tax taxprof` run with the right taxonomy metadata (-t)?"
        )

    n_no_taxonomy = int(strains["reference"].str.startswith("NO_TAXONOMY").sum())
    if n_no_taxonomy:
        logging.warning(
            f"{n_no_taxonomy} strain-level rows have no taxonomy (NO_TAXONOMY). "
            "Check that the taxonomy metadata (-t) used for `sylph-tax taxprof` matches the Sylph database."
        )

    n_duplicates = int(strains.duplicated(subset=["sample", "reference"]).sum())
    if n_duplicates:
        logging.warning(
            f"{n_duplicates} strain-level rows occur more than once for the same sample "
            "(e.g. the same genome in two databases). All rows are kept; they are not summed."
        )
    return strains


# ----------------------------------------------------------------------
# sylph profile TSV files (confidence intervals)
# ----------------------------------------------------------------------
def _genome_file_to_gcf_acc(file_name: str) -> str:
    if "ASM" in file_name:
        return file_name.split("/")[-1].split("_ASM")[0]
    return file_name.split("/")[-1].split("_genomic")[0]


def _accession_candidates(genome_file: str, contig_name: str) -> Set[str]:
    """
    Return all identifiers that sylph-tax may have used as the genome accession ('t__<accession>')
    for a row of a sylph profile. This mirrors how sylph-tax links Sylph genomes to its taxonomy
    metadata, so a .sylphmpa strain row can be linked to its row in the Sylph profile.
    """
    contig_name = contig_name.strip()
    candidates: Set[str] = set()

    trimmed_genome = genome_file.split("/")[-1]
    trimmed_contig = contig_name.split()[0] if contig_name else ""

    for genome, contig, genome_glob in (
        (_genome_file_to_gcf_acc(genome_file), contig_name.split(" ")[0].split("|")[0], trimmed_genome.split(".")[0]),
        (trimmed_genome, trimmed_contig, trimmed_genome),
    ):
        candidates.update({
            genome,
            genome + ".gz",
            genome_glob,
            contig,
            genome.split(".fa")[0],
            genome.split(".fasta")[0],
            genome.split(".fna")[0],
        })

    # Fallback accession used by sylph-tax when no taxonomy was found
    candidates.add(f"{trimmed_genome}:{trimmed_contig}")
    candidates.discard("")
    return candidates


def _parse_ci(value: str):
    """Parse a 'low-high' confidence interval string. Returns (nan, nan) for 'NA-NA'."""
    match = _CI_PATTERN.match(str(value))
    if not match:
        return np.nan, np.nan
    return float(match.group(1)), float(match.group(2))


def read_sylph_profiles(profile_files: List[str]) -> pd.DataFrame:
    """
    Read one or more `sylph profile` output TSV files (single- or multi-sample).

    Files without the expected Sylph columns are skipped with a warning.

    Returns
    -------
    pandas.DataFrame
        Columns 'sample' (Sample_file), 'candidates' (set of possible genome accessions),
        'coverage' (Eff_cov or True_cov), 'ci_low' and 'ci_high' (5th and 95th percentile of
        the lambda estimate, on the scale of the coverage), 'source_file'.
    """
    tables = []
    for path in profile_files:
        with _open_text(path) as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
        coverage_column = next((c for c in ("Eff_cov", "True_cov") if c in header), None)
        if coverage_column is None or any(c not in header for c in PROFILE_REQUIRED_COLUMNS):
            logging.warning(f"Skipping '{path}': not a `sylph profile` output file with taxonomic abundances.")
            continue

        # Restrict to the header's columns (a tab in a contig name would otherwise add fields)
        table = pd.read_csv(
            path, sep="\t", usecols=range(len(header)), dtype=str, keep_default_na=False,
        )
        lambda_value = pd.to_numeric(table["Eff_lambda"], errors="coerce")  # NaN for LOW/HIGH
        coverage = pd.to_numeric(table[coverage_column], errors="coerce")
        ci = table["Lambda_5-95_percentile"].map(_parse_ci)
        ci_low = ci.str[0].astype(float)
        ci_high = ci.str[1].astype(float)

        # With `sylph profile -u` the reported coverage is the true coverage, while the
        # confidence interval is on the effective-coverage (lambda) scale: rescale it.
        if coverage_column == "True_cov":
            scale = coverage / lambda_value
            ci_low, ci_high = ci_low * scale, ci_high * scale

        tables.append(pd.DataFrame({
            "sample": table["Sample_file"].astype(str),
            "candidates": [
                _accession_candidates(g, c) for g, c in zip(table["Genome_file"], table["Contig_name"])
            ],
            "coverage": coverage,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "source_file": os.path.basename(path),
        }))
        logging.info(f"Read {len(table)} rows from Sylph profile: {os.path.basename(path)}")

    if not tables:
        raise ValueError("None of the supplied files is a `sylph profile` output file.")
    return pd.concat(tables, ignore_index=True)


def add_confidence_intervals(strains: pd.DataFrame, profiles: pd.DataFrame) -> pd.DataFrame:
    """
    Add the Sylph confidence intervals to the strain-level rows of the .sylphmpa files.

    A .sylphmpa row is linked to its row in the Sylph profile by sample name (the
    '#SampleID' of the .sylphmpa equals 'Sample_file' of the profile) and genome
    accession (the 't__<accession>' of the .sylphmpa). The coverage values of both files
    are compared as a consistency check. The interval is only available where Sylph could
    estimate lambda (not for LOW/HIGH coverage); other rows get NaN.

    Returns
    -------
    pandas.DataFrame
        'strains' with the extra columns 'lower_ci' and 'upper_ci'.
    """
    index: Dict[tuple, List[int]] = {}
    for row_number, (sample, candidates) in enumerate(zip(profiles["sample"], profiles["candidates"])):
        for candidate in candidates:
            index.setdefault((sample, candidate), []).append(row_number)

    coverage = profiles["coverage"].to_numpy()
    ci_low = profiles["ci_low"].to_numpy()
    ci_high = profiles["ci_high"].to_numpy()

    lower, upper = [], []
    n_unmatched = n_inconsistent = 0
    for sample, accession, depth in zip(strains["sample"], strains["accession"], strains["depth"]):
        rows = index.get((sample, accession), [])
        # Several profile rows can share an identifier: keep those with the same coverage
        rows = [r for r in rows if np.isclose(coverage[r], depth, rtol=1e-3, atol=1e-3)] or rows
        if not rows:
            n_unmatched += 1
            lower.append(np.nan)
            upper.append(np.nan)
            continue
        row = rows[0]
        if not np.isclose(coverage[row], depth, rtol=1e-3, atol=1e-3):
            n_inconsistent += 1
            lower.append(np.nan)
            upper.append(np.nan)
            continue
        lower.append(ci_low[row])
        upper.append(ci_high[row])

    result = strains.copy()
    result["lower_ci"] = lower
    result["upper_ci"] = upper

    n_with_ci = int(result["lower_ci"].notna().sum())
    logging.info(
        f"Added Sylph confidence intervals to {n_with_ci} of {len(result)} strain-level rows. "
        "Rows without interval have a coverage that is too low or too high for Sylph to estimate one."
    )
    if n_unmatched:
        logging.warning(
            f"{n_unmatched} strain-level rows could not be linked to a row in the Sylph profile(s) "
            "(sample or genome not found); their confidence interval is NA. "
            "Check that the profile files are the sylph output that the .sylphmpa files were made from."
        )
    if n_inconsistent:
        logging.warning(
            f"{n_inconsistent} strain-level rows have a different coverage in the Sylph profile than in the "
            ".sylphmpa file; their confidence interval is NA. "
            "The files probably come from different Sylph runs."
        )
    return result


# ----------------------------------------------------------------------
# Conversion to the MGCalibrator depth table
# ----------------------------------------------------------------------
def sylph_to_depth_table(strains: pd.DataFrame, sample_map: Dict[str, str]) -> pd.DataFrame:
    """
    Convert strain-level Sylph rows into the depth table used by MGCalibrator.

    The Sylph sample names are replaced by the MGCalibrator sample names. 'reference' is the
    full taxonomic string, 'depth' the coverage, and 'lower_ci' / 'upper_ci' the confidence
    interval (NaN if no profile was supplied or no interval is available).

    Returns
    -------
    pandas.DataFrame
        Columns 'sample', 'reference', 'depth', 'lower_ci' and 'upper_ci'.
    """
    table = pd.DataFrame({
        "sample": strains["sample"].map(sample_map),
        "reference": strains["reference"],
        "depth": strains["depth"],
        "lower_ci": strains["lower_ci"] if "lower_ci" in strains else np.nan,
        "upper_ci": strains["upper_ci"] if "upper_ci" in strains else np.nan,
    })
    return table.reset_index(drop=True)
