"""
sylph.py

Functions for reading Sylph results for MGCalibrator.

Two kinds of files are used:

* ``.sylphmpa`` files written by ``sylph-tax taxprof`` (required). These contain the
  full taxonomic string of every detected genome and, for strain-level rows
  (``t__<genome>``), the coverage (column "Coverage (if strain-level)"). ``sylph-tax
  taxprof`` writes one ``.sylphmpa`` file per sample; the sample name is stored in the
  ``#SampleID`` header line of each file.
* ``sylph profile`` output TSV files (optional). They are only used to obtain the
  estimated number of reads per genome (column ``Sequence_abundance``), which is used to
  calculate a confidence interval around the coverage. The confidence intervals that
  Sylph itself reports (``Lambda_5-95_percentile``) are NOT used.

Sylph has to be run with ``--estimate-unknown`` (``-u``): the reported coverage is then the
true coverage (``True_cov``), which is much more accurate than the effective coverage
(``Eff_cov``). To obtain confidence intervals, Sylph also has to be run with
``--estimate-read-counts`` (Sylph >= 0.9.0; this implies ``-u``), which makes Sylph write
the estimated read counts instead of percentages to the ``Sequence_abundance`` column.

Confidence interval of the coverage
-----------------------------------
The number of reads of a genome is Poisson distributed, and the coverage is directly
proportional to the number of reads. The relative 95% interval that is simulated around
the read count (``numpy.random.poisson``) is therefore also the relative interval around
the coverage: for a coverage of 0.5 with 1,000 reads the interval is about 0.47-0.53, with
100,000 reads about 0.497-0.503.

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
PROFILE_REQUIRED_COLUMNS = ["Sample_file", "Genome_file", "Contig_name", "Sequence_abundance"]

STRAIN_PREFIX = "t__"

# Confidence interval of the coverage: percentiles of the simulated Poisson read counts
CI_LOWER_PERCENTILE = 2.5
CI_UPPER_PERCENTILE = 97.5
DEFAULT_N_SIMULATIONS = 500

# Without --estimate-read-counts, Sequence_abundance is a percentage and sums to at most 100
# per sample (exactly 100 without -u, the percentage of classified reads with -u).
PERCENTAGE_SUM_THRESHOLD = 100.5


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
# sylph profile TSV files (estimated read counts)
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


def read_sylph_profiles(profile_files: List[str]) -> pd.DataFrame:
    """
    Read one or more `sylph profile` output TSV files (single- or multi-sample).

    Files without the expected Sylph columns are skipped with a warning. The
    'Sequence_abundance' column is only usable as number of reads if Sylph was run with
    `--estimate-read-counts`; samples whose values look like percentages (sum <= 100) get
    NaN as read count.

    Returns
    -------
    pandas.DataFrame
        Columns 'sample' (Sample_file), 'candidates' (set of possible genome accessions),
        'coverage' (True_cov, or Eff_cov), 'read_counts' (estimated number of reads of the
        genome) and 'source_file'.

    Raises
    ------
    ValueError
        If none of the files is a `sylph profile` output file, or if none of the samples
        has read counts (Sylph was run without `--estimate-read-counts`).
    """
    tables = []
    for path in profile_files:
        with _open_text(path) as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
        coverage_column = next((c for c in ("True_cov", "Eff_cov") if c in header), None)
        if coverage_column is None or any(c not in header for c in PROFILE_REQUIRED_COLUMNS):
            logging.warning(f"Skipping '{path}': not a `sylph profile` output file with taxonomic abundances.")
            continue
        if coverage_column == "Eff_cov":
            logging.warning(
                f"'{os.path.basename(path)}' contains the effective coverage (Eff_cov): Sylph was run without "
                "--estimate-unknown (-u). Always run Sylph with --estimate-unknown (or --estimate-read-counts, "
                "which implies it) to get the true coverage (True_cov), which is much more accurate."
            )

        # Restrict to the header's columns (a tab in a contig name would otherwise add fields)
        table = pd.read_csv(
            path, sep="\t", usecols=range(len(header)), dtype=str, keep_default_na=False,
        )
        sample = table["Sample_file"].astype(str)
        read_counts = pd.to_numeric(table["Sequence_abundance"], errors="coerce")
        looks_like_percentage = read_counts.groupby(sample).transform("sum") <= PERCENTAGE_SUM_THRESHOLD

        tables.append(pd.DataFrame({
            "sample": sample,
            "candidates": [
                _accession_candidates(g, c) for g, c in zip(table["Genome_file"], table["Contig_name"])
            ],
            "coverage": pd.to_numeric(table[coverage_column], errors="coerce"),
            "read_counts": read_counts.mask(looks_like_percentage),
            "looks_like_percentage": looks_like_percentage,
            "source_file": os.path.basename(path),
        }))
        logging.info(f"Read {len(table)} rows from Sylph profile: {os.path.basename(path)}")

    if not tables:
        raise ValueError("None of the supplied files is a `sylph profile` output file.")
    profiles = pd.concat(tables, ignore_index=True)

    if profiles["looks_like_percentage"].all():
        raise ValueError(
            "The 'Sequence_abundance' column of the Sylph profile(s) contains percentages instead of read counts. "
            "Confidence intervals are calculated from the number of reads, so run `sylph profile` with "
            "--estimate-read-counts (Sylph >= 0.9.0, implies --estimate-unknown), or leave out "
            "--sylph_profile_folder."
        )
    percentage_samples = sorted(profiles.loc[profiles["looks_like_percentage"], "sample"].unique())
    if percentage_samples:
        logging.warning(
            "The 'Sequence_abundance' of the following sample(s) looks like percentages instead of read counts, "
            f"so no confidence interval is calculated for them: {', '.join(percentage_samples)}. "
            "Was `sylph profile` run with --estimate-read-counts for all samples?"
        )
    return profiles.drop(columns="looks_like_percentage")


def add_read_counts(strains: pd.DataFrame, profiles: pd.DataFrame) -> pd.DataFrame:
    """
    Add the estimated number of reads per genome from the Sylph profile(s) to the
    strain-level rows of the .sylphmpa files.

    A .sylphmpa row is linked to its row in the Sylph profile by sample name (the
    '#SampleID' of the .sylphmpa equals 'Sample_file' of the profile) and genome
    accession (the 't__<accession>' of the .sylphmpa). The coverage values of both files
    are compared as a consistency check; rows that cannot be linked, or whose coverage
    differs, get NaN.

    Returns
    -------
    pandas.DataFrame
        'strains' with the extra column 'read_counts'.
    """
    index: Dict[tuple, List[int]] = {}
    for row_number, (sample, candidates) in enumerate(zip(profiles["sample"], profiles["candidates"])):
        for candidate in candidates:
            index.setdefault((sample, candidate), []).append(row_number)

    coverage = profiles["coverage"].to_numpy()
    read_counts = profiles["read_counts"].to_numpy()

    counts = []
    n_unmatched = n_inconsistent = 0
    for sample, accession, depth in zip(strains["sample"], strains["accession"], strains["depth"]):
        rows = index.get((sample, accession), [])
        # Several profile rows can share an identifier: keep those with the same coverage
        rows = [r for r in rows if np.isclose(coverage[r], depth, rtol=1e-3, atol=1e-3)] or rows
        if not rows:
            n_unmatched += 1
            counts.append(np.nan)
            continue
        row = rows[0]
        if not np.isclose(coverage[row], depth, rtol=1e-3, atol=1e-3):
            n_inconsistent += 1
            counts.append(np.nan)
            continue
        counts.append(read_counts[row])

    result = strains.copy()
    result["read_counts"] = counts

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


def poisson_coverage_ci(
    coverage,
    read_counts,
    n_simulations: int = DEFAULT_N_SIMULATIONS,
    rng: Optional[np.random.Generator] = None,
    batch_size: int = 5000,
):
    """
    95% confidence interval around the coverage, derived from the Poisson noise of the read count.

    The number of reads is simulated `n_simulations` times from a Poisson distribution with
    the estimated read count as mean. The 2.5th and 97.5th percentile of the simulated
    counts, relative to the estimated count, are applied to the coverage (the coverage is
    proportional to the number of reads). Because this interval only depends on the number
    of reads, it is small for genomes with many reads and wide for genomes with few reads.

    The interval is returned on the scale of the coverage, before any Qubit error is
    added: the caller turns it into relative errors and combines these with the Qubit
    error (see cli.py). The interval always contains the coverage itself.

    Parameters
    ----------
    coverage : array-like
        Coverage (depth) of every genome.
    read_counts : array-like
        Estimated number of reads of every genome (same length as `coverage`).
    n_simulations : int
        Number of Poisson draws per genome.
    rng : numpy.random.Generator, optional
        Random number generator (pass a seeded one for reproducible results).
    batch_size : int
        Number of genomes simulated at once (limits memory use).

    Returns
    -------
    (lower, upper) : tuple of numpy.ndarray
        Lower and upper bound of the coverage; NaN where the coverage or the read count
        is missing or the read count is not positive.
    """
    coverage = np.asarray(coverage, dtype=float)
    read_counts = np.asarray(read_counts, dtype=float)
    if coverage.shape != read_counts.shape:
        raise ValueError("coverage and read_counts must have the same length.")
    if n_simulations < 1:
        raise ValueError("n_simulations must be at least 1.")
    rng = np.random.default_rng() if rng is None else rng

    lower = np.full(coverage.shape, np.nan)
    upper = np.full(coverage.shape, np.nan)
    valid = np.flatnonzero(
        np.isfinite(coverage) & np.isfinite(read_counts) & (read_counts > 0) & (coverage >= 0)
    )
    for start in range(0, valid.size, batch_size):
        rows = valid[start:start + batch_size]
        expected = read_counts[rows]
        simulated = rng.poisson(expected[:, None], size=(rows.size, n_simulations))
        q_low, q_high = np.percentile(simulated, [CI_LOWER_PERCENTILE, CI_UPPER_PERCENTILE], axis=1)
        lower[rows] = coverage[rows] * q_low / expected
        upper[rows] = coverage[rows] * q_high / expected

    # Guard against simulation noise for very small read counts
    return np.minimum(lower, coverage), np.maximum(upper, coverage)


def add_poisson_confidence_intervals(
    strains: pd.DataFrame,
    n_simulations: int = DEFAULT_N_SIMULATIONS,
    rng: Optional[np.random.Generator] = None,
) -> pd.DataFrame:
    """
    Add the confidence interval of the coverage ('lower_ci' and 'upper_ci') to the
    strain-level rows, calculated from their 'read_counts' (see add_read_counts and
    poisson_coverage_ci). Rows without read count get NaN.
    """
    result = strains.copy()
    result["lower_ci"], result["upper_ci"] = poisson_coverage_ci(
        result["depth"], result["read_counts"], n_simulations=n_simulations, rng=rng
    )
    n_with_ci = int(result["lower_ci"].notna().sum())
    logging.info(
        f"Calculated Poisson confidence intervals for {n_with_ci} of {len(result)} strain-level rows "
        f"({n_simulations} simulations per row). Rows without an estimated read count have no interval."
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
