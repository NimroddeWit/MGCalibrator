"""
singlem.py

Functions for reading SingleM output tables and for linking the sample names
used by SingleM to the sample names derived from the FASTQ file names.
"""

import logging
import os
from typing import Dict, Iterable, List, Optional

import pandas as pd

from .fileutils import get_sample_name

REQUIRED_COLUMNS = ["sample", "full_coverage", "taxonomy"]

# Characters that may separate a sample name from the rest of a longer name
_SEPARATORS = "_.-"


def read_singlem_tables(singlem_files: List[str]) -> pd.DataFrame:
    """
    Read one or more SingleM output TSV files and combine them.

    Only the columns 'sample', 'full_coverage' and 'taxonomy' are used.
    A single file may contain multiple samples.

    Parameters
    ----------
    singlem_files : list of str
        Paths to SingleM output files (tab-separated, optionally gzip-compressed).

    Returns
    -------
    pandas.DataFrame
        Table with columns 'sample', 'full_coverage' (float) and 'taxonomy'.
    """
    tables = []
    for path in singlem_files:
        table = pd.read_csv(
            path,
            sep="\t",
            usecols=lambda column: column in REQUIRED_COLUMNS,
            dtype=str,
            keep_default_na=False,
            na_values=[""],
        )
        missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
        if missing:
            raise ValueError(
                f"SingleM file '{path}' is missing required column(s): {', '.join(missing)}."
            )
        logging.info(f"Read {len(table)} rows from SingleM file: {os.path.basename(path)}")
        tables.append(table[REQUIRED_COLUMNS])

    combined = pd.concat(tables, ignore_index=True)

    if combined.isna().any().any():
        raise ValueError(
            "SingleM output contains empty values in the 'sample', 'full_coverage' "
            "or 'taxonomy' column."
        )
    combined["full_coverage"] = pd.to_numeric(combined["full_coverage"], errors="raise")
    return combined


def read_sample_map(sample_map_csv: str) -> Dict[str, str]:
    """
    Read a manual mapping of SingleM sample names to MGCalibrator sample names.

    Parameters
    ----------
    sample_map_csv : str
        CSV file with the columns 'singlem_sample' and 'sample_id'.

    Returns
    -------
    dict
        Mapping of SingleM sample name -> MGCalibrator sample name.
    """
    sample_map_df = pd.read_csv(sample_map_csv, dtype=str)
    if "singlem_sample" not in sample_map_df.columns or "sample_id" not in sample_map_df.columns:
        raise ValueError("Sample map file must contain 'singlem_sample' and 'sample_id' columns.")
    return dict(zip(sample_map_df["singlem_sample"], sample_map_df["sample_id"]))


def _match_sample(
    singlem_sample: str,
    fastq_samples: Iterable[str],
    manual_map: Dict[str, str],
) -> Optional[str]:
    """
    Find the MGCalibrator sample name that belongs to a SingleM sample name.

    The following strategies are tried in order:
      1. The manual sample map.
      2. An exact match with a FASTQ-derived sample name.
      3. Applying the FASTQ naming rule (first two underscore-separated parts of
         the name, ignoring the file extension) to the SingleM sample name.
      4. A FASTQ-derived sample name that the SingleM sample name starts with,
         followed by a separator ('_', '.' or '-'). The longest such name wins.

    Returns None if no match is found.
    """
    fastq_samples = set(fastq_samples)

    if singlem_sample in manual_map:
        target = manual_map[singlem_sample]
        if target not in fastq_samples:
            raise ValueError(
                f"Sample map links SingleM sample '{singlem_sample}' to '{target}', "
                "but no FASTQ files were found for that sample."
            )
        return target

    if singlem_sample in fastq_samples:
        return singlem_sample

    derived = get_sample_name(singlem_sample)
    if derived in fastq_samples:
        return derived

    candidates = [
        name for name in fastq_samples
        if singlem_sample.startswith(name)
        and len(singlem_sample) > len(name)
        and singlem_sample[len(name)] in _SEPARATORS
    ]
    if candidates:
        return max(candidates, key=len)

    return None


def build_sample_map(
    singlem_samples: Iterable[str],
    fastq_samples: Iterable[str],
    manual_map: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Link every SingleM sample name to a sample name derived from the FASTQ files.

    Raises a ValueError listing all SingleM samples that could not be matched.

    Returns
    -------
    dict
        Mapping of SingleM sample name -> MGCalibrator sample name.
    """
    manual_map = manual_map or {}
    fastq_samples = set(fastq_samples)

    mapping = {}
    unmatched = []
    for singlem_sample in sorted(set(singlem_samples)):
        target = _match_sample(singlem_sample, fastq_samples, manual_map)
        if target is None:
            unmatched.append(singlem_sample)
        else:
            mapping[singlem_sample] = target
            logging.info(f"Matched SingleM sample '{singlem_sample}' -> '{target}'")

    if unmatched:
        raise ValueError(
            "Could not match the following SingleM sample(s) to FASTQ files: "
            f"{', '.join(unmatched)}.\n"
            f"Sample names derived from the FASTQ files: {', '.join(sorted(fastq_samples))}.\n"
            "Use --sample_map (CSV with columns 'singlem_sample' and 'sample_id') to link them manually."
        )
    return mapping


def singlem_to_depth_table(singlem_df: pd.DataFrame, sample_map: Dict[str, str]) -> pd.DataFrame:
    """
    Convert a SingleM table into the depth table used by MGCalibrator.

    SingleM 'sample' names are replaced by the MGCalibrator sample names,
    'taxonomy' becomes 'reference' and 'full_coverage' becomes 'depth'.
    If multiple rows end up with the same sample and reference (for example
    because R1 and R2 were analysed by SingleM separately), their depths are summed.

    Returns
    -------
    pandas.DataFrame
        Table with columns 'sample', 'reference' and 'depth'.
    """
    depth_table = pd.DataFrame({
        "sample": singlem_df["sample"].map(sample_map),
        "reference": singlem_df["taxonomy"],
        "depth": singlem_df["full_coverage"],
    })

    n_duplicates = int(depth_table.duplicated(subset=["sample", "reference"]).sum())
    if n_duplicates:
        logging.warning(
            f"{n_duplicates} rows share a sample and taxonomy with another row; "
            "their full_coverage values are summed."
        )

    return (
        depth_table
        .groupby(["sample", "reference"], sort=False, as_index=False)["depth"]
        .sum()
    )