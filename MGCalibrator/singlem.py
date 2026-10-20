"""
singlem.py

Functions for reading SingleM output tables (standard taxonomic profile format,
as written by `singlem pipe --taxonomic-profile` or `singlem condense`).
Sample name matching is implemented in samples.py.
"""

import logging
import os
from typing import Dict, List

import pandas as pd

# Re-exported for backwards compatibility
from .samples import build_sample_map, read_sample_map  # noqa: F401

REQUIRED_COLUMNS = ["sample", "coverage", "taxonomy"]


def _is_species(taxonomy: str) -> bool:
    """Return True if the last rank of a SingleM taxonomy string is a species ('s__<name>')."""
    last_rank = taxonomy.split(";")[-1].strip()
    return last_rank.startswith("s__") and len(last_rank) > len("s__")


def read_singlem_tables(singlem_files: List[str]) -> pd.DataFrame:
    """
    Read one or more SingleM taxonomic profiles (TSV) and combine them.

    Only the columns 'sample', 'coverage' and 'taxonomy' are used, and only rows at
    species level (taxonomy ending in 's__<name>') are kept. At species level the
    (unfilled) coverage is equal to the full coverage, because a species has no
    descendants. A single file may contain multiple samples.

    Parameters
    ----------
    singlem_files : list of str
        Paths to SingleM output files (tab-separated, optionally gzip-compressed).

    Returns
    -------
    pandas.DataFrame
        Table with columns 'sample', 'coverage' (float) and 'taxonomy'.
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
                f"SingleM file '{path}' is missing required column(s): {', '.join(missing)}. "
                "Expected the standard SingleM taxonomic profile with the columns "
                "sample, coverage and taxonomy."
            )
        logging.info(f"Read {len(table)} rows from SingleM file: {os.path.basename(path)}")
        tables.append(table[REQUIRED_COLUMNS])

    combined = pd.concat(tables, ignore_index=True)

    if combined.isna().any().any():
        raise ValueError(
            "SingleM output contains empty values in the 'sample', 'coverage' "
            "or 'taxonomy' column."
        )
    combined["coverage"] = pd.to_numeric(combined["coverage"], errors="raise")

    is_species = combined["taxonomy"].map(_is_species)
    logging.info(
        f"Keeping {int(is_species.sum())} species-level rows "
        f"(dropping {int((~is_species).sum())} rows at higher taxonomic levels)."
    )
    species = combined.loc[is_species].reset_index(drop=True)
    if species.empty:
        raise ValueError("No species-level rows (taxonomy ending in 's__<name>') found in the SingleM output.")
    return species


def singlem_to_depth_table(singlem_df: pd.DataFrame, sample_map: Dict[str, str]) -> pd.DataFrame:
    """
    Convert a SingleM table into the depth table used by MGCalibrator.

    SingleM 'sample' names are replaced by the MGCalibrator sample names,
    'taxonomy' (the full taxonomic string) becomes 'reference' and 'coverage'
    becomes 'depth'. If multiple rows end up with the same sample and reference
    (for example because R1 and R2 were analysed by SingleM separately), their
    depths are summed.

    Returns
    -------
    pandas.DataFrame
        Table with columns 'sample', 'reference' and 'depth'.
    """
    depth_table = pd.DataFrame({
        "sample": singlem_df["sample"].map(sample_map),
        "reference": singlem_df["taxonomy"],
        "depth": singlem_df["coverage"],
    })

    n_duplicates = int(depth_table.duplicated(subset=["sample", "reference"]).sum())
    if n_duplicates:
        logging.warning(
            f"{n_duplicates} rows share a sample and taxonomy with another row; "
            "their coverage values are summed."
        )

    return (
        depth_table
        .groupby(["sample", "reference"], sort=False, as_index=False)["depth"]
        .sum()
    )
