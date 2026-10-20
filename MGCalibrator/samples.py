"""
samples.py

Functions for linking the sample names used in profiler output files (SingleM,
Sylph) to the sample names derived from the FASTQ file names.
"""

import logging
import os
from typing import Dict, Iterable, Optional

import pandas as pd

from .fileutils import get_sample_name

# Accepted names for the first column of a --sample_map file. 'singlem_sample' and
# 'sylph_sample' are equivalent: both hold the sample name as it appears in the
# profiler output.
SAMPLE_MAP_NAME_COLUMNS = ["singlem_sample", "sylph_sample"]

# Characters that may separate a sample name from the rest of a longer name
_SEPARATORS = "_.-"


def read_sample_map(sample_map_csv: str) -> Dict[str, str]:
    """
    Read a manual mapping of profiler sample names to MGCalibrator sample names.

    Parameters
    ----------
    sample_map_csv : str
        CSV file with the columns 'singlem_sample' (or 'sylph_sample') and 'sample_id'.

    Returns
    -------
    dict
        Mapping of profiler sample name -> MGCalibrator sample name.
    """
    sample_map_df = pd.read_csv(sample_map_csv, dtype=str)
    name_column = next((c for c in SAMPLE_MAP_NAME_COLUMNS if c in sample_map_df.columns), None)
    if name_column is None or "sample_id" not in sample_map_df.columns:
        raise ValueError(
            "Sample map file must contain the columns 'singlem_sample' (or 'sylph_sample') and 'sample_id'."
        )
    return dict(zip(sample_map_df[name_column], sample_map_df["sample_id"]))


def _match_sample(
    profiler_sample: str,
    fastq_samples: Iterable[str],
    manual_map: Dict[str, str],
) -> Optional[str]:
    """
    Find the MGCalibrator sample name that belongs to a profiler sample name.

    Directory parts of the name are ignored (Sylph reports the path of the reads).
    The following strategies are tried in order:
      1. The manual sample map (full name or base name).
      2. An exact match with a FASTQ-derived sample name.
      3. Applying the FASTQ naming rule (first two underscore-separated parts of
         the name, ignoring the file extension) to the profiler sample name.
      4. A FASTQ-derived sample name that the profiler sample name starts with,
         followed by a separator ('_', '.' or '-'). The longest such name wins.

    Returns None if no match is found.
    """
    fastq_samples = set(fastq_samples)
    base = os.path.basename(profiler_sample)

    for key in (profiler_sample, base):
        if key in manual_map:
            target = manual_map[key]
            if target not in fastq_samples:
                raise ValueError(
                    f"Sample map links sample '{key}' to '{target}', "
                    "but no FASTQ files were found for that sample."
                )
            return target

    if base in fastq_samples:
        return base

    derived = get_sample_name(base)
    if derived in fastq_samples:
        return derived

    candidates = [
        name for name in fastq_samples
        if base.startswith(name)
        and len(base) > len(name)
        and base[len(name)] in _SEPARATORS
    ]
    if candidates:
        return max(candidates, key=len)

    return None


def build_sample_map(
    profiler_samples: Iterable[str],
    fastq_samples: Iterable[str],
    manual_map: Optional[Dict[str, str]] = None,
    source: str = "SingleM",
) -> Dict[str, str]:
    """
    Link every profiler sample name to a sample name derived from the FASTQ files.

    Raises a ValueError listing all profiler samples that could not be matched.

    Parameters
    ----------
    profiler_samples : iterable of str
        Sample names as they appear in the profiler (SingleM/Sylph) output.
    fastq_samples : iterable of str
        Sample names derived from the FASTQ file names.
    manual_map : dict, optional
        Manual mapping (see read_sample_map).
    source : str
        Name of the profiler, used in log and error messages.

    Returns
    -------
    dict
        Mapping of profiler sample name -> MGCalibrator sample name.
    """
    manual_map = manual_map or {}
    fastq_samples = set(fastq_samples)

    mapping = {}
    unmatched = []
    for profiler_sample in sorted(set(profiler_samples)):
        target = _match_sample(profiler_sample, fastq_samples, manual_map)
        if target is None:
            unmatched.append(profiler_sample)
        else:
            mapping[profiler_sample] = target
            logging.info(f"Matched {source} sample '{profiler_sample}' -> '{target}'")

    if unmatched:
        name_column = "sylph_sample" if source == "Sylph" else "singlem_sample"
        raise ValueError(
            f"Could not match the following {source} sample(s) to FASTQ files: "
            f"{', '.join(unmatched)}.\n"
            f"Sample names derived from the FASTQ files: {', '.join(sorted(fastq_samples))}.\n"
            f"Use --sample_map (CSV with columns '{name_column}' and 'sample_id') to link them manually."
        )
    return mapping
