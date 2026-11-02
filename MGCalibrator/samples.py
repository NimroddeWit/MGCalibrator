"""
samples.py

Functions for linking input names (FASTQ/BAM file names and the sample names used in
SingleM/Sylph output) to the sample IDs of the DNA mass table.

The ``sample_id`` column of the DNA mass table is the single source of truth for what a
sample is called. Nothing is guessed from the *structure* of a file name (such as "the
first two underscore-separated parts"); instead, every input name is matched against
the known sample IDs:

  1. The manual sample map (``--sample_map``), if the name occurs in it.
  2. The name equals a sample ID (ignoring directory and file extension).
  3. The name starts with a sample ID, followed by a separator ('_', '.' or '-'). If
     several IDs match, the longest one wins ('S1_rep_R1.fastq.gz' belongs to 'S1_rep'
     rather than 'S1' when both are known).

Because only the *prefix* is compared, sample IDs can contain any number of underscores,
dots or hyphens, and mate/lane/trimming tags after the ID ('_1', '_R2', '_filtered97',
'.sorted.bam', ...) never have to be interpreted.
"""

import logging
import os
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from .fileutils import strip_extension

# Accepted names for the first column of a --sample_map file. 'name' can hold any input
# name (a FASTQ/BAM file name or a sample name from SingleM/Sylph output);
# 'singlem_sample' and 'sylph_sample' are equivalent aliases kept for compatibility.
SAMPLE_MAP_NAME_COLUMNS = ["name", "singlem_sample", "sylph_sample"]

# Characters that may separate a sample ID from the rest of a longer name
SEPARATORS = "_.-"


def read_sample_map(sample_map_csv: str) -> Dict[str, str]:
    """
    Read a manual mapping of input names to sample IDs.

    Parameters
    ----------
    sample_map_csv : str
        CSV file with the columns 'name' (or 'singlem_sample' / 'sylph_sample') and
        'sample_id'. A name can be a FASTQ/BAM file name (with or without extension) or a
        sample name as written in the SingleM/Sylph output.

    Returns
    -------
    dict
        Mapping of input name -> sample ID.
    """
    sample_map_df = pd.read_csv(sample_map_csv, dtype=str)
    name_column = next((c for c in SAMPLE_MAP_NAME_COLUMNS if c in sample_map_df.columns), None)
    if name_column is None or "sample_id" not in sample_map_df.columns:
        raise ValueError(
            "Sample map file must contain the columns 'name' (or 'singlem_sample'/'sylph_sample') "
            "and 'sample_id'."
        )
    return dict(zip(sample_map_df[name_column].str.strip(), sample_map_df["sample_id"].str.strip()))


def check_sample_ids(sample_ids: Iterable[str]) -> List[str]:
    """
    Validate the sample IDs of the DNA mass table and return them as a list.

    Raises a ValueError for empty or duplicated IDs, which would make the matching of
    input files to samples ambiguous.
    """
    ids = list(sample_ids)
    if any(pd.isna(i) or not str(i).strip() for i in ids):
        raise ValueError("The DNA mass file contains an empty 'sample_id'.")
    duplicated = sorted({i for i in ids if ids.count(i) > 1})
    if duplicated:
        raise ValueError(f"Duplicated 'sample_id' in the DNA mass file: {', '.join(duplicated)}.")
    return ids


def match_sample_id(
    name: str,
    known_ids: Iterable[str],
    manual_map: Optional[Dict[str, str]] = None,
) -> Optional[str]:
    """
    Find the sample ID that an input name (file name, path or profiler sample name) belongs to.

    Directory parts of the name are ignored (Sylph reports the path of the reads).
    See the module docstring for the matching rules. Returns None if there is no match.

    Raises
    ------
    ValueError
        If the manual sample map links the name to an ID that is not in `known_ids`.
    """
    manual_map = manual_map or {}
    known_ids = set(known_ids)
    base = os.path.basename(name)

    for key in (name, base, strip_extension(base)):
        if key in manual_map:
            target = manual_map[key]
            if target not in known_ids:
                raise ValueError(
                    f"Sample map links '{key}' to '{target}', but '{target}' is not a sample_id "
                    "in the DNA mass file."
                )
            return target

    candidates = [
        sample_id for sample_id in known_ids
        if base == sample_id
        or (base.startswith(sample_id) and base[len(sample_id)] in SEPARATORS)
    ]
    return max(candidates, key=len) if candidates else None


def assign_files_to_samples(
    files: Iterable[str],
    known_ids: Iterable[str],
    manual_map: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, str], List[str]]:
    """
    Assign input files (FASTQ or BAM) to sample IDs.

    Returns
    -------
    (sample_by_file, unmatched)
        `sample_by_file` maps file path -> sample ID; `unmatched` lists the files that
        could not be assigned to any sample ID.
    """
    known_ids = set(known_ids)
    sample_by_file: Dict[str, str] = {}
    unmatched: List[str] = []
    for path in files:
        sample_id = match_sample_id(path, known_ids, manual_map)
        if sample_id is None:
            unmatched.append(path)
        else:
            sample_by_file[path] = sample_id
            logging.info(f"Assigned '{os.path.basename(path)}' -> sample '{sample_id}'")
    return sample_by_file, unmatched


def group_by_sample(sample_by_file: Dict[str, str]) -> Dict[str, List[str]]:
    """Invert a file -> sample mapping into sample -> sorted list of files."""
    grouped: Dict[str, List[str]] = {}
    for path, sample_id in sample_by_file.items():
        grouped.setdefault(sample_id, []).append(path)
    return {sample_id: sorted(paths) for sample_id, paths in grouped.items()}


def unmatched_files_message(unmatched: List[str], known_ids: Iterable[str], kind: str) -> str:
    """Build a helpful message for files that could not be assigned to a sample."""
    return (
        f"Could not assign the following {kind} file(s) to a sample: "
        f"{', '.join(os.path.basename(p) for p in unmatched)}.\n"
        "A file name must start with a sample_id of the DNA mass file, followed by '_', '.' or '-' "
        f"(known sample_ids: {', '.join(sorted(known_ids))}).\n"
        "Rename the files, or use --sample_map (CSV with the columns 'name' and 'sample_id') to link them manually."
    )


def build_sample_map(
    profiler_samples: Iterable[str],
    known_ids: Iterable[str],
    manual_map: Optional[Dict[str, str]] = None,
    source: str = "SingleM",
) -> Dict[str, str]:
    """
    Link every profiler sample name (SingleM/Sylph) to a sample ID of the DNA mass table.

    Raises a ValueError listing all profiler samples that could not be matched.

    Parameters
    ----------
    profiler_samples : iterable of str
        Sample names as they appear in the profiler (SingleM/Sylph) output.
    known_ids : iterable of str
        The sample IDs of the DNA mass table.
    manual_map : dict, optional
        Manual mapping (see read_sample_map).
    source : str
        Name of the profiler, used in log and error messages.

    Returns
    -------
    dict
        Mapping of profiler sample name -> sample ID.
    """
    manual_map = manual_map or {}
    known_ids = set(known_ids)

    mapping = {}
    unmatched = []
    for profiler_sample in sorted(set(profiler_samples)):
        target = match_sample_id(profiler_sample, known_ids, manual_map)
        if target is None:
            unmatched.append(profiler_sample)
        else:
            mapping[profiler_sample] = target
            logging.info(f"Matched {source} sample '{profiler_sample}' -> '{target}'")

    if unmatched:
        name_column = "sylph_sample" if source == "Sylph" else "singlem_sample"
        raise ValueError(
            f"Could not match the following {source} sample(s) to a sample_id of the DNA mass file: "
            f"{', '.join(unmatched)}.\n"
            f"Known sample_ids: {', '.join(sorted(known_ids))}.\n"
            "A sample name must start with a sample_id, followed by '_', '.' or '-'.\n"
            f"Use --sample_map (CSV with columns '{name_column}' and 'sample_id') to link them manually."
        )
    return mapping
