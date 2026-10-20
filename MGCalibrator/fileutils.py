import os
from glob import glob
from typing import Dict, List, Optional


def list_bam_files(
    folder: str,
    extensions: Optional[List[str]] = None,
    suffixes: Optional[List[str]] = None
) -> List[str]:
    """
    List BAM files in a specified folder based on given extensions and suffixes.

    This function searches for BAM files (or compressed BAM files) within a folder.
    You can filter by specific file suffixes (e.g., '_R1', '_R2') after removing
    the file extension.

    Parameters
    ----------
    folder : str
        The path to the directory containing BAM files.
    extensions : list of str, optional
        List of file extensions to include (default: [".sorted.bam", ".sort.bam"]).
        Compressed variants (".gz", ".gzip") are also automatically included.
    suffixes : list of str, optional
        List of suffixes to match at the end of the filename (before the extension).
        If provided, only files whose base names end with one of the given suffixes
        are returned.
    singleton_files : list of str, optional
        Reserved for future use; currently not used in filtering.

    Returns
    -------
    list of str
        Sorted list of matching file paths.
    """
    if extensions is None:
        extensions = [".sorted.bam", ".sort.bam"]

    # Expand extensions to include compressed versions
    all_extensions = []
    for ext in extensions:
        all_extensions.extend([ext, f"{ext}.gz", f"{ext}.gzip"])

    # Find all files with matching extensions
    files = []
    for ext in all_extensions:
        files.extend(glob(os.path.join(folder, f"*{ext}")))

    # Remove duplicates and sort for consistency
    files = sorted(set(files))

    # If suffix filters are provided, filter the file list
    if suffixes:
        filtered_files = []
        # Sort extensions by length to remove the longest first (e.g., '.bam.gz' before '.bam')
        sorted_exts = sorted(all_extensions, key=len, reverse=True)

        for file_path in files:
            base_name = file_path

            # Strip the extension (longest possible match first)
            for ext in sorted_exts:
                if base_name.endswith(ext):
                    base_name = base_name[:-len(ext)]
                    break

            # Include the file if it ends with any of the provided suffixes
            if any(base_name.endswith(suffix) for suffix in suffixes):
                filtered_files.append(file_path)

        files = filtered_files

    return files


# ----------------------------------------------------------------------
# FASTQ / SingleM input discovery
# ----------------------------------------------------------------------
DEFAULT_FASTQ_EXTENSIONS = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
DEFAULT_SINGLEM_EXTENSIONS = [".tsv", ".tsv.gz"]


def get_sample_name(path_or_name: str) -> str:
    """
    Derive the MGCalibrator sample name from a file name or a SingleM sample name.

    The sample name consists of the first two underscore-separated parts of the
    base name, ignoring everything after the first dot. For example,
    'sample_1_input_R1.fastq.gz' -> 'sample_1'.

    Parameters
    ----------
    path_or_name : str
        File path, file name or sample name.

    Returns
    -------
    str
        The derived sample name.
    """
    base_name = os.path.basename(path_or_name)
    name_without_ext = base_name.split(".")[0]
    return "_".join(name_without_ext.split("_")[0:2])


def list_files_by_extension(
    folder: str,
    extensions: List[str],
    suffixes: Optional[List[str]] = None,
) -> List[str]:
    """
    List files in a folder (non-recursive) that end with one of the given extensions.

    Parameters
    ----------
    folder : str
        Directory to search.
    extensions : list of str
        File extensions to include (e.g. ['.fastq.gz', '.fq.gz']).
    suffixes : list of str, optional
        If provided, only files whose name ends with one of these suffixes
        (after removing the extension) are returned, e.g. ['_1', '_2'] to
        skip singleton files.

    Returns
    -------
    list of str
        Sorted list of matching file paths.
    """
    if not os.path.isdir(folder):
        raise NotADirectoryError(f"Folder does not exist: {folder}")

    # Remove the longest matching extension first (e.g. '.fastq.gz' before '.gz')
    sorted_exts = sorted(extensions, key=len, reverse=True)

    files = []
    for name in sorted(os.listdir(folder)):
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            continue

        matching_ext = next((ext for ext in sorted_exts if name.endswith(ext)), None)
        if matching_ext is None:
            continue

        if suffixes:
            base_name = name[: -len(matching_ext)]
            if not any(base_name.endswith(suffix) for suffix in suffixes):
                continue

        files.append(path)

    return files


def group_files_by_sample(files: List[str]) -> Dict[str, List[str]]:
    """
    Group file paths by the sample name derived from their file names.

    Parameters
    ----------
    files : list of str
        File paths (e.g. FASTQ files).

    Returns
    -------
    dict
        Mapping of sample name -> sorted list of file paths.
    """
    grouped: Dict[str, List[str]] = {}
    for path in files:
        grouped.setdefault(get_sample_name(path), []).append(path)
    return {sample: sorted(paths) for sample, paths in grouped.items()}

DEFAULT_SYLPH_EXTENSIONS = [".sylphmpa", ".sylphmpa.gz"]
DEFAULT_SYLPH_PROFILE_EXTENSIONS = [".tsv", ".tsv.gz"]


def collect_files(paths: List[str], extensions: List[str]) -> List[str]:
    """
    Collect files from a list of folders and/or individual files.

    Folders are searched (non-recursively) for files with one of the given
    extensions. Individual files are always included, whatever their extension.
    Duplicates are removed while keeping the order.

    Parameters
    ----------
    paths : list of str
        Folders and/or file paths.
    extensions : list of str
        File extensions to search for inside folders.

    Returns
    -------
    list of str
        File paths.
    """
    files: List[str] = []
    for path in paths:
        if os.path.isdir(path):
            files.extend(list_files_by_extension(path, extensions))
        elif os.path.isfile(path):
            files.append(path)
        else:
            raise FileNotFoundError(f"Path does not exist: {path}")

    seen = set()
    return [f for f in files if not (f in seen or seen.add(f))]
