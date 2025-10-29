import os
from glob import glob
from typing import List, Optional


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
