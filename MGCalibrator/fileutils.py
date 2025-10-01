import os
from glob import glob

def list_bam_files(folder, extensions=None, suffixes=None, singleton_files=None):
    if extensions is None:
        extensions = [".sorted.bam", ".sort.bam"]

    files = []
    all_extensions = []
    for ext in extensions:
        all_extensions.append(ext)
        all_extensions.append(f"{ext}.gz")
        all_extensions.append(f"{ext}.gzip")

    # Find all files with the specified extensions
    for ext in all_extensions:
        files.extend(glob(os.path.join(folder, f"*{ext}")))
    
    # Remove duplicates that might arise from overlapping patterns
    files = sorted(list(set(files)))

    # Filter by suffixes if provided
    if suffixes:
        filtered_files = []
        for f in files:
            # Check if the file name, stripped of its extension, ends with a given suffix
            fn_no_ext = f
            # Sort extensions by length to remove longest match first (e.g., .fastq.gz before .fastq)
            for ext_to_strip in sorted(all_extensions, key=len, reverse=True):
                if fn_no_ext.endswith(ext_to_strip):
                    fn_no_ext = fn_no_ext[:-len(ext_to_strip)]
                    break
            
            if any(fn_no_ext.endswith(s) for s in suffixes):
                filtered_files.append(f)
        files = filtered_files
        
    return files
