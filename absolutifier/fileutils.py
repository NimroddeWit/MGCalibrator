import os
from glob import glob

def list_fastq_files(folder, extension=".fastq", suffixes=None, singleton_files=None):
    files = glob(os.path.join(folder, f"*{extension}"))
    if suffixes:
        files = [f for f in files if any(f.endswith(s + extension) for s in suffixes)]
    if singleton_files:
        files += singleton_files
    return files
