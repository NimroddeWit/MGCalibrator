"""
MGCalibrator: A command-line tool for calculating absolute abundances from
metagenomic data (BAM files, SingleM output + FASTQ files, or Sylph output + FASTQ files).

This package is intended to be used exclusively via the command line.
The main entry point is defined in `cli.py`.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("MGCalibrator")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "unknown"
