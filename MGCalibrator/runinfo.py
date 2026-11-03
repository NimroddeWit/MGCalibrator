"""
runinfo.py

Writes a small JSON file next to the output with everything that is needed to reproduce a run:
the command line, all parameters, the software versions and a few numbers about the data.
"""

import json
import logging
import os
import platform
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from typing import Any, Dict, Optional

from . import __version__

PYTHON_PACKAGES = ["numpy", "pandas", "scipy", "pysam"]
EXTERNAL_TOOLS = ["samtools", "coverm"]  # only used in BAM mode


def run_info_path(output_csv: str) -> str:
    """Path of the run info file that belongs to an output file: '<output>_run_info.json'."""
    return os.path.splitext(output_csv)[0] + "_run_info.json"


def _package_version(name: str) -> Optional[str]:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _tool_version(tool: str) -> Optional[str]:
    """First line of '<tool> --version', or None if the tool cannot be run."""
    try:
        result = subprocess.run([tool, "--version"], capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    # Decode leniently: some tools print non-UTF-8 characters (e.g. a copyright sign)
    output = (result.stdout or result.stderr).decode("utf-8", errors="replace")
    lines = output.strip().splitlines()
    return lines[0].strip() if lines else None


def write_run_info(
    output_csv: str,
    args,
    mode_name: str,
    extra: Optional[Dict[str, Any]] = None,
    argv: Optional[list] = None,
) -> str:
    """
    Write the run info file next to `output_csv` and return its path.

    Parameters
    ----------
    output_csv : str
        Path of the main output file.
    args : argparse.Namespace
        All parsed command line parameters (saved as they are, including the defaults).
    mode_name : str
        Input mode of the run (BAM, SingleM + FASTQ, Sylph + FASTQ).
    extra : dict, optional
        Additional information about the data, such as the number of samples and rows.
    argv : list, optional
        Command line (default: sys.argv).
    """
    info: Dict[str, Any] = {
        "mgcalibrator_version": __version__,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": " ".join(shlex.quote(a) for a in (sys.argv if argv is None else argv)),
        "input_mode": mode_name,
        "seed": getattr(args, "seed", None),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {name: _package_version(name) for name in PYTHON_PACKAGES},
        "parameters": vars(args),
    }
    if mode_name == "BAM":
        info["external_tools"] = {tool: _tool_version(tool) for tool in EXTERNAL_TOOLS}
    info.update(extra or {})

    path = run_info_path(output_csv)
    with open(path, "w") as handle:
        json.dump(info, handle, indent=2, default=str)
        handle.write("\n")
    logging.info(f"Run information saved to: {path}")
    return path
