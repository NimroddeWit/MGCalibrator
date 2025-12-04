"""
MGCalibrator Command-Line Interface (CLI)

This module provides the command-line interface for MGCalibrator, a tool for
calculating absolute microbial abundances with associated uncertainty from
BAM files. It handles file discovery, filtering, depth computation,
Monte Carlo error estimation, and calibration.
"""

import argparse
import logging
import os
import subprocess

import numpy as np
import pandas as pd

from .fileutils import list_bam_files
from .processor import (
    run_coverm_filter,
    compute_depths_with_error,
    calculate_scaling_factors,
)

# ----------------------------------------------------------------------
# Startup Banner Function
# ----------------------------------------------------------------------
def print_startup_message():
    """Display a styled startup banner for MGCalibrator."""
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    RESET = "\033[0m"

    print("\n" + GREEN + "=" * 104 + RESET)
    print(CYAN + r'  888b     d888  .d8888b.   .d8888b.           888 d8b 888                      888                    ' + RESET)
    print(CYAN + r'  8888b   d8888 d88P  Y88b d88P  Y88b          888 Y8P 888                      888                    ' + RESET)
    print(CYAN + r'  88888b.d88888 888    888 888    888          888     888                      888                    ' + RESET)
    print(CYAN + r'  888Y88888P888 888        888         8888b.  888 888 88888b.  888d888 8888b.  888888 .d88b.  888d888 ' + RESET)
    print(CYAN + r'  888 Y888P 888 888  88888 888            "88b 888 888 888 "88b 888P"      "88b 888   d88""88b 888P"   ' + RESET)
    print(CYAN + r'  888  Y8P  888 888    888 888    888 .d888888 888 888 888  888 888    .d888888 888   888  888 888     ' + RESET)
    print(CYAN + r'  888   "   888 Y88b  d88P Y88b  d88P 888  888 888 888 888 d88P 888    888  888 Y88b. Y88..88P 888     ' + RESET)
    print(CYAN + r'  888       888  "Y8888P88  "Y8888P"  "Y888888 888 888 88888P"  888    "Y888888  "Y888 "Y88P"  888     ' + RESET)
    print(GREEN + "=" * 104 + RESET + "\n")
    print(CYAN + "                          MGCalibrator – Metagenomic abundance calibration tool" + RESET + "\n")

# -------------------------------------------------------------------------
# Main entry point
# -------------------------------------------------------------------------
def main() -> None:
    """Main function for the MGCalibrator command-line interface."""

    # ----------------------
    # Argument parser setup
    # ----------------------
    parser = argparse.ArgumentParser(
    prog="MGCalibrator",
    formatter_class=argparse.RawTextHelpFormatter,
    description=(
        "MGCalibrator — Calculate absolute microbial abundances with standard error "
        "from metagenomic BAM files.\n\n"
        "Example usage:\n"
        "  mgcalibrator \\\n"
        "    --bam_folder data/bams \\\n"
        "    --dna_mass data/dna_mass.csv \\\n"
        "    --scaling_factors results/scaling_factors.txt \\\n"
        "    --output results/calibrated_depths.csv \\\n"
        "    --threads 8"
        ),
    )
    
    # --- Input/output options ---
    parser.add_argument(
        "--bam_folder",
        required=True,
        help="Folder containing BAM files.",
    )
    parser.add_argument(
        "--extensions",
        nargs="*",
        default=[".sorted.bam", ".sort.bam"],
        help="File extensions to search for (default: .sorted.bam, .sort.bam).",
    )
    parser.add_argument(
        "--suffixes",
        nargs="*",
        help="Optional list of suffixes to filter BAM files.",
    )
    parser.add_argument(
        "--reference_clusters",
        help="CSV file mapping reference sequences to clusters (e.g., stack multiple alleles into one gene).",
    )
    parser.add_argument(
        "--reference_bins",
        help="CSV file mapping reference sequences to bins (e.g., concatenate contigs from the same MAG).",
    )
    parser.add_argument(
        "--dna_mass",
        required=True,
        help="CSV file containing measured DNA mass (in ng).",
    )
    parser.add_argument(
        "--scaling_factors",
        required=True,
        help="Path to a file storing precomputed scaling factors (if it exists, it will be reused).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output CSV file.",
    )

    # --- Performance options ---
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of parallel threads to use (default: all available cores).",
    )

    # --- Filtering options ---
    parser.add_argument(
        "--skip_filtering",
        action="store_true",
        help="Skip filtering based on mapping identity.",
    )
    parser.add_argument(
        "--filter_percent",
        type=float,
        default=97.0,
        help="Minimum read percent identity for filtering (default: 97).",
    )

    # --- Calculation options ---
    parser.add_argument(
        "--depth_percent",
        type=float,
        default=50.0,
        help="Percentage of middle depth values to include for depth calculation (default: 50).",
    )
    parser.add_argument(
        "--qubit_error",
        type=float,
        default=0.0,
        help="Qubit error (default: 0.0).",
    )

    # --- Monte Carlo error propagation ---
    parser.add_argument(
        "--mc_samples",
        type=int,
        default=100,
        help="Number of Monte Carlo samples for confidence intervals (default: 100).",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=500,
        help="Number of simulations per batch (default: 500).",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=1.0,
        help="Pseudocount for read-mapping probability at zero-depth positions (default: 1.0).",
    )
    parser.add_argument(
        "--skip_error_calculation",
        action="store_true",
        help="Skip error calculation by Monte Carlo simulation.",
    )

    args = parser.parse_args()

    # ----------------------
    # Logging configuration
    # ----------------------
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    print_startup_message()
    logging.info("Starting MGCalibrator...")
    logging.info(f"Using up to {args.threads or 'all available'} threads.")

    # ----------------------
    # Discover BAM files
    # ----------------------
    bam_files = list_bam_files(
        folder=args.bam_folder,
        extensions=args.extensions,
        suffixes=args.suffixes,
    )
    if not bam_files:
        raise FileNotFoundError(
            f"No BAM files found in {args.bam_folder} with the specified extensions."
        )
    logging.info(f"Found {len(bam_files)} BAM files to process.")

    # ----------------------
    # Load metadata
    # ----------------------
    dna_mass_df = pd.read_csv(args.dna_mass)
    if "sample_id" not in dna_mass_df.columns or "DNA_mass" not in dna_mass_df.columns:
        raise ValueError("DNA mass file must contain 'sample_id' and 'DNA_mass' columns.")

    initial_dna_mass = dict(zip(dna_mass_df.sample_id, dna_mass_df.DNA_mass))

    reference_clusters_csv = args.reference_clusters
    reference_bins_csv = args.reference_bins
    min_read_perc_identity = args.filter_percent
    depth_percent = args.depth_percent
    qubit_error = args.qubit_error
    scaling_factors_file = args.scaling_factors
    output_dir = os.path.dirname(args.output)

    # ----------------------
    # Filter BAM files (optional)
    # ----------------------
    if not args.skip_filtering:
        logging.info(f"Filtering BAM files at ≥{min_read_perc_identity}% identity...")
        bam_files_filtered = run_coverm_filter(
            bam_files=bam_files,
            min_read_perc_identity=min_read_perc_identity,
            output_dir=output_dir,
        )
        logging.info("BAM filtering completed.")
    else:
        bam_files_filtered = bam_files
        # Ensure BAM index files exist
        for bam_file in bam_files_filtered:
            bai_file = f"{bam_file}.bai"
            if not os.path.exists(bai_file):
                try:
                    subprocess.run(["samtools", "index", bam_file], check=True)
                except subprocess.CalledProcessError as e:
                    logging.warning(f"Failed to index {bam_file}: {e}")

    # ----------------------
    # Compute raw depths with uncertainty
    # ----------------------
    logging.info("Computing raw depths with Monte Carlo error estimation...")
    calculated_depths = compute_depths_with_error(
        bam_files_filtered,
        reference_clusters_csv=reference_clusters_csv,
        reference_bins_csv=reference_bins_csv,
        depth_percent=depth_percent,
        skip_error_calculation=args.skip_error_calculation,
        n_simulations=args.mc_samples,
        pseudocount=args.pseudocount,
        batch_size=args.batch_size,
        n_jobs=args.threads,
    )

    # ----------------------
    # Compute scaling factors
    # ----------------------
    samples = set(calculated_depths["sample"])
    scaling_factors = calculate_scaling_factors(
        samples=samples,
        bam_files=bam_files,
        initial_dna_mass=initial_dna_mass,
        scaling_factors_file=scaling_factors_file,
        n_workers=args.threads,
    )

    # ----------------------
    # Apply Qubit error and calibrate depth estimates
    # ----------------------
    logging.info("Applying Qubit error and calibrating depth estimates...")

    # Apply Qubit error and calibration
    calibrated_depths = (
        calculated_depths
        .assign(
            # Asymmetric relative errors (before Qubit correction)
            relative_error_down=lambda df: (df["depth"] - df["lower_ci"]) / df["depth"],
            relative_error_up=lambda df: (df["upper_ci"] - df["depth"]) / df["depth"],

            # Raw absolute errors
            error_down=lambda df: df["depth"] * df["relative_error_down"],
            error_up=lambda df: df["depth"] * df["relative_error_up"],

            # Combine uncertainty with Qubit measurement error
            total_relative_error_down=lambda df: np.sqrt(df["relative_error_down"]**2 + qubit_error**2),
            total_relative_error_up=lambda df: np.sqrt(df["relative_error_up"]**2 + qubit_error**2),

            # Map scaling factors and apply calibration
            scaling_factor=lambda df: df["sample"].map(scaling_factors),
            calibrated_depth=lambda df: df["depth"] * df["scaling_factor"],

            # Propagate errors through calibration
            calibrated_error_down=lambda df: df["calibrated_depth"] * df["total_relative_error_down"],
            calibrated_error_up=lambda df: df["calibrated_depth"] * df["total_relative_error_up"],
        )
        # Keep only relevant columns for output
        .loc[
            :,
            [
                "sample",
                "reference",
                "depth",
                "error_down",
                "error_up",
                "calibrated_depth",
                "calibrated_error_down",
                "calibrated_error_up",
            ],
        ]
    )

    logging.info("Qubit error propagation and calibration complete.")

    # ----------------------
    # Save output
    # ----------------------
    os.makedirs(output_dir, exist_ok=True)
    calibrated_depths.to_csv(args.output, index=False)
    logging.info(f"Results saved to: {args.output}")

    # ----------------------
    # Final message
    # ----------------------
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    RESET = "\033[0m"

    logging.info(GREEN + "=" * 49 + RESET)
    logging.info(CYAN + "  🌟 MGCalibrator has finished successfully! 🌟" + RESET)
    logging.info(GREEN + "=" * 49 + RESET)
    logging.info(CYAN + "       Thank you for using MGCalibrator!\n" + RESET)
