"""
MGCalibrator Command-Line Interface (CLI)

This module provides the command-line interface for MGCalibrator, a tool for
calculating absolute microbial abundances. It supports three (mutually exclusive)
input modes:

  * BAM mode: BAM files as input. Handles file discovery, filtering, depth
    computation, Monte Carlo error estimation, and calibration.
  * SingleM mode: SingleM output files + the raw FASTQ files as input. The
    SingleM species-level coverage values are used as depths, so no filtering
    or Monte Carlo simulation is needed.
  * Sylph mode: Sylph results (.sylphmpa files, optionally with the Sylph profile
    TSV files for confidence intervals) + the raw FASTQ files as input. The
    strain-level Sylph coverage values are used as depths.

All modes produce an output file with the same columns.
"""

import argparse
import logging
import os
import subprocess
from typing import Tuple

import numpy as np
import pandas as pd

from . import __version__
from .fileutils import (
    DEFAULT_FASTQ_EXTENSIONS,
    DEFAULT_SINGLEM_EXTENSIONS,
    DEFAULT_SYLPH_EXTENSIONS,
    DEFAULT_SYLPH_PROFILE_EXTENSIONS,
    collect_files,
    list_bam_files,
    list_files_by_extension,
)
from .processor import (
    filtered_bam_path,
    run_coverm_filter,
    compute_depths_with_error,
    calculate_scaling_factors,
)
from .runinfo import write_run_info
from .samples import (
    assign_files_to_samples,
    build_sample_map,
    check_sample_ids,
    group_by_sample,
    read_sample_map,
    unmatched_files_message,
)
from .singlem import read_singlem_tables, singlem_to_depth_table
from .sylph import (
    add_poisson_confidence_intervals,
    add_read_counts,
    read_sylph_profiles,
    read_sylphmpa_files,
    sylph_to_depth_table,
)

def apply_min_depth(depths: pd.DataFrame, min_depth: float) -> Tuple[pd.DataFrame, int]:
    """
    Remove the rows with a depth below `min_depth` (0 keeps everything).

    Returns the remaining rows and the number of removed rows. Raises a ValueError if no row is left.
    """
    if min_depth <= 0:
        return depths, 0
    keep = ~(depths["depth"] < min_depth)
    n_removed = int((~keep).sum())
    if n_removed:
        logging.info(f"Removed {n_removed} of {len(depths)} rows with a depth below --min_depth {min_depth}.")
    if not keep.any():
        raise ValueError(f"No rows left: all depths are below --min_depth {min_depth}.")
    return depths[keep].reset_index(drop=True), n_removed


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
        "MGCalibrator — Calculate absolute microbial abundances from metagenomic data.\n\n"
        "Three input modes are supported (use exactly one of them):\n"
        "  1. BAM mode:      --bam_folder\n"
        "  2. SingleM mode:  --singlem_folder together with --fastq_folder\n"
        "  3. Sylph mode:    --sylph_folder together with --fastq_folder\n"
        "                    (optionally --sylph_profile_folder for confidence intervals;\n"
        "                    run Sylph with --estimate-unknown, and with --estimate-read-counts\n"
        "                    if you want confidence intervals)\n\n"
        "Example usage (BAM mode):\n"
        "  mgcalibrator \\\n"
        "    --bam_folder data/bams \\\n"
        "    --dna_mass data/dna_mass.csv \\\n"
        "    --scaling_factors results/scaling_factors.txt \\\n"
        "    --output results/calibrated_depths.csv \\\n"
        "    --threads 8\n\n"
        "Example usage (SingleM mode):\n"
        "  mgcalibrator \\\n"
        "    --singlem_folder data/singlem \\\n"
        "    --fastq_folder data/fastq \\\n"
        "    --dna_mass data/dna_mass.csv \\\n"
        "    --scaling_factors results/scaling_factors.txt \\\n"
        "    --output results/calibrated_depths.csv\n\n"
        "Example usage (Sylph mode, two databases, with confidence intervals):\n"
        "  mgcalibrator \\\n"
        "    --sylph_folder results/sylph_gtdb results/sylph_virus \\\n"
        "    --sylph_profile_folder results/sylph_gtdb_profiles results/sylph_virus_profiles \\\n"
        "    --fastq_folder data/fastq \\\n"
        "    --dna_mass data/dna_mass.csv \\\n"
        "    --scaling_factors results/scaling_factors.txt \\\n"
        "    --output results/calibrated_depths.csv"
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # --- Input/output options ---
    parser.add_argument(
        "--bam_folder",
        help="BAM mode: folder containing BAM files.",
    )
    parser.add_argument(
        "--singlem_folder",
        help=(
            "SingleM mode: folder containing the SingleM output files (standard SingleM taxonomic\n"
            "profile: TSV with the columns sample, coverage and taxonomy; only species-level rows\n"
            "are used). Requires --fastq_folder."
        ),
    )
    parser.add_argument(
        "--sylph_folder",
        nargs="+",
        help=(
            "Sylph mode: one or more folders (or files) with .sylphmpa files, as written by\n"
            "`sylph-tax taxprof` (one file per sample). Results of separate Sylph runs, e.g. with\n"
            "different databases, can be given together; files of the same sample are combined\n"
            "using the sample name in the file header. Only strain-level rows are used.\n"
            "Sylph must have been run with --estimate-unknown (-u), so that the coverage is the\n"
            "true coverage (True_cov). Requires --fastq_folder."
        ),
    )
    parser.add_argument(
        "--sylph_profile_folder",
        nargs="+",
        help=(
            "Sylph mode (optional): one or more folders (or files) with the `sylph profile` output\n"
            "TSV files that the .sylphmpa files were made from (single- or multi-sample TSVs).\n"
            "Used to calculate confidence intervals around the coverage (error_down and error_up):\n"
            "the exact Poisson 95% interval of the estimated number of reads per genome (column\n"
            "Sequence_abundance) is applied to the coverage. Requires `sylph profile\n"
            "--estimate-read-counts` (Sylph >= 0.9.0).\n"
            "The confidence intervals that Sylph reports itself are not used."
        ),
    )
    parser.add_argument(
        "--fastq_folder",
        help=(
            "SingleM/Sylph mode: folder containing the raw FASTQ files (optionally gzip-compressed),\n"
            "1 file per sample for Nanopore reads or 2 files per sample for paired-end\n"
            "Illumina reads. Requires --singlem_folder or --sylph_folder."
        ),
    )
    parser.add_argument(
        "--extensions",
        nargs="*",
        default=[".sorted.bam", ".sort.bam"],
        help="BAM mode: file extensions to search for (default: .sorted.bam, .sort.bam).",
    )
    parser.add_argument(
        "--fastq_extensions",
        nargs="*",
        default=DEFAULT_FASTQ_EXTENSIONS,
        help="SingleM mode: FASTQ file extensions to search for (default: .fastq.gz .fq.gz .fastq .fq).",
    )
    parser.add_argument(
        "--singlem_extensions",
        nargs="*",
        default=DEFAULT_SINGLEM_EXTENSIONS,
        help="SingleM mode: SingleM output file extensions to search for (default: .tsv .tsv.gz).",
    )
    parser.add_argument(
        "--sylphmpa_extensions",
        nargs="*",
        default=DEFAULT_SYLPH_EXTENSIONS,
        help="Sylph mode: .sylphmpa file extensions to search for in folders (default: .sylphmpa .sylphmpa.gz).",
    )
    parser.add_argument(
        "--sylph_profile_extensions",
        nargs="*",
        default=DEFAULT_SYLPH_PROFILE_EXTENSIONS,
        help="Sylph mode: Sylph profile file extensions to search for in folders (default: .tsv .tsv.gz).",
    )
    parser.add_argument(
        "--suffixes",
        nargs="*",
        help=(
            "Optional list of suffixes to filter BAM files (BAM mode) or FASTQ files (SingleM/Sylph mode),\n"
            "matched at the end of the file name before the extension (e.g. '_1 _2' to skip\n"
            "singleton FASTQ files)."
        ),
    )
    parser.add_argument(
        "--sample_map",
        help=(
            "Optional CSV file that links input names to the sample_id of the DNA mass file, with the\n"
            "columns 'name' (or 'singlem_sample' / 'sylph_sample') and 'sample_id'. A name can be a\n"
            "BAM/FASTQ file name or a sample name in the SingleM/Sylph output. Only needed for names\n"
            "that do not start with their sample_id (followed by '_', '.' or '-')."
        ),
    )
    parser.add_argument(
        "--cluster",
        help="CSV file mapping reference sequences to clusters (e.g., stack multiple alleles into one gene).",
    )
    parser.add_argument(
        "--bin",
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
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help=(
            "BAM mode: random seed to make the Monte Carlo simulations reproducible (default: not set, so\n"
            "the intervals differ slightly between runs). The seed is saved in the run info file."
        ),
    )
    parser.add_argument(
        "--min_depth",
        type=float,
        default=0.0,
        help=(
            "Remove references/taxa with a depth below this value before the calibration, e.g. to leave out\n"
            "noise at very low coverage (default: 0, keep everything)."
        ),
    )
    parser.add_argument(
        "--overwrite_scaling_factors",
        action="store_true",
        help=(
            "Ignore an existing --scaling_factors file and calculate all scaling factors again. Without this\n"
            "option, the factors in an existing file are reused (the sequence files are not counted again)."
        ),
    )

    args = parser.parse_args()
    if args.min_depth < 0:
        parser.error("--min_depth cannot be negative.")

    # ----------------------
    # Determine input mode
    # ----------------------
    bam_mode = args.bam_folder is not None
    singlem_mode = args.singlem_folder is not None
    sylph_mode = args.sylph_folder is not None
    fastq_mode = singlem_mode or sylph_mode

    selected_modes = [name for name, active in
                      (("--bam_folder", bam_mode), ("--singlem_folder", singlem_mode), ("--sylph_folder", sylph_mode))
                      if active]
    if len(selected_modes) > 1:
        parser.error(
            f"Provide input for only one mode, but got {' and '.join(selected_modes)}. Use either a BAM "
            "folder (--bam_folder), SingleM output (--singlem_folder) or Sylph output (--sylph_folder), "
            "with --fastq_folder for the latter two."
        )
    if not selected_modes:
        parser.error(
            "No input provided. Use --bam_folder (BAM mode), --singlem_folder with --fastq_folder "
            "(SingleM mode) or --sylph_folder with --fastq_folder (Sylph mode)."
        )
    if fastq_mode and args.fastq_folder is None:
        parser.error(f"{selected_modes[0]} requires --fastq_folder.")
    if bam_mode and args.fastq_folder is not None:
        parser.error("--fastq_folder can only be used together with --singlem_folder or --sylph_folder.")
    if args.sylph_profile_folder is not None and not sylph_mode:
        parser.error("--sylph_profile_folder can only be used together with --sylph_folder.")

    mode_name = "BAM" if bam_mode else ("SingleM + FASTQ" if singlem_mode else "Sylph + FASTQ")

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
    logging.info(f"Input mode: {mode_name}.")
    logging.info(f"Using up to {args.threads or 'all available'} threads.")

    # Warn about options that have no effect in the selected mode
    calculation_options = [
        "extensions", "cluster", "bin", "skip_filtering", "filter_percent", "depth_percent",
        "mc_samples", "batch_size", "pseudocount", "skip_error_calculation", "seed",
    ]
    if bam_mode:
        mode_specific_options = [
            "fastq_extensions", "singlem_extensions", "sylphmpa_extensions",
            "sylph_profile_extensions",
        ]
    elif singlem_mode:
        # SingleM provides no confidence intervals, so the Qubit error cannot be propagated
        mode_specific_options = calculation_options + ["qubit_error", "sylphmpa_extensions", "sylph_profile_extensions"]
    else:
        mode_specific_options = calculation_options + ["singlem_extensions"]
        if args.sylph_profile_folder is None:
            mode_specific_options.append("qubit_error")
    ignored_options = [
        f"--{name}" for name in mode_specific_options
        if getattr(args, name) != parser.get_default(name)
    ]
    if ignored_options:
        logging.warning(
            f"Ignoring option(s) that are not used in {mode_name} mode: "
            f"{', '.join(ignored_options)}"
        )

    # ----------------------
    # Load metadata
    # ----------------------
    dna_mass_df = pd.read_csv(args.dna_mass, dtype={"sample_id": str})
    if "sample_id" not in dna_mass_df.columns or "DNA_mass" not in dna_mass_df.columns:
        raise ValueError("DNA mass file must contain 'sample_id' and 'DNA_mass' columns.")

    dna_mass_df["sample_id"] = dna_mass_df["sample_id"].str.strip()
    # The sample_id column defines the sample names; all input files are matched against it.
    known_sample_ids = check_sample_ids(dna_mass_df["sample_id"])
    initial_dna_mass = dict(zip(dna_mass_df.sample_id, dna_mass_df.DNA_mass))
    manual_map = read_sample_map(args.sample_map) if args.sample_map else {}

    reference_clusters_csv = args.cluster
    reference_bins_csv = args.bin
    min_read_perc_identity = args.filter_percent
    depth_percent = args.depth_percent
    qubit_error = args.qubit_error
    scaling_factors_file = args.scaling_factors
    output_dir = os.path.dirname(args.output)

    if bam_mode:
        # ==============================================================
        # BAM mode
        # ==============================================================
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
        # Assign BAM files to samples (sample_id of the DNA mass file)
        # ----------------------
        sample_by_bam, unmatched_bams = assign_files_to_samples(bam_files, known_sample_ids, manual_map)
        if unmatched_bams:
            raise ValueError(unmatched_files_message(unmatched_bams, known_sample_ids, "BAM"))
        bams_by_sample = group_by_sample(sample_by_bam)
        multiple_bams = {sample: paths for sample, paths in bams_by_sample.items() if len(paths) > 1}
        if multiple_bams:
            details = "; ".join(
                f"{sample}: {', '.join(os.path.basename(p) for p in paths)}"
                for sample, paths in sorted(multiple_bams.items())
            )
            raise ValueError(
                "BAM mode expects exactly one BAM file per sample, but several BAM files belong to the "
                f"same sample ({details}). Merge the BAM files per sample, use --suffixes to select one, "
                "or give the samples different sample_ids."
            )

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
            # The filtered BAM files inherit the sample of the BAM file they were made from
            sample_by_filtered_bam = {
                filtered_bam_path(bam, output_dir, min_read_perc_identity): sample
                for bam, sample in sample_by_bam.items()
            }
            sample_by_filtered_bam = {
                bam: sample_by_filtered_bam[bam] for bam in bam_files_filtered
            }
        else:
            bam_files_filtered = bam_files
            sample_by_filtered_bam = sample_by_bam
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
            sample_by_bam=sample_by_filtered_bam,
            reference_clusters_csv=reference_clusters_csv,
            reference_bins_csv=reference_bins_csv,
            depth_percent=depth_percent,
            skip_error_calculation=args.skip_error_calculation,
            n_simulations=args.mc_samples,
            pseudocount=args.pseudocount,
            batch_size=args.batch_size,
            n_jobs=args.threads,
            seed=args.seed,
        )
        calculated_depths, n_removed_min_depth = apply_min_depth(calculated_depths, args.min_depth)

        # ----------------------
        # Compute scaling factors
        # ----------------------
        samples = set(calculated_depths["sample"])
        scaling_factors, reused_scaling_factors = calculate_scaling_factors(
            samples=samples,
            sample_by_bam=sample_by_bam,
            initial_dna_mass=initial_dna_mass,
            scaling_factors_file=scaling_factors_file,
            n_workers=args.threads,
            overwrite=args.overwrite_scaling_factors,
        )

    else:
        # ==============================================================
        # SingleM / Sylph mode (profiler output + raw FASTQ files)
        # ==============================================================
        # ----------------------
        # Discover FASTQ files and group them per sample
        # ----------------------
        fastq_files = list_files_by_extension(
            folder=args.fastq_folder,
            extensions=args.fastq_extensions,
            suffixes=args.suffixes,
        )
        if not fastq_files:
            raise FileNotFoundError(
                f"No FASTQ files found in {args.fastq_folder} with the specified extensions."
            )
        sample_by_fastq, unmatched_fastq = assign_files_to_samples(fastq_files, known_sample_ids, manual_map)
        if unmatched_fastq:
            logging.warning(
                "Ignoring FASTQ file(s) that do not belong to a sample_id of the DNA mass file: "
                f"{', '.join(os.path.basename(p) for p in unmatched_fastq)}."
            )
        if not sample_by_fastq:
            raise ValueError(unmatched_files_message(unmatched_fastq, known_sample_ids, "FASTQ"))
        fastq_files_by_sample = group_by_sample(sample_by_fastq)
        logging.info(
            f"Found {len(sample_by_fastq)} FASTQ files for {len(fastq_files_by_sample)} samples."
        )
        for sample, paths in fastq_files_by_sample.items():
            if len(paths) > 2:
                logging.warning(
                    f"Sample '{sample}' has {len(paths)} FASTQ files (expected 1 for Nanopore or 2 for "
                    "paired-end Illumina). The bases of all files are counted; use --suffixes to "
                    "select the right files if this is not intended."
                )

        if singlem_mode:
            # ----------------------
            # Read SingleM output (species level)
            # ----------------------
            singlem_files = list_files_by_extension(
                folder=args.singlem_folder,
                extensions=args.singlem_extensions,
            )
            if not singlem_files:
                raise FileNotFoundError(
                    f"No SingleM output files found in {args.singlem_folder} with the specified extensions."
                )
            logging.info(f"Found {len(singlem_files)} SingleM output files.")
            singlem_df = read_singlem_tables(singlem_files)

            sample_map = build_sample_map(
                profiler_samples=singlem_df["sample"].unique(),
                known_ids=known_sample_ids,
                manual_map=manual_map,
                source="SingleM",
            )

            # SingleM coverage -> depth, taxonomy -> reference. SingleM provides no
            # confidence interval, so the error columns are left as NaN (written as NA).
            calculated_depths = singlem_to_depth_table(singlem_df, sample_map).assign(
                lower_ci=np.nan,
                upper_ci=np.nan,
                read_counts=np.nan,
                error_note="no_error_estimate_for_singlem",
            )

        else:
            # ----------------------
            # Read Sylph results (strain level)
            # ----------------------
            sylphmpa_files = collect_files(args.sylph_folder, args.sylphmpa_extensions)
            if not sylphmpa_files:
                raise FileNotFoundError(
                    f"No .sylphmpa files found in {' '.join(args.sylph_folder)} with the specified extensions."
                )
            logging.info(f"Found {len(sylphmpa_files)} .sylphmpa files.")
            strains = read_sylphmpa_files(sylphmpa_files)

            logging.info(
                "Sylph mode: the coverage in the .sylphmpa files is assumed to be the true coverage "
                "(Sylph run with --estimate-unknown)."
            )
            if args.sylph_profile_folder is not None:
                profile_files = collect_files(args.sylph_profile_folder, args.sylph_profile_extensions)
                if not profile_files:
                    raise FileNotFoundError(
                        "No Sylph profile files found in "
                        f"{' '.join(args.sylph_profile_folder)} with the specified extensions."
                    )
                logging.info(f"Found {len(profile_files)} Sylph profile files.")
                # Confidence intervals: exact Poisson interval of the estimated read count of each genome
                strains = add_read_counts(strains, read_sylph_profiles(profile_files))
                strains = add_poisson_confidence_intervals(strains)

            sample_map = build_sample_map(
                profiler_samples=strains["sample"].unique(),
                known_ids=known_sample_ids,
                manual_map=manual_map,
                source="Sylph",
            )
            calculated_depths = sylph_to_depth_table(strains, sample_map)

        calculated_depths, n_removed_min_depth = apply_min_depth(calculated_depths, args.min_depth)
        samples = set(calculated_depths["sample"])
        unused_fastq_samples = sorted(set(fastq_files_by_sample) - samples)
        if unused_fastq_samples:
            logging.warning(
                f"No {'SingleM' if singlem_mode else 'Sylph'} results found for FASTQ sample(s): "
                f"{', '.join(unused_fastq_samples)}. They are not included in the output."
            )

        # ----------------------
        # Compute scaling factors (bases counted in the FASTQ files)
        # ----------------------
        scaling_factors, reused_scaling_factors = calculate_scaling_factors(
            samples=samples,
            sample_by_bam=None,
            initial_dna_mass=initial_dna_mass,
            scaling_factors_file=scaling_factors_file,
            n_workers=args.threads,
            fastq_files_by_sample=fastq_files_by_sample,
            overwrite=args.overwrite_scaling_factors,
        )

    unused_ids = sorted(set(known_sample_ids) - set(calculated_depths["sample"]))
    if unused_ids:
        logging.warning(
            "No results found for sample_id(s) in the DNA mass file: "
            f"{', '.join(unused_ids)}. They are not included in the output."
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
                "scaling_factor",
                "read_counts",
                "error_note",
            ],
        ]
    )

    logging.info("Qubit error propagation and calibration complete.")

    # ----------------------
    # Save output
    # ----------------------
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    # In SingleM/Sylph mode, missing errors (no confidence interval available) are written as NA.
    calibrated_depths.to_csv(args.output, index=False, na_rep="NA" if fastq_mode else "")
    logging.info(f"Results saved to: {args.output}")

    write_run_info(
        args.output,
        args,
        mode_name,
        extra={
            "data": {
                "n_samples": int(calibrated_depths["sample"].nunique()),
                "n_output_rows": int(len(calibrated_depths)),
                "n_rows_removed_by_min_depth": n_removed_min_depth,
                "n_rows_without_error_estimate": int(calibrated_depths["error_down"].isna().sum()),
            },
            "scaling_factors": {
                "file": scaling_factors_file,
                "overwritten": bool(args.overwrite_scaling_factors),
                "reused_from_file": reused_scaling_factors,
                "calculated_in_this_run": sorted(set(calibrated_depths["sample"]) - set(reused_scaling_factors)),
            },
        },
    )

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