import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import os
import subprocess
from .processor import run_coverm_filter, compute_raw_depths_with_error, calculate_scaling_factors
from .fileutils import list_bam_files

def main():
    parser = argparse.ArgumentParser(description="MGCalibrator - calculate absolute abundances with standard error")
    
    # Input/output options
    parser.add_argument("--bam_folder", required=True, help="Folder containing BAM files")
    parser.add_argument("--extensions", nargs='*', default=[".sorted.bam", ".sort.bam"], 
                        help="List of sequence file extensions to search for (default: .sorted.bam, .sort.bam)")
    parser.add_argument("--suffixes", nargs="*", help="List of suffixes to filter in files")
    parser.add_argument("--reference_clusters", help="CSV file indicating to which cluster belongs each reference sequence")
    parser.add_argument("--reference_bins", help="CSV file indicating to which bin belongs each reference sequence")
    parser.add_argument("--dna_mass", required=True, help="CSV file with measured DNA mass (ng)") 
    parser.add_argument("--scaling_factors", required=True, help="Dictionary containing calculated scaling factors (loaded when file exists)")
    parser.add_argument("--output", required=True, help="Output CSV file")
    
    # Performance options
    parser.add_argument("--threads", type=int, default=None,
                        help="Number of threads to use for parallel processing (default: all available cores)")
    
    # Depth calculation options
    parser.add_argument("--skip_filtering", action='store_true', help="Skip filtering based on mapping identity")
    parser.add_argument("--perc_ident", default=97, help="Minimal mapping identity (%) for calculating depth")
    
    # Error propagation options
    parser.add_argument("--mc_samples", type=int, default=100, 
                       help="Number of Monte Carlo samples for confidence intervals (default: 100)")
    parser.add_argument("--batch_size", type=int, default=500, 
                       help="Number of simulations to run per batch (default: 500)")
    parser.add_argument("--pseudocount", type=float, default=1.0,
                       help="Pseudocount for read-mapping probability of zero-depth positions during MC simulation (default: 1.0)")
    
    # Plotting options
    parser.add_argument("--plot", action="store_true",
                       help="Generate visualization plots of features across samples with confidence intervals")
    parser.add_argument("--top_features", type=int, default=20,
                       help="Number of top abundant features to show in plots (default: 20)")
    parser.add_argument("--plot_format", choices=["png", "pdf", "svg"], default="png",
                       help="Output format for plots (default: png)")
    parser.add_argument("--figsize", nargs=2, type=float, default=[12, 8],
                       help="Figure size as width height (default: 12 8)")
    
    args = parser.parse_args()

    # --- Setup Logging ---
    logging.basicConfig(level=logging.DEBUG, 
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info("Starting MGCalibrator run")
    logging.info(f"Using up to {args.threads or 'all available'} threads.")

    # Get list of BAM files
    bam_files = list_bam_files(
        folder=args.bam_folder,
        extensions=args.extensions,
        suffixes=args.suffixes,
    )
    if not bam_files:
        raise ValueError(f"No BAM files found in {args.bam_folder} with specified extensions.")
    
    logging.info(f"Found {len(bam_files)} BAM files.")

    # Get variables from arguments
    dna_mass_df = pd.read_csv(args.dna_mass)
    initial_dna_mass = dict(zip(dna_mass_df.sample_id, dna_mass_df.DNA_mass))
    reference_clusters_csv = args.reference_clusters
    reference_bins_csv = args.reference_bins
    min_read_perc_identity = args.perc_ident
    scaling_factors_dict = args.scaling_factors
    output_dir = os.path.dirname(args.output)
    
    if not args.skip_filtering:
        # Filter BAM files
        logging.info(f"Filtering BAM files with {min_read_perc_identity}% identity...")
        
        bam_files_filtered = run_coverm_filter(bam_files, min_read_perc_identity, output_dir)
        
        logging.info(f"Filtering completed.")

    else:
        bam_files_filtered = bam_files
        # Index BAM files that are not indexed
        for bam_file in bam_files_filtered:
            if not os.path.exists(f"{bam_file}.bai"):
                try:
                    # Index the filtered BAM file using samtools
                    index_cmd = ["samtools", "index", bam_file]
                    subprocess.run(index_cmd, check=True)
                except subprocess.CalledProcessError as e:
                    # Handle errors if the samtools command fails
                    logging.info(f"Error occurred while processing {bam_file}: {e}")

    # Compute raw depths with error
    logging.info(f"Start computing raw depths with errors...")

    depths_with_errors = compute_raw_depths_with_error(bam_files_filtered, 
                                                        reference_clusters_csv=reference_clusters_csv,
                                                        reference_bins_csv=reference_bins_csv, 
                                                        n_simulations=args.mc_samples, 
                                                        pseudocount=args.pseudocount,
                                                        batch_size=args.batch_size,
                                                        n_jobs=args.threads)

    # Get sample names
    samples = set(depths_with_errors.Sample)

    # Load or calculate scaling factors    
    scaling_factors = calculate_scaling_factors(samples, bam_files, initial_dna_mass, scaling_factors_dict, n_workers=args.threads)

    logging.info("Calibrating depths...")

    # Maak een Series van de scaling_factors met dezelfde index als depths_with_errors
    depths_with_errors["scaling_factor"] = depths_with_errors["Sample"].map(scaling_factors)

    # Calibreer direct, zonder loops:
    depths_with_errors["calibrated_depth"] = depths_with_errors["IQM_mean"] * depths_with_errors["scaling_factor"]
    depths_with_errors["calibrated_lower_ci"] = depths_with_errors["IQM_lower_ci"] * depths_with_errors["scaling_factor"]
    depths_with_errors["calibrated_upper_ci"] = depths_with_errors["IQM_upper_ci"] * depths_with_errors["scaling_factor"]

    # Save final output to CSV
    depths_with_errors.to_csv(args.output, index=False)

    logging.info(f"Results saved to '{args.output}'")
    
    # Final message
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    RESET = "\033[0m"
    logging.info(GREEN + "="*45 + RESET)
    logging.info(CYAN + "🌟 MGCalibrator has finished successfully! 🌟" + RESET)
    logging.info(GREEN + "="*45 + RESET)
    logging.info(CYAN + "    Thank you for choosing MGCalibrator!\n" + RESET)

#     # Generate plots if requested
#     if args.plot:
#         plot_output_base = args.output.replace('.csv', '')
#         plot_absolute_abundances_with_ci(
#             raw_depths, absolute, lower_ci, upper_ci,
#             plot_output_base, args.top_features, args.plot_format, args.figsize
#         )
#         logging.info(f"  - Plots saved: {plot_output_base}_*.{args.plot_format}")

# def create_consolidated_output(raw_depths_df, absolute_df, lower_ci_df, upper_ci_df, zero_replaced_df, scaling_factors):
#     """
#     Create a consolidated DataFrame with original raw_depths, absolute abundances, confidence intervals,
#     zero-replaced raw_depths, and scaling factors for full transparency.
    
#     Output format:
#     - Rows: features
#     - Columns: sample_name_raw_depths, sample_name_zero_replaced, sample_name_absolute, 
#                sample_name_lower_95ci, sample_name_upper_95ci, sample_name_scaling_factor
#     """
#     consolidated_columns = []
#     consolidated_data = []
    
#     for sample in raw_depths_df.columns:
#         # Add columns for this sample
#         consolidated_columns.extend([
#             f"{sample}_raw_depths",
#             f"{sample}_zero_replaced",
#             f"{sample}_absolute", 
#             f"{sample}_lower_95ci",
#             f"{sample}_upper_95ci",
#             f"{sample}_scaling_factor"
#         ])
        
#         # Create scaling factor column (same value for all features in this sample)
#         scaling_factor_column = np.full(len(raw_depths_df), scaling_factors[sample])
        
#         # Add data for this sample
#         if len(consolidated_data) == 0:
#             # Initialize with first sample
#             consolidated_data = [
#                 raw_depths_df[sample].values,
#                 zero_replaced_df[sample].values,
#                 absolute_df[sample].values,
#                 lower_ci_df[sample].values, 
#                 upper_ci_df[sample].values,
#                 scaling_factor_column
#             ]
#         else:
#             # Append additional samples
#             consolidated_data.extend([
#                 raw_depths_df[sample].values,
#                 zero_replaced_df[sample].values,
#                 absolute_df[sample].values,
#                 lower_ci_df[sample].values,
#                 upper_ci_df[sample].values,
#                 scaling_factor_column
#             ])
    
#     # Transpose to get correct shape
#     consolidated_array = np.array(consolidated_data).T
    
#     return pd.DataFrame(
#         consolidated_array,
#         index=raw_depths_df.index,
#         columns=consolidated_columns
#     )

# def plot_absolute_abundances_with_ci(raw_depths_df, absolute_df, lower_ci_df, upper_ci_df, 
#                                    output_base, top_features=20, plot_format="png", figsize=[12, 8]):
#     """
#     Create publication-quality plots showing absolute abundances with confidence intervals.
    
#     Parameters:
#     -----------
#     raw_depths_df : pd.DataFrame
#         Original raw_depths
#     absolute_df : pd.DataFrame  
#         Absolute abundances (point estimates)
#     lower_ci_df : pd.DataFrame
#         Lower 95% confidence intervals
#     upper_ci_df : pd.DataFrame
#         Upper 95% confidence intervals
#     output_base : str
#         Base filename for output plots
#     top_features : int
#         Number of top abundant features to display
#     plot_format : str
#         Output format (png, pdf, svg)
#     figsize : list
#         Figure size [width, height]
#     """
#     logging.info(f"\n=== GENERATING PUBLICATION-QUALITY PLOTS ===")
    
#     # Set publication-style plot parameters
#     plt.style.use('default')
#     sns.set_palette("husl")
#     plt.rcParams.update({
#         'font.size': 12,
#         'axes.titlesize': 14,
#         'axes.labelsize': 12,
#         'xtick.labelsize': 10,
#         'ytick.labelsize': 10,
#         'legend.fontsize': 10,
#         'figure.titlesize': 16
#     })
    
#     # Get top features by mean absolute abundance
#     mean_abundance = absolute_df.mean(axis=1).sort_values(ascending=False)
#     top_feature_names = mean_abundance.head(top_features).index
    
#     logging.info(f"Plotting top {len(top_feature_names)} features by mean absolute abundance")
    
#     # Plot 1: Bar plot with error bars for each sample
#     n_samples = len(absolute_df.columns)
#     fig, axes = plt.subplots(1, n_samples, figsize=(figsize[0] * n_samples/2, figsize[1]), 
#                             squeeze=False)
#     axes = axes.flatten()
    
#     for i, sample in enumerate(absolute_df.columns):
#         ax = axes[i]
        
#         # Get data for this sample and top features
#         sample_absolute = absolute_df.loc[top_feature_names, sample]
#         sample_lower = lower_ci_df.loc[top_feature_names, sample]
#         sample_upper = upper_ci_df.loc[top_feature_names, sample]
        
#         # Calculate error bars
#         lower_err = sample_absolute - sample_lower
#         upper_err = sample_upper - sample_absolute
        
#         # Create bar plot
#         x_pos = np.arange(len(top_feature_names))
#         bars = ax.bar(x_pos, sample_absolute, yerr=[lower_err, upper_err], 
#                      capsize=5, alpha=0.8, color=sns.color_palette("husl", len(top_feature_names)))
        
#         ax.set_title(f'Sample: {sample}', fontweight='bold')
#         ax.set_xlabel('Features')
#         ax.set_ylabel('Absolute Abundance')
#         ax.set_xticks(x_pos)
#         ax.set_xticklabels([name[:15] + '...' if len(name) > 15 else name 
#                            for name in top_feature_names], rotation=45, ha='right')
#         ax.grid(True, alpha=0.3)
        
#         # Add value labels on bars
#         for j, (bar, val) in enumerate(zip(bars, sample_absolute)):
#             if val > 0:
#                 ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + upper_err.iloc[j]*0.1,
#                        f'{val:.1e}', ha='center', va='bottom', fontsize=8)
    
#     plt.tight_layout()
#     plt.savefig(f"{output_base}_barplot_by_sample.{plot_format}", dpi=300, bbox_inches='tight')
#     plt.close()
    
#     # Plot 2: Heatmap of absolute abundances
#     fig, ax = plt.subplots(figsize=figsize)
    
#     heatmap_data = absolute_df.loc[top_feature_names].T  # Samples as rows, features as columns
    
#     # Use log scale for better visualization if values span many orders of magnitude
#     log_data = np.log10(heatmap_data + 1e-10)  # Add small constant to handle zeros
    
#     sns.heatmap(log_data, annot=False, cmap='viridis', 
#                 xticklabels=[name[:20] + '...' if len(name) > 20 else name for name in top_feature_names],
#                 yticklabels=heatmap_data.index, ax=ax, cbar_kws={'label': 'log10(Absolute Abundance)'})
    
#     ax.set_title('Heatmap of Absolute Abundances (log10 scale)', fontweight='bold', pad=20)
#     ax.set_xlabel('Features')
#     ax.set_ylabel('Samples')
#     plt.xticks(rotation=45, ha='right')
#     plt.yticks(rotation=0)
    
#     plt.tight_layout()
#     plt.savefig(f"{output_base}_heatmap.{plot_format}", dpi=300, bbox_inches='tight')
#     plt.close()
    
#     # Plot 3: Confidence interval comparison across samples for top 5 features
#     top_5_features = top_feature_names[:5]
#     fig, ax = plt.subplots(figsize=figsize)
    
#     x_offset = np.linspace(-0.3, 0.3, len(absolute_df.columns))
#     colors = sns.color_palette("husl", len(absolute_df.columns))
    
#     for i, sample in enumerate(absolute_df.columns):
#         sample_pos = np.arange(len(top_5_features)) + x_offset[i]
#         sample_absolute = absolute_df.loc[top_5_features, sample]
#         sample_lower = lower_ci_df.loc[top_5_features, sample]
#         sample_upper = upper_ci_df.loc[top_5_features, sample]
        
#         lower_err = sample_absolute - sample_lower
#         upper_err = sample_upper - sample_absolute
        
#         ax.errorbar(sample_pos, sample_absolute, yerr=[lower_err, upper_err],
#                    fmt='o', capsize=5, label=sample, color=colors[i], markersize=8)
    
#     ax.set_xlabel('Features')
#     ax.set_ylabel('Absolute Abundance')
#     ax.set_title('Top 5 Features: Absolute Abundances with 95% Confidence Intervals', fontweight='bold')
#     ax.set_xticks(range(len(top_5_features)))
#     ax.set_xticklabels([name[:25] + '...' if len(name) > 25 else name for name in top_5_features],
#                       rotation=45, ha='right')
#     ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
#     ax.grid(True, alpha=0.3)
#     ax.set_yscale('log')
    
#     plt.tight_layout()
#     plt.savefig(f"{output_base}_top5_comparison.{plot_format}", dpi=300, bbox_inches='tight')
#     plt.close()
    
#     logging.info(f"✓ Generated 3 plots:")
#     logging.info(f"  - Bar plots by sample: {output_base}_barplot_by_sample.{plot_format}")
#     logging.info(f"  - Heatmap: {output_base}_heatmap.{plot_format}")
#     logging.info(f"  - Top 5 comparison: {output_base}_top5_comparison.{plot_format}")
