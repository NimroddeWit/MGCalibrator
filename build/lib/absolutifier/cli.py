import argparse
import pandas as pd
import numpy as np
from .parser import read_fastq
from .processor import compute_absolute_abundance, compute_absolute_abundance_with_error
from .fileutils import list_fastq_files

def main():
    parser = argparse.ArgumentParser(description="Absolutifier - convert relative abundances to absolute")
    parser.add_argument("--counts", required=True, help="CSV file with counts")
    parser.add_argument("--meta", required=True, help="CSV file with concentrations")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--volume", type=float, required=True, help="DNA volume (µl), used for all samples")
    parser.add_argument("--fastq_folder", required=True, help="Folder containing FASTQ files (mandatory)")
    parser.add_argument("--extension", default=".fastq", help="FASTQ file extension (default: .fastq)")
    parser.add_argument("--suffixes", nargs="*", help="List of suffixes to filter in files")
    parser.add_argument("--singleton", nargs="*", help="List of singleton files to include")
    
    # Error propagation options
    parser.add_argument("--error_bars", action="store_true", 
                       help="Calculate 95% confidence intervals using Dirichlet multinomial Monte Carlo")
    parser.add_argument("--mc_samples", type=int, default=1000, 
                       help="Number of Monte Carlo instances for confidence intervals (default: 1000)")
    parser.add_argument("--alpha", type=float, default=0.1, 
                       help="Prior pseudocount for Dirichlet sampling (default: 0.1 for numerical stability)")
    
    args = parser.parse_args()

    # FASTQ files are now mandatory
    fastq_files = list_fastq_files(
        folder=args.fastq_folder,
        extension=args.extension,
        suffixes=args.suffixes,
        singleton_files=args.singleton
    )
    
    if not fastq_files:
        raise ValueError(f"No FASTQ files found in {args.fastq_folder} with extension {args.extension}")
    
    print("FASTQ files found:", fastq_files)

    counts = pd.read_csv(args.counts, index_col=0).T
    meta = pd.read_csv(args.meta)
    dna_conc = dict(zip(meta.sample_id, meta.DNA_conc))
    volume = {sample: float(args.volume) for sample in meta.sample_id}

    if args.error_bars:
        print(f"Calculating 95% confidence intervals with {args.mc_samples} Monte Carlo instances...")
        print(f"Using Dirichlet multinomial sampling with alpha={args.alpha} for numerical stability...")
        print("Zero-count taxa will have confidence intervals that include zero.")
        
        absolute, lower_ci, upper_ci = compute_absolute_abundance_with_error(
            counts, dna_conc, volume, fastq_files,
            n_monte_carlo=args.mc_samples,
            alpha=args.alpha
        )
        
        # Create consolidated output with all information
        consolidated_df = create_consolidated_output(counts, absolute, lower_ci, upper_ci)
        consolidated_df.to_csv(args.output)
        
        print(f"Results saved:")
        print(f"  - Consolidated results (counts, absolute, CI): {args.output}")
        
    else:
        absolute = compute_absolute_abundance(counts, dna_conc, volume, fastq_files)
        
        # Create simple consolidated output without confidence intervals
        consolidated_df = create_simple_consolidated_output(counts, absolute)
        consolidated_df.to_csv(args.output)
        print(f"Absolute abundances saved to: {args.output}")

def create_consolidated_output(counts_df, absolute_df, lower_ci_df, upper_ci_df):
    """
    Create a consolidated DataFrame with original counts, absolute abundances, and confidence intervals.
    
    Output format:
    - Rows: taxa
    - Columns: sample_name_counts, sample_name_absolute, sample_name_lower_95ci, sample_name_upper_95ci
    """
    consolidated_columns = []
    consolidated_data = []
    
    for sample in counts_df.columns:
        # Add columns for this sample
        consolidated_columns.extend([
            f"{sample}_counts",
            f"{sample}_absolute", 
            f"{sample}_lower_95ci",
            f"{sample}_upper_95ci"
        ])
        
        # Add data for this sample
        if len(consolidated_data) == 0:
            # Initialize with first sample
            consolidated_data = [
                counts_df[sample].values,
                absolute_df[sample].values,
                lower_ci_df[sample].values, 
                upper_ci_df[sample].values
            ]
        else:
            # Append additional samples
            consolidated_data.extend([
                counts_df[sample].values,
                absolute_df[sample].values,
                lower_ci_df[sample].values,
                upper_ci_df[sample].values
            ])
    
    # Transpose to get correct shape
    consolidated_array = np.array(consolidated_data).T
    
    return pd.DataFrame(
        consolidated_array,
        index=counts_df.index,
        columns=consolidated_columns
    )

def create_simple_consolidated_output(counts_df, absolute_df):
    """
    Create a simple consolidated DataFrame with original counts and absolute abundances.
    
    Output format:
    - Rows: taxa  
    - Columns: sample_name_counts, sample_name_absolute
    """
    consolidated_columns = []
    consolidated_data = []
    
    for sample in counts_df.columns:
        # Add columns for this sample
        consolidated_columns.extend([
            f"{sample}_counts",
            f"{sample}_absolute"
        ])
        
        # Add data for this sample
        if len(consolidated_data) == 0:
            # Initialize with first sample
            consolidated_data = [
                counts_df[sample].values,
                absolute_df[sample].values
            ]
        else:
            # Append additional samples
            consolidated_data.extend([
                counts_df[sample].values,
                absolute_df[sample].values
            ])
    
    # Transpose to get correct shape
    consolidated_array = np.array(consolidated_data).T
    
    return pd.DataFrame(
        consolidated_array,
        index=counts_df.index,
        columns=consolidated_columns
    )
