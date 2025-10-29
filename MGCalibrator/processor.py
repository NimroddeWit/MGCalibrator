"""
processor.py

Functions for processing BAM files, computing depths, clustering, binning,
and Monte Carlo simulations.
"""

# Standard library imports
import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

# Third-party imports
import numpy as np
import pandas as pd
import pysam

# Local package imports
from .parser import calculate_total_base_pairs


def _get_depth_M98(depth_list):
    depths = np.sort(np.array(depth_list))
    n = len(depths)
    q1_idx = int(0.01 * n)
    q3_idx = int(0.99 * n)
    iqr_depths = depths[q1_idx:q3_idx+1]
    iqm = np.mean(iqr_depths)
    return iqm

def run_coverm_filter(
    bam_files: List[str],
    min_read_perc_identity: float,
    output_dir: str,
    max_workers: int = 4
) -> List[str]:
    """
    Run CoverM filtering and BAM indexing in parallel for multiple BAM files.

    Each input BAM file is filtered using the given minimum read percent identity
    and saved into the specified output directory with a modified filename.

    Parameters
    ----------
    bam_files : list of str
        Paths to the input BAM files.
    min_read_perc_identity : float
        Minimum read percent identity threshold for CoverM filtering.
    output_dir : str
        Directory where the filtered BAM files will be written.
    max_workers : int, optional
        Maximum number of parallel threads to use (default: 4).

    Returns
    -------
    list of str
        List of successfully processed output BAM file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_bam_files = []

    def process_bam(bam_file: str) -> str:
        """Run CoverM and samtools index on a single BAM file."""
        filename = os.path.basename(bam_file)
        output_bam = os.path.join(
            output_dir,
            filename.replace(
                ".sorted.bam",
                f"_filtered{min_read_perc_identity}.bam"
            )
        )

        coverm_cmd = [
            "coverm", "filter",
            "-b", bam_file,
            "-o", output_bam,
            "--min-read-percent-identity", str(min_read_perc_identity)
        ]
        index_cmd = ["samtools", "index", output_bam]

        try:
            subprocess.run(coverm_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subprocess.run(index_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            logging.info(f"Successfully processed {bam_file}")
            return output_bam

        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {bam_file}: {e}")
            return None

    # Run in parallel (I/O bound external calls → threads are efficient)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_bam, bam): bam for bam in bam_files}
        for future in as_completed(futures):
            result = future.result()
            if result:
                output_bam_files.append(result)

    return output_bam_files

def get_depth_dict_with_samtools(bam_path: str) -> Dict[str, np.ndarray]:
    """
    Compute per-reference depth arrays from a BAM file using `samtools depth`.

    Parameters
    ----------
    bam_path : str
        Path to the input BAM file.

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping reference names (contigs) to NumPy arrays of read depth
        at each position (1-based positions converted to 0-based indices).

    Notes
    -----
    - This function requires `samtools` to be available in the system PATH.
    - The `-a` flag ensures all positions, including zero-coverage, are reported.
    - Uses pandas for fast parsing and vectorized assignment for speed.
    """
    cmd = ["samtools", "depth", "-a", bam_path]
    
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise RuntimeError("samtools is not installed or not in PATH.")

    # Read samtools output into a DataFrame in one shot
    df = pd.read_csv(proc.stdout, sep='\t', names=['ref', 'pos', 'depth'], dtype={'ref': str, 'pos': int, 'depth': np.int16})

    proc.stdout.close()
    proc.wait()
    if proc.returncode != 0:
        logging.error(f"samtools depth failed for {bam_path}")
        raise subprocess.CalledProcessError(proc.returncode, cmd)

    # Convert DataFrame to dictionary of NumPy arrays, one per reference
    depth_dict: Dict[str, np.ndarray] = {}
    for ref, group in df.groupby('ref'):
        arr = np.zeros(group['pos'].max(), dtype=np.int16)
        arr[group['pos'].values - 1] = group['depth'].values
        depth_dict[ref] = arr

    return depth_dict

def get_reads_dict_from_bam(bam_path: str) -> Dict[str, np.ndarray]:
    """
    Extract per-reference read start positions and lengths from a BAM file.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file to read.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary mapping each reference (contig/chromosome name) to a
        NumPy array of shape (N, 2), where each row is:
        [read_start_position, read_length].

    Notes
    -----
    - Requires `pysam` to be installed and the BAM file to be indexed.
    - Uses memory-efficient buffering suitable for large BAM files.
    """
    reads_dict: Dict[str, list] = {}

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            # Initialize lists for all references only when needed
            for read in bamfile.fetch(until_eof=True):
                ref_id = read.reference_id
                if ref_id < 0:
                    # Unmapped reads have reference_id == -1
                    continue

                ref_name = bamfile.get_reference_name(ref_id)
                # Lazily create list for each reference
                if ref_name not in reads_dict:
                    reads_dict[ref_name] = []

                reads_dict[ref_name].append(
                    (read.reference_start, read.query_length)
                )

        # Convert lists to compact NumPy arrays
        for ref, reads in reads_dict.items():
            if reads:
                reads_dict[ref] = np.array(reads, dtype=np.int32)
            else:
                reads_dict[ref] = np.empty((0, 2), dtype=np.int32)

    except (OSError, ValueError, pysam.SamtoolsError) as e:
        logging.error(f"Failed to process BAM file '{bam_path}': {e}")
        raise

    return reads_dict

def apply_clustering_to_dicts(
    reads_dict: Dict[str, np.ndarray],
    depth_dict: Dict[str, np.ndarray],
    clusters_csv_path: str
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Merge reference-level read and depth data into clusters as defined by a CSV file.

    Parameters
    ----------
    reads_dict : dict[str, np.ndarray]
        Dictionary mapping reference names to arrays of read information
        (e.g., start position and read length per read).
    depth_dict : dict[str, np.ndarray]
        Dictionary mapping reference names to depth arrays.
    clusters_csv_path : str
        Path to a CSV file defining clusters. The file must contain
        two columns: 'sequence' and 'cluster'.

    Returns
    -------
    tuple(dict[str, np.ndarray], dict[str, np.ndarray])
        Updated dictionaries with individual references replaced
        by their corresponding merged cluster entries.

    Notes
    -----
    - Each cluster groups one or more reference sequences.
    - The resulting depth array for a cluster is the sum of the depths
      of its member references (up to the maximum sequence length).
    - The read arrays are vertically stacked for all member sequences.
    """
    try:
        ref_clusters = pd.read_csv(clusters_csv_path, usecols=[0, 1])
    except Exception as e:
        logging.error(f"Failed to read cluster CSV: {clusters_csv_path}")
        raise RuntimeError(f"Error reading {clusters_csv_path}: {e}")

    # Standardize column names
    ref_clusters.columns = ["sequence", "cluster"]

    # Pre-group sequences by cluster
    cluster_map = (
        ref_clusters.groupby("cluster")["sequence"]
        .apply(list)
        .to_dict()
    )

    # Process clusters
    for cluster_name, sequences in cluster_map.items():
        # Check that all sequences exist in both dicts
        missing_reads = [s for s in sequences if s not in reads_dict]
        missing_depths = [s for s in sequences if s not in depth_dict]
        if missing_reads or missing_depths:
            logging.warning(
                f"Skipping cluster '{cluster_name}' due to missing data: "
                f"{missing_reads or missing_depths}"
            )
            continue

        # Combine depth arrays (vectorized summation)
        max_len = max(len(depth_dict[s]) for s in sequences)
        summed_depth = np.zeros(max_len, dtype=np.float32)
        for s in sequences:
            d = depth_dict[s]
            summed_depth[:len(d)] += d

        # Combine read arrays efficiently
        reads_combined = np.vstack(
            [reads_dict[s] for s in sequences if len(reads_dict[s]) > 0]
        ) if any(len(reads_dict[s]) > 0 for s in sequences) else np.empty((0, 2), dtype=np.int32)

        # Remove individual sequences to save memory
        for s in sequences:
            reads_dict.pop(s, None)
            depth_dict.pop(s, None)

        # Add cluster-level data
        depth_dict[cluster_name] = summed_depth
        reads_dict[cluster_name] = reads_combined

        logging.debug(f"Cluster '{cluster_name}' merged ({len(sequences)} refs)")

    return reads_dict, depth_dict

def apply_binning_to_dicts(
    reads_dict: Dict[str, np.ndarray],
    depth_dict: Dict[str, np.ndarray],
    bam_path: str,
    bins_csv_path: str
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Merge reference-level read and depth data into genome bins defined by a CSV file.

    Parameters
    ----------
    reads_dict : dict[str, np.ndarray]
        Dictionary mapping reference names to arrays of read information
        (columns: [start, length]).
    depth_dict : dict[str, np.ndarray]
        Dictionary mapping reference names to depth arrays.
    bam_path : str
        Path to the input BAM file.
    bins_csv_path : str
        Path to a CSV file defining bins. Must contain two columns: 'sequence' and 'bin'.

    Returns
    -------
    tuple(dict[str, np.ndarray], dict[str, np.ndarray])
        Updated dictionaries where references have been merged into their corresponding bins.
    """
    try:
        ref_bins = pd.read_csv(bins_csv_path, usecols=[0, 1])
        ref_bins.columns = ["sequence", "bin"]
        ref_bins["ref_length"] = np.nan
    except Exception as e:
        logging.error(f"Failed to read bin CSV: {bins_csv_path}")
        raise RuntimeError(f"Error reading {bins_csv_path}: {e}")

    # Add reference lengths from BAM file
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for ref in bamfile.references:
            if ref in set(ref_bins["sequence"]):
                ref_length = bamfile.get_reference_length(ref)
                ref_bins.loc[ref_bins["sequence"] == ref, "ref_length"] = ref_length

    # Precompute mapping for fast lookup
    ref_length_map = {
        row.sequence: int(row.ref_length)
        for _, row in ref_bins.dropna(subset=["ref_length"]).iterrows()
    }

    # Group sequences by bin for efficient iteration
    bins_map = (
        ref_bins.groupby("bin")["sequence"]
        .apply(list)
        .to_dict()
    )

    # Process each bin
    for bin_name, sequences in bins_map.items():
        bin_created = False
        pos_displacement = 0

        for sequence in sequences:
            if sequence not in reads_dict or sequence not in depth_dict:
                continue

            ref_length = ref_length_map.get(sequence)
            if ref_length is None:
                logging.warning(f"No reference length found for '{sequence}'; skipping.")
                continue

            if not bin_created:
                # Initialize new bin entries
                reads_dict[bin_name] = reads_dict[sequence]
                depth_dict[bin_name] = depth_dict[sequence]
                pos_displacement = ref_length
                bin_created = True
            else:
                # Shift read positions by cumulative displacement
                reads_array = reads_dict[sequence]
                if reads_array.size > 0:
                    reads_array_shifted = reads_array.copy()
                    reads_array_shifted[:, 0] += pos_displacement
                    reads_dict[bin_name] = np.vstack([reads_dict[bin_name], reads_array_shifted])

                # Concatenate depth arrays
                depth_dict[bin_name] = np.concatenate(
                    [depth_dict[bin_name], depth_dict[sequence]]
                )

                # Update position offset
                pos_displacement += ref_length

            # Remove individual reference entries to free memory
            reads_dict.pop(sequence, None)
            depth_dict.pop(sequence, None)

        if bin_created:
            logging.debug(f"Created bin '{bin_name}' with {len(sequences)} references.")

    return reads_dict, depth_dict

def MC_simulation_for_M98_depth(
    depths: np.ndarray,
    reads: np.ndarray,
    n_simulations: int = 1000,
    pseudocount: int = 1,
    batch_size: int = 500,
) -> tuple[float, float, float]:
    """
    Monte Carlo simulation to estimate M98 depth (mean of 98% of inner values) from read coverage.

    Parameters
    ----------
    depths : np.ndarray
        Array of depth per genomic position (shape: genome_length,).
    reads : np.ndarray
        Array of mapped reads with columns [start, length] (shape: n_reads, 2).
    n_simulations : int, optional
        Number of simulations to run, by default 1000.
    pseudocount : int, optional
        Small constant added for sampling probabilities, by default 1.
    batch_size : int, optional
        Number of simulations to process in a batch, by default 500.

    Returns
    -------
    tuple[float, float, float]
        (mean M98, lower 2.5% CI, upper 97.5% CI)
    """
    n_mapped_reads = reads.shape[0]
    sim_total_reads = np.random.poisson(n_mapped_reads, n_simulations)
    reads_per_sim = sim_total_reads - n_mapped_reads

    M98_list = []
    n_batches = int(np.ceil(n_simulations / batch_size))

    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_simulations)
        batch_reads_per_sim = reads_per_sim[start_idx:end_idx]
        batch_size_actual = end_idx - start_idx

        # Create a batch matrix of depths (int16 to save memory)
        batch_depths = np.tile(depths, (batch_size_actual, 1)).astype(np.int16)

        for i, n in enumerate(batch_reads_per_sim):
            sim_depths = batch_depths[i]

            n_reads = abs(n)
            if n_reads == 0:
                M98_list.append(_get_depth_M98(sim_depths))
                continue

            if n < 0:
                # Remove random reads
                n_reads = min(n_reads, reads.shape[0])
                indices_to_remove = np.random.choice(reads.shape[0], size=n_reads, replace=False)
                remove_starts = reads[indices_to_remove, 0]
                remove_lengths = reads[indices_to_remove, 1]
                for start, length in zip(remove_starts, remove_lengths):
                    sim_depths[start:start+length] -= 1
                M98_list.append(_get_depth_M98(sim_depths))

            else:
                # Add reads
                if len(depths) > 151:
                    valid_indices = np.arange(len(depths) - 151)
                    probs = np.full(len(depths), pseudocount, dtype=float)
                    np.add.at(probs, reads[:, 0], 1)
                    probs = probs[:-151]
                    probs /= probs.sum()
                    positions_to_add = np.random.choice(valid_indices, size=n_reads, p=probs, replace=True)
                else:
                    positions_to_add = np.zeros(n_reads, dtype=int)

                read_lengths = reads[:, 1]
                lengths_to_add = np.random.choice(read_lengths, size=n_reads, replace=n_reads > len(read_lengths))
                for pos, length in zip(positions_to_add, lengths_to_add):
                    sim_depths[pos:pos+length] += 1

                M98_list.append(_get_depth_M98(sim_depths))

    m98_mean = float(np.mean(M98_list))
    m98_lower_ci = float(np.percentile(M98_list, 2.5))
    m98_upper_ci = max(m98_mean, float(np.percentile(M98_list, 97.5)))

    return m98_mean, m98_lower_ci, m98_upper_ci

def process_bam_file(
    bam_file: str,
    reference_clusters_csv: str | None,
    reference_bins_csv: str | None,
    n_simulations: int,
    pseudocount: int,
    batch_size: int
) -> list[dict]:
    """
    Process a BAM file: generate reads and depth dictionaries,
    apply clustering and binning, and run Monte Carlo simulations
    to compute M98 statistics.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.
    reference_clusters_csv : str | None
        Path to CSV defining clusters. If None, clustering is skipped.
    reference_bins_csv : str | None
        Path to CSV defining bins. If None, binning is skipped.
    n_simulations : int
        Number of Monte Carlo simulations per reference.
    pseudocount : int
        Pseudocount added for probability calculations.
    batch_size : int
        Number of simulations per batch.

    Returns
    -------
    list[dict]
        List of dictionaries with M98 results per reference.
    """
    filename = os.path.basename(bam_file)
    name_without_ext = filename.split('.')[0]
    sample_name = "_".join(name_without_ext.split('_')[0:2])

    logging.info(f"Creating reads_dict and depth_dict for {sample_name}...")
    depth_dict = get_depth_dict_with_samtools(bam_file)
    reads_dict = get_reads_dict_from_bam(bam_file)
    logging.info(f"Reads_dict and depth_dict created for {sample_name}.")

    if reference_clusters_csv:
        num_refs_before = len(reads_dict)
        reads_dict, depth_dict = apply_clustering_to_dicts(
            reads_dict, depth_dict, reference_clusters_csv
        )
        num_refs_after = len(reads_dict)
        logging.info(f"Clustering applied. Number of references: {num_refs_before} --> {num_refs_after}")

    if reference_bins_csv:
        num_refs_before = len(reads_dict)
        reads_dict, depth_dict = apply_binning_to_dicts(
            reads_dict, depth_dict, bam_file, reference_bins_csv
        )
        num_refs_after = len(reads_dict)
        logging.info(f"Binning applied. Number of references: {num_refs_before} --> {num_refs_after}")

    logging.info(f"Starting {n_simulations} MC simulations in batches of {batch_size} for {sample_name}...")

    result_rows = []
    for ref, reads_array in reads_dict.items():
        if reads_array.shape[0] == 0:
            continue
        m98_mean, m98_lower_ci, m98_upper_ci = MC_simulation_for_M98_depth(
            depth_dict[ref],
            reads_array,
            n_simulations=n_simulations,
            pseudocount=pseudocount,
            batch_size=batch_size
        )
        result_rows.append({
            "Sample": sample_name,
            "Reference": ref,
            "M98_mean": m98_mean,
            "M98_lower_ci": m98_lower_ci,
            "M98_upper_ci": m98_upper_ci,
        })

    logging.info(f"Monte Carlo simulation for {sample_name} done.")
    return result_rows

def compute_raw_depths_with_error(
    bam_files_filtered: List[str],
    reference_clusters_csv: Optional[str] = None,
    reference_bins_csv: Optional[str] = None,
    n_simulations: int = 100,
    pseudocount: int = 1,
    batch_size: int = 500,
    n_jobs: int = 4
) -> pd.DataFrame:
    """
    Compute raw depths with Monte Carlo error simulations for a list of filtered BAM files.

    Parameters
    ----------
    bam_files_filtered : List[str]
        List of paths to filtered BAM files.
    reference_clusters_csv : Optional[str]
        CSV defining clusters for references. If None, clustering is skipped.
    reference_bins_csv : Optional[str]
        CSV defining bins for references. If None, binning is skipped.
    n_simulations : int
        Number of Monte Carlo simulations per reference.
    pseudocount : int
        Pseudocount for probability adjustments in MC simulations.
    batch_size : int
        Number of simulations per batch.
    n_jobs : int
        Number of parallel processes to use.

    Returns
    -------
    pd.DataFrame
        DataFrame with raw depths and Monte Carlo statistics per reference per sample.
    """
    result_rows = []

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = [
            executor.submit(
                process_bam_file,
                bam_file,
                reference_clusters_csv,
                reference_bins_csv,
                n_simulations,
                pseudocount,
                batch_size
            )
            for bam_file in bam_files_filtered
        ]

        for future in as_completed(futures):
            result_rows.extend(future.result())

    return pd.DataFrame(result_rows)

def calculate_scaling_factors(samples, bam_files, initial_dna_mass, scaling_factors_file, n_workers=None):
    """
    Calculate scaling factors for each sample based on BAM files and initial DNA mass.
    If a scaling factors file exists, it is loaded. Otherwise, scaling factors are computed
    and saved to the file.

    Args:
        samples (iterable): Sample names to calculate scaling factors for.
        bam_files (list[str]): Paths to BAM files.
        initial_dna_mass (dict): Initial DNA mass per sample.
        scaling_factors_file (str): Path to file to read/write scaling factors.
        n_workers (int, optional): Number of parallel workers for base pair calculation.

    Returns:
        dict: Scaling factors keyed by sample name.
    """
    scaling_factors = {}

    # Load existing scaling factors if file exists
    if os.path.exists(scaling_factors_file):
        with open(scaling_factors_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or ':' not in line:
                    continue
                sample, factor = line.split(':', 1)
                scaling_factors[sample.strip()] = float(factor.strip())

    # Check if all samples are covered
    missing_samples = set(samples) - scaling_factors.keys()
    if missing_samples:
        logging.info("Calculating missing scaling factors...")

        # Calculate final DNA mass from BAM files
        sample_base_pairs = calculate_total_base_pairs(bam_files=bam_files, n_workers=n_workers)
        final_dna_mass = {sample: sample_base_pairs.get(sample, 0) * 1.079e-12  # Convert to ng
                          for sample in missing_samples}

        # Check for samples with no BAM data
        for sample, mass in final_dna_mass.items():
            if mass == 0:
                raise ValueError(f"No base pairs found for sample '{sample}' in BAM files.")

        # Compute scaling factors for missing samples
        for sample in missing_samples:
            scaling_factors[sample] = initial_dna_mass[sample] / final_dna_mass[sample]

        logging.info("Scaling factors calculated. Saving to file...")

        # Save/update scaling factors file
        with open(scaling_factors_file, "w") as f:
            for sample, factor in scaling_factors.items():
                f.write(f"{sample}: {factor}\n")

    else:
        logging.info("All necessary scaling factors loaded from file.")

    return scaling_factors
