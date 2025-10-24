import pandas as pd
import numpy as np
import logging
import os
import subprocess
import io
import tempfile
import pysam
import concurrent.futures
from functools import reduce
from .parser import calculate_total_base_pairs

def _get_depth_IQM(depth_list):
    depths = np.sort(np.array(depth_list))
    n = len(depths)
    q1_idx = int(0.25 * n)
    q3_idx = int(0.75 * n)
    iqr_depths = depths[q1_idx:q3_idx+1]
    iqm = np.mean(iqr_depths)
    return iqm

def run_coverm_filter(bam_files, min_read_perc_identity, output_dir):
    # List to store the output BAM filenames
    output_bam_files = []

    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each BAM file
    for bam_file in bam_files:
        # Determine the output BAM filename by appending '_filtered' to the original filename
        filename = os.path.basename(bam_file)  # Get the filename from the full path
        output_bam_name = os.path.join(output_dir, filename.replace(".sorted.bam", f"_filtered{min_read_perc_identity}.bam"))

        try:
            # Run CoverM filter with the new output BAM file
            coverm_cmd = [
                "coverm", "filter", 
                "-b", bam_file,       # Input BAM file
                "-o", output_bam_name,  # Output BAM file
                "--min-read-percent-identity", str(min_read_perc_identity)  # Minimum read identity threshold
            ]
            # Execute the CoverM command and check for errors
            subprocess.run(coverm_cmd, check=True)

            # Index the filtered BAM file using samtools
            index_cmd = ["samtools", "index", output_bam_name]
            subprocess.run(index_cmd, check=True)

            # Add the output BAM filename to the list of results
            output_bam_files.append(output_bam_name)

        except subprocess.CalledProcessError as e:
            # Handle errors if the CoverM or samtools command fails
            print(f"Error occurred while processing {bam_file}: {e}")
            # You can choose to either re-raise the error or handle it accordingly

    # Return the list of new BAM files
    return output_bam_files

def get_depth_dict_with_samtools(bam_path):
    # Roep samtools depth aan
    cmd = ["samtools", "depth", "-a", bam_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # Lees de uitvoer direct in als dataframe
    df = pd.read_csv(proc.stdout, sep='\t', names=['ref', 'pos', 'depth'])
    
    # Zet om naar dict van numpy arrays per referentie
    depth_dict = {}
    for ref, group in df.groupby('ref'):
        arr = np.zeros(group['pos'].max(), dtype=np.int16)
        arr[group['pos'].values - 1] = group['depth'].values
        depth_dict[ref] = arr
    return depth_dict
    
def get_reads_dict_from_bam(bam_path):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    reads_dict = {ref: [] for ref in bamfile.references}
    for read in bamfile.fetch():
        ref = bamfile.get_reference_name(read.reference_id)
        start = read.reference_start
        length = read.query_length
        reads_dict[ref].append((start, length))
    # Zet naar numpy arrays
    for ref in reads_dict:
        if reads_dict[ref]:
            reads_dict[ref] = np.array(reads_dict[ref], dtype=np.int32)
        else:
            reads_dict[ref] = np.empty((0, 2), dtype=np.int32)
    bamfile.close()
    return reads_dict

def apply_clustering_to_dicts(reads_dict, depth_dict, clusters_csv_path):
    """
    Voegt referenties samen tot clusters volgens een csv-bestand.
    """
    ref_clusters = pd.read_csv(clusters_csv_path)
    ref_clusters.columns = ["sequence", "cluster"]

    # Loop over bins
    for cluster_name in ref_clusters["cluster"].unique():
        # Initiate lists
        list_of_depth_lists = []
        list_of_reads_arrays = []
        # Loop over sequences in deze cluster
        for i, sequence in enumerate(ref_clusters.loc[ref_clusters["cluster"] == cluster_name, "sequence"]):
            # Add depth_list to list_of_depth_lists
            list_of_depth_lists.append(depth_dict[sequence])
            list_of_reads_arrays.append(reads_dict[sequence])
            # Verwijder individuele sequence uit dicts
            del reads_dict[sequence]
            del depth_dict[sequence]
        # Sum depths of all depth lists and add to depth dictionary
        summed_depths = np.zeros(max(len(depth_list) for depth_list in list_of_depth_lists))
        for depth_list in list_of_depth_lists:
            summed_depths[:len(depth_list)] += depth_list
        depth_dict[cluster_name] = summed_depths
        # Stack reads arrays and add to dictionary
        reads_dict[cluster_name] = np.vstack(list_of_reads_arrays)
    return reads_dict, depth_dict

def apply_binning_to_dicts(reads_dict, depth_dict, bam_path, bins_csv_path):
    """
    Voegt referenties samen tot bins volgens een csv-bestand.
    """
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    ref_bins = pd.read_csv(bins_csv_path)
    ref_bins.columns = ["sequence", "bin"]
    ref_bins["ref_length"] = ""

    # Voeg lengte van referenties toe aan ref_bins
    for ref in bamfile.references:
        if ref in list(ref_bins["sequence"]):
            ref_length = bamfile.get_reference_length(ref)
            ref_bins.loc[ref_bins["sequence"] == ref, "ref_length"] = ref_length

    # Loop over bins
    for bin_name in ref_bins["bin"].unique():
        # Loop over sequences in deze bin
        for i, sequence in enumerate(ref_bins.loc[ref_bins["bin"] == bin_name, "sequence"]):
            if i == 0:
                # READS_DICT: eerste sequence wordt basis van bin
                reads_dict[bin_name] = reads_dict[sequence]
                pos_displacement = int(ref_bins.loc[ref_bins["sequence"] == sequence, "ref_length"].values[0])
                # DEPTH_DICT
                depth_dict[bin_name] = depth_dict[sequence]
            else:
                # READS_DICT: voeg displacement toe aan starts
                reads_dict[sequence] = [(int(read[0] + pos_displacement), read[1]) for read in reads_dict[sequence]]
                reads_array = np.array(reads_dict[sequence], dtype=np.int32)
                if reads_array.ndim == 1:
                    reads_array = reads_array.reshape((0, 2))
                reads_dict[bin_name] = np.vstack([reads_dict[bin_name], reads_array])
                pos_displacement += int(ref_bins.loc[ref_bins["sequence"] == sequence, "ref_length"].values[0])
                # DEPTH_DICT: concateneer depths
                depth_dict[bin_name] = np.concat([depth_dict[bin_name], depth_dict[sequence]])
            # Verwijder individuele sequence uit dicts
            del reads_dict[sequence]
            del depth_dict[sequence]
    bamfile.close()
    return reads_dict, depth_dict

def MC_simulation_for_IQM_depth(depths, reads, n_simulations=1000, pseudocount=1, batch_size=500):
    """
    depths: np.array, shape = (genoomlengte,)
    reads: np.array van shape (n_reads, 2), met columns [start, length]
    """

    n_mapped_reads = reads.shape[0]
    sim_total_reads = np.random.poisson(n_mapped_reads, n_simulations)
    reads_per_simulation = sim_total_reads - n_mapped_reads

    IQM_list = []
    n_batches = int(np.ceil(n_simulations / batch_size))

    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_simulations)
        batch_reads_per_sim = reads_per_simulation[start_idx:end_idx]
        batch_size_actual = end_idx - start_idx

        # Maak in één keer een matrix met kopieën van depths (int16 om geheugen te besparen)
        batch_depths = np.tile(depths, (batch_size_actual, 1)).astype(np.int16)
        # reads is meestal klein genoeg om gewoon per simulatie te kopiëren
        for i, n in enumerate(batch_reads_per_sim):
            sim_reads = reads.copy()
            sim_depths = batch_depths[i]

            n_reads = abs(n)

            if n_reads == 0:
                IQM_list.append(_get_depth_IQM(sim_depths))
                continue

            if n < 0:
                # Randomly select n_reads reads to remove
                if n_reads > sim_reads.shape[0]:
                    n_reads = sim_reads.shape[0]
                indices_to_remove = np.random.choice(sim_reads.shape[0], size=n_reads, replace=False)
                reads_to_remove = sim_reads[indices_to_remove]
                
                # Update depths: subtract for each removed read
                for start, length in reads_to_remove:
                    sim_depths[start:start+length] -= 1

                iqm = _get_depth_IQM(sim_depths)
                IQM_list.append(iqm)

            if n > 0:
                if len(depths) > 151:
                    valid_indices = np.arange(0, len(depths) - 151)
                    probs = np.full(len(depths), pseudocount, dtype=float)

                    starts = sim_reads[:, 0].astype(np.intp)
                    np.add.at(probs, starts, 1)
                    probs = probs[:-151]
                    probs = probs / probs.sum()

                    positions_to_add = np.random.choice(valid_indices, size=n_reads, p=probs, replace=True)
                
                else:
                    # Only use first position to add reads
                    positions_to_add = np.full(n_reads, 0, dtype=int)
                    
                lengths_list = sim_reads[:, 1]
                if n_reads > len(lengths_list):
                    lengths_to_add = np.random.choice(lengths_list, size=n_reads, replace=True)
                else:
                    lengths_to_add = np.random.choice(lengths_list, size=n_reads, replace=False)

                for pos_to_add, length_to_add in zip(positions_to_add, lengths_to_add):
                    region = slice(pos_to_add, pos_to_add + length_to_add)
                    sim_depths[region] += 1

                IQM_list.append(_get_depth_IQM(sim_depths))

    iqm_mean = np.mean(IQM_list)
    iqm_lower_ci = np.percentile(IQM_list, 2.5)
    iqm_upper_ci = np.percentile(IQM_list, 97.5)

    return iqm_mean, iqm_lower_ci, iqm_upper_ci

def process_bam_file(bam_file, reference_clusters_csv, reference_bins_csv, n_simulations, pseudocount, batch_size):
    filename = os.path.basename(bam_file)
    name_without_ext = filename.split('.')[0]
    sample_name = "_".join(name_without_ext.split('_')[0:2])
    
    depth_dict = get_depth_dict_with_samtools(bam_file)
    reads_dict = get_reads_dict_from_bam(bam_file)

    logging.info(f"Reads_dict and depth_dict created for {sample_name}")

    if reference_clusters_csv:
        num_refs_before = len(reads_dict.keys())
        reads_dict, depth_dict = apply_clustering_to_dicts(reads_dict, depth_dict, reference_clusters_csv)
        num_refs_after = len(reads_dict.keys())

        logging.info(f"Clustering applied. Number of references: {num_refs_before} --> {num_refs_after}")

    if reference_bins_csv:
        num_refs_before = len(reads_dict.keys())
        reads_dict, depth_dict = apply_binning_to_dicts(reads_dict, depth_dict, bam_file, reference_bins_csv)
        num_refs_after = len(reads_dict.keys())

        logging.info(f"Binning applied. Number of references: {num_refs_before} --> {num_refs_after}")

    logging.info(f"Starting {n_simulations} MC simulations in batches of {batch_size} for {sample_name}...")

    result_rows = []
    for ref in reads_dict.keys():
        if np.shape(reads_dict[ref])[0] > 0:
            IQM_mean, IQM_lower_ci, IQM_upper_ci = MC_simulation_for_IQM_depth(
                depth_dict[ref], reads_dict[ref], 
                n_simulations=n_simulations, 
                pseudocount=pseudocount, 
                batch_size=batch_size)
            result_rows.append({
                "Sample": sample_name,
                "Reference": ref,
                "IQM_mean": IQM_mean,
                "IQM_lower_ci": IQM_lower_ci,
                "IQM_upper_ci": IQM_upper_ci
            })

    logging.info(f"Monte Carlo simulation for {sample_name} done.")

    return result_rows

def compute_raw_depths_with_error(
    bam_files_filtered, reference_clusters_csv=None, reference_bins_csv=None, n_simulations=100, pseudocount=1, batch_size=500, n_jobs=4):

    result_rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = []
        for bam_file in bam_files_filtered:
            futures.append(executor.submit(
                process_bam_file, bam_file, reference_clusters_csv, reference_bins_csv, n_simulations, pseudocount, batch_size))
        for future in concurrent.futures.as_completed(futures):
            result_rows.extend(future.result())
    raw_depths_df = pd.DataFrame(result_rows)

    return raw_depths_df

def calculate_scaling_factors(samples, bam_files, initial_dna_mass, scaling_factors_dict, n_workers=None):
    scaling_factors = {}
    # Check if scaling_factors_dict exists, if True: load into dictionary
    if os.path.exists(scaling_factors_dict):
        with open(scaling_factors_dict, 'r') as f:
            for line in f:
                line = line.strip()
                # Split at the first colon
                sample, scaling_factor = line.split(':', 1)
                sample = sample.strip()
                scaling_factor = float(scaling_factor.strip())
                scaling_factors[sample] = scaling_factor

    # If not all samples are represented in the scaling_factors dictionary: calculate new scaling_factors and save to file
    if not samples <= set(scaling_factors.keys()):

        logging.info(f"Calculating scaling factors...")

        # Calculate final DNA mass from BAM files
        sample_base_pairs = calculate_total_base_pairs(bam_files=bam_files, n_workers=n_workers)
        final_dna_mass = {sample: sample_base_pairs.get(sample, 0) * 1.079e-12  # Convert to ng
                        for sample in samples}
            
        logging.debug(f"sample_base_pairs: {sample_base_pairs}")
        logging.debug(f"samples: {samples}")
        logging.debug(f"final_dna_mass: {final_dna_mass}")

        # Check for samples with no BAM data
        for sample in samples:
            if final_dna_mass[sample] == 0:
                raise ValueError(f"No base pairs found for {sample} in BAM files.")

        # Calculate scaling factors
        scaling_factors = {sample: initial_dna_mass[sample] / final_dna_mass[sample] 
                        for sample in samples}

        logging.info(f"Done calculating scaling factors.")

        # Save it to a file in the expected format:
        with open(scaling_factors_dict, "w") as file:
            for key, value in scaling_factors.items():
                file.write(f"{key}: {value}\n")

    else:
        logging.info(f"All necessary scaling factors loaded from file.")

    return scaling_factors
