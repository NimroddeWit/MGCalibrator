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
    cmd = ["samtools", "depth", "-aa", bam_path]
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

def apply_binning_to_dicts(reads_dict, depth_dict, bam_path, bin_csv_path):
    """
    Voegt referenties samen tot bins volgens een csv-bestand.
    """
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    ref_bins = pd.read_csv(bin_csv_path)
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
                valid_indices = np.arange(0, len(depths) - 151)
                probs = np.full(len(depths), pseudocount, dtype=float)

                starts = sim_reads[:, 0].astype(np.intp)
                np.add.at(probs, starts, 1)
                probs = probs[:-151]
                probs = probs / probs.sum()

                positions_to_add = np.random.choice(valid_indices, size=n_reads, p=probs, replace=True)

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

def process_bam_file(bam_file, reference_bins_csv, n_simulations, pseudocount, batch_size):
    filename = os.path.basename(bam_file)
    name_without_ext = filename.split('.')[0]
    sample_name = "_".join(name_without_ext.split('_')[0:2])
    
    depth_dict = get_depth_dict_with_samtools(bam_file)
    reads_dict = get_reads_dict_from_bam(bam_file)

    logging.info(f"Reads_dict and depth_dict created for {sample_name}")

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
        else:
            result_rows.append({
                "Sample": sample_name,
                "Reference": ref,
                "IQM_mean": 0,
                "IQM_lower_ci": 0,
                "IQM_upper_ci": 0,
            })
    
    logging.info(f"Monte Carlo simulation for {sample_name} done.")

    return result_rows

def compute_raw_depths_with_error(
    bam_files_filtered, reference_bins_csv=None, n_simulations=100, pseudocount=1, batch_size=500, n_jobs=4):

    result_rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = []
        for bam_file in bam_files_filtered:
            futures.append(executor.submit(
                process_bam_file, bam_file, reference_bins_csv, n_simulations, pseudocount, batch_size))
        for future in concurrent.futures.as_completed(futures):
            result_rows.extend(future.result())
    raw_depths_df = pd.DataFrame(result_rows)

    return raw_depths_df

# def compute_raw_depths_with_error(
#         bam_files_filtered, reference_bins_csv=None, n_simulations=100, pseudocount=1, batch_size=500):

#     result_rows = []
    
#     for bam_file in bam_files_filtered:

#         filename = os.path.basename(bam_file)
#         sample_name = "_".join(filename.split('_')[0:2])

#         depth_dict = get_depth_dict_with_samtools(bam_file)
#         reads_dict = get_reads_dict_from_bam(bam_file)

#         logging.info(f"Reads_dict and depth_dict created for {sample_name}")

#         if reference_bins_csv:
#             num_refs_before = len(reads_dict.keys())
            
#             reads_dict, depth_dict = apply_binning_to_dicts(reads_dict, depth_dict, bam_file, reference_bins_csv)
            
#             num_refs_after = len(reads_dict.keys())

#             logging.info(f"Binning applied. Number of references: {num_refs_before} --> {num_refs_after}")
        
#         logging.info(f"Starting {n_simulations} MC simulations in batches of {batch_size} for {sample_name}...")

#         for ref in reads_dict.keys():
            
#             if np.shape(reads_dict[ref])[0] > 0:
                
#                 IQM_mean, IQM_lower_ci, IQM_upper_ci = MC_simulation_for_IQM_depth(depth_dict[ref], reads_dict[ref], 
#                                                                                    n_simulations=n_simulations, 
#                                                                                    pseudocount=pseudocount, 
#                                                                                    batch_size=batch_size)
                
#                 result_rows.append({
#                     "Sample": sample_name,
#                     "Reference": ref,
#                     "IQM_mean": IQM_mean,
#                     "IQM_lower_ci": IQM_lower_ci,
#                     "IQM_upper_ci": IQM_upper_ci
#                 })

#             else:
                
#                 result_rows.append({
#                     "Sample": sample_name,
#                     "Reference": ref,
#                     "IQM_mean": 0,
#                     "IQM_lower_ci": 0,
#                     "IQM_upper_ci": 0,
#                 })

#         logging.info(f"Monte Carlo simulation for {sample_name} done.")

#     raw_depths_df = pd.DataFrame(result_rows)

#     return raw_depths_df

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

        logging.info(f"No scaling factors provided or not all samples represented: Calculating scaling factors...")

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

        logging.info(f"Calculated scaling factors: {scaling_factors}")

        # Save it to a file in the expected format:
        with open(scaling_factors_dict, "w") as file:
            for key, value in scaling_factors.items():
                file.write(f"{key}: {value}\n")

    else:
        logging.info(f"All necessary scaling factors loaded from file: {scaling_factors}")

    return scaling_factors

def function_to_collapse_all_comments():
    
# def compute_raw_depths(bam_files, min_read_perc_identity=97, depth_calculation_method="IQM", reference_bins_csv=None):
    # """Some description here..."""
    
    

    # raw_depth_dfs = []
    
    # for bam_file in bam_files:
    #     with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as temp_bam:
    #         temp_bam_name = temp_bam.name

    #     try: 
    #         # Run CoverM filter
    #         coverm_cmd = [
    #             "coverm", "filter", 
    #             "-b", bam_file,
    #             "-o", temp_bam_name,
    #             "--min-read-percent-identity", str(min_read_perc_identity)
    #         ]
    #         subprocess.run(coverm_cmd, check=True)

    #         # Run samtools depth
    #         samtools_cmd = [
    #             # Use argument -aa to include all references in the depth file, also references without any reads mapped to them
    #             "samtools", "depth", "-aa", temp_bam_name
    #         ]
    #         samtools_result = subprocess.run(samtools_cmd, capture_output=True, text=True, check=True)

    #         # Read output to pandas DataFrame and drop 'pos' column
    #         # samtools depth output: sequence pos depth 
    #         depth_df = pd.read_csv(io.StringIO(samtools_result.stdout), sep='\t', header=None, names=['sequence', 'pos', 'depth'])
    #         depth_df.drop(columns="pos", inplace=True)

    #         logging.debug(f"depth_df['depth'].mean(): {depth_df['depth'].mean()}")

    #         # Bin reference sequences
    #         if reference_bins_csv:
    #             reference_bins_df = pd.read_csv(reference_bins_csv)
    #             reference_bins = dict(zip(reference_bins_df.iloc[:, 0], reference_bins_df.iloc[:, 1]))

    #             try:
    #                 depth_df["sequence"] = depth_df["sequence"].apply(lambda x: reference_bins[x])
    #             except Exception as exc:
    #                 logging.error(f"Applying reference bins failed: {exc}")

    #         # Group by sequences and apply depth calculation
    #         # Interquartile mean
    #         if depth_calculation_method == "IQM":
    #             raw_depth_df = pd.DataFrame(depth_df.groupby(by="sequence").apply(_get_depth_IQM)).reset_index()
    #         # Median    
    #         elif depth_calculation_method == "M":
    #             raw_depth_df = pd.DataFrame(depth_df.groupby(by="sequence").apply(np.median)).reset_index()
    #         # Poisson-Gamma
    #         # elif depth_calculation_method == "PG":
    #         #     raw_depth_df = pd.DataFrame(depth_df.groupby(by="sequence").apply(_get_depth_PG)).reset_index()
            
    #         # Set sample name as column name
    #         bam_filename = os.path.basename(bam_file)
    #         sample_name = "_".join(bam_filename.split('_')[0:2])
    #         raw_depth_df.columns = ["sequence", sample_name]

    #         raw_depth_dfs.append(raw_depth_df)

    #     finally:
    #         # Delete temp BAM file
    #         if os.path.exists(temp_bam_name):
    #             os.remove(temp_bam_name)

    # # Combine all raw_depth dataframes into single dataframe and set sequences to index
    # raw_depth_df = reduce(lambda left, right: pd.merge(left, right, on="sequence", how="outer"), raw_depth_dfs)
    # raw_depth_df.set_index("sequence", inplace=True)
    # raw_depth_df.fillna(0, inplace=True)

    # return raw_depth_df

# def _perform_mc_simulation_for_sample(args):
#     """Helper function to run MC simulation for a single sample in a separate process."""
#     sample_idx, sample_name, sample_counts, scaling_factor, total_counts, prior, n_monte_carlo = args
    
#     logging.debug(f"Starting MC simulation for sample: {sample_name}")

#     if total_counts == 0:
#         # Return index and zero array
#         return sample_idx, np.zeros((n_monte_carlo, len(sample_counts)))

#     # Special case for a single feature to avoid degenerate Dirichlet distribution.
#     # Model the count uncertainty using a Gamma distribution, which is the conjugate
#     # prior for the rate of a Poisson distribution. This is the correct model for a single count.
#     if len(sample_counts) == 1:
#         count = sample_counts[0]
#         # The posterior for a Poisson rate, with a Gamma(alpha, beta) prior, is
#         # Gamma(alpha + count, beta + 1). We use a non-informative prior where beta -> 0,
#         # giving a posterior of Gamma(count + prior, 1).
#         sampled_counts = np.random.gamma(shape=count + prior, scale=1.0, size=n_monte_carlo)
        
#         # Apply scaling factor and reshape for consistency
#         sampled_absolute = (sampled_counts * scaling_factor).reshape(-1, 1)
        
#         logging.debug(f"Finished MC simulation for single-feature sample: {sample_name}")

#         return sample_idx, sampled_absolute

#     # --- Standard multi-feature case ---
    
#     # Posterior is Dirichlet(counts + prior)
#     dirichlet_params = sample_counts + prior

#     # Monte Carlo samples from the posterior distribution
#     # Sample from a Dirichlet distribution `n_monte_carlo` times
#     sampled_proportions = np.random.dirichlet(dirichlet_params, size=n_monte_carlo)

#     # For each set of proportions, sample from a multinomial to get counts
#     # This adds the sampling noise from the sequencing process
#     sampled_counts = np.array([np.random.multinomial(int(total_counts), p) for p in sampled_proportions])

#     # Apply scaling factor to get absolute abundance
#     sampled_absolute = sampled_counts * scaling_factor
    
#     logging.debug(f"Finished MC simulation for sample: {sample_name}")
    
#     return sample_idx, sampled_absolute

# def compute_absolute_abundance_with_error(counts_df, dna_mass, bam_files, scaling_factors_dict,
#                                         n_monte_carlo=1000,
#                                         alpha=0.5,
#                                         n_workers=None):
#     """
#     Compute absolute abundance with 95% confidence intervals using a Bayesian approach.

#     The method adds a Dirichlet prior to the counts and then uses Monte Carlo sampling
#     from the posterior distribution to estimate uncertainty.

#     SCALING FACTOR ERROR PROPAGATION:
#     The scaling factor directly affects the width of confidence intervals. Samples with
#     higher scaling factors (i.e., higher initial DNA concentration relative to what was
#     sequenced) will have proportionally larger confidence intervals, reflecting
#     the increased uncertainty when extrapolating from a smaller sequenced fraction.

#     Framework:
#     1. Model counts with a Multinomial distribution.
#     2. Use a Dirichlet prior on the proportions, which is equivalent to adding a pseudocount (alpha).
#     3. The posterior is a Dirichlet distribution. We draw Monte Carlo samples of compositions from it.
#     4. These are converted to counts and scaled to get absolute abundance samples.
#     5. Confidence intervals are derived from the distribution of these absolute abundance samples.

#     Parameters:
#     -----------
#     counts_df : pd.DataFrame
#         DataFrame with features (genes/taxa) as rows and samples as columns
#     dna_mass : dict
#         Dictionary mapping sample names to DNA mass (ng)
#     bam_files : list
#         List of BAM files to calculate final DNA weight (mandatory)
#     n_monte_carlo : int
#         Number of Monte Carlo samples for confidence interval estimation
#     alpha : float
#         Dirichlet prior parameter (pseudocount to add to each feature).
#     n_workers : int, optional
#         Number of worker processes to use for parallel computation.

#     Returns:
#     --------
#     tuple: (absolute_abundance, lower_ci, upper_ci, zero_replaced_counts, scaling_factors)
#         - absolute_abundance: DataFrame with point estimates (posterior mean)
#         - lower_ci: DataFrame with lower 95% confidence bounds
#         - upper_ci: DataFrame with upper 95% confidence bounds
#         - zero_replaced_counts: DataFrame with posterior mean counts for transparency
#         - scaling_factors: Dictionary of scaling factors for each sample
#     """

#     if os.path.exists(scaling_factors_dict):
#         # Load the dictionary from a text file
#         scaling_factors_loaded_dict = {}

#         with open(scaling_factors_dict, "r") as file:
#             for line in file:
#                 # Strip any extra whitespace and split the line into key and value
#                 key, value = line.strip().split(": ", 1)

#                 scaling_factors_loaded_dict[key] = float(value)
        
#         scaling_factors = scaling_factors_loaded_dict

#         logging.debug(f"Scaling factors are loaded from file.")

#     else:
#         if not bam_files:
#             raise ValueError("BAM files are required to calculate DNA weight from base pairs")

#         # Calculate initial DNA weight for each sample
#         initial_dna_weight = dna_mass

#         # Calculate final DNA weight from BAM files
#         sample_base_pairs = calculate_total_base_pairs(bam_files, n_workers=n_workers)
#         final_dna_weight = {sample: sample_base_pairs.get(sample, 0) * 1.079e-12  # Convert to ng
#                         for sample in counts_df.columns}
        
#         logging.debug(f"sample_base_pairs: {sample_base_pairs}")
#         logging.debug(f"counts_df.columns: {counts_df.columns}")
#         logging.debug(f"final_dna_weight: {final_dna_weight}")

#         # Check for samples with no BAM data
#         for sample in counts_df.columns:
#             if final_dna_weight[sample] == 0:
#                 raise ValueError(f"No base pairs found for sample {sample} in BAM files")

#         # Calculate scaling factors
#         scaling_factors = {sample: initial_dna_weight[sample] / final_dna_weight[sample]
#                         for sample in counts_df.columns}

#         # Save it to a file in the expected format:
#         with open(scaling_factors_dict, "w") as file:
#             for key, value in scaling_factors.items():
#                 file.write(f"{key}: {value}\n")

#         logging.info(f"Calculated scaling factors: {scaling_factors}")

#     logging.info(f"Note: Samples with higher scaling factors will have proportionally larger confidence intervals")
#     logging.info(f"Using Dirichlet prior for Monte Carlo sampling with alpha={alpha}")

#     if n_workers is None:
#         n_workers = os.cpu_count()

#     logging.info(f"Running Monte Carlo simulation with {n_monte_carlo} samples using up to {n_workers} workers...")
    
#     counts_matrix = counts_df.values.astype(float)
#     prior = alpha

#     # Initialize arrays to store Monte Carlo results
#     mc_results = np.zeros((n_monte_carlo, counts_df.shape[0], counts_df.shape[1]))

#     # --- Parallel Monte Carlo Simulation ---
#     with ProcessPoolExecutor(max_workers=n_workers) as executor:
#         # Create a list of arguments for each sample's simulation
#         tasks = []
#         for sample_idx, sample in enumerate(counts_df.columns):
#             sample_counts = counts_matrix[:, sample_idx]
#             task_args = (
#                 sample_idx,
#                 sample,
#                 sample_counts,
#                 scaling_factors[sample],
#                 sample_counts.sum(),
#                 prior,
#                 n_monte_carlo
#             )
#             tasks.append(task_args)

#         logging.debug(f"tasks: {tasks}")
        
#         # Submit tasks and collect results
#         futures = [executor.submit(_perform_mc_simulation_for_sample, task) for task in tasks]
        
#         for future in as_completed(futures):
#             try:
#                 sample_idx, result_matrix = future.result()
#                 mc_results[:, :, sample_idx] = result_matrix
#             except Exception as e:
#                 logging.error(f"A task in Monte Carlo simulation failed: {e}")

#     logging.info("Monte Carlo simulation finished.")
    
#     # Compute point estimate and CIs from the MC samples.
#     logging.info("Calculating point estimates and confidence intervals from simulation results...")
#     absolute_point_values = np.mean(mc_results, axis=0)
#     lower_ci_values = np.percentile(mc_results, 2.5, axis=0)
#     upper_ci_values = np.percentile(mc_results, 97.5, axis=0)

#     # Convert back to DataFrames
#     absolute_point = pd.DataFrame(absolute_point_values, index=counts_df.index, columns=counts_df.columns)
#     lower_ci = pd.DataFrame(lower_ci_values, index=counts_df.index, columns=counts_df.columns)
#     upper_ci = pd.DataFrame(upper_ci_values, index=counts_df.index, columns=counts_df.columns)

#     # provide the posterior mean of the counts
#     # This can be seen as the "zero-replaced" or "denoised" counts.
#     posterior_mean_counts_matrix = np.zeros_like(counts_matrix)
#     for sample_idx, sample in enumerate(counts_df.columns):
#         sample_counts = counts_matrix[:, sample_idx]
#         total_counts = sample_counts.sum()
#         if total_counts > 0:
#             dirichlet_params = sample_counts + prior
#             post_mean_props = dirichlet_params / dirichlet_params.sum()
#             posterior_mean_counts_matrix[:, sample_idx] = post_mean_props * total_counts

#     zero_replaced_df = pd.DataFrame(posterior_mean_counts_matrix, index=counts_df.index, columns=counts_df.columns)

#     logging.info(f"Confidence intervals computed with proper scaling factor error propagation.")

#     if len(scaling_factors) > 1:
#         max_sf = max(scaling_factors.values())
#         min_sf = min(scaling_factors.values())
#         if min_sf > 0:
#             logging.info(f"Check: Samples with scaling factors {max_sf:.2f} vs {min_sf:.2f} "
#                          f"should have ~{max_sf/min_sf:.1f}x wider confidence intervals")

#     return absolute_point, lower_ci, upper_ci, zero_replaced_df, scaling_factors
    print("Hello world!")