import os
import logging
import pysam
import concurrent.futures

def _init_worker_logging():
    """
    Initialize logging configuration for worker processes.
    This is necessary to ensure logging works properly in child processes,
    especially on platforms like Windows that use 'spawn' instead of 'fork'.
    """
    logging.basicConfig(
        level=logging.INFO,  # Adjust level as needed (DEBUG/INFO)
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def _count_base_pairs_in_bam(bam_path):
    """
    Count the total number of base pairs in a BAM file, excluding
    secondary and supplementary alignments.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.

    Returns
    -------
    int
        Total number of base pairs in the BAM file.
    """
    total_base_pairs = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            if not read.is_secondary and not read.is_supplementary:
                total_base_pairs += read.query_length
    return total_base_pairs

def _process_bam_file(bam_path):
    """
    Helper function to process a BAM file and count base pairs.
    Used for parallel execution.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.

    Returns
    -------
    tuple
        (sample_name, total_base_pairs)
    """
    filename = os.path.basename(bam_path)
    logging.info(f"Processing BAM file: {filename}")

    # Extract sample name from filename (first two underscore-separated parts)
    name_without_ext = filename.split('.')[0]
    sample_name = "_".join(name_without_ext.split('_')[0:2])
    total_base_pairs = _count_base_pairs_in_bam(bam_path)

    logging.debug(f"Sample '{sample_name}': {total_base_pairs} base pairs counted in '{filename}'")
    return sample_name, total_base_pairs

def calculate_total_base_pairs(bam_files, n_workers=None):
    """
    Calculate total base pairs for each sample from a list of BAM files in parallel.

    Parameters
    ----------
    bam_files : list of str
        List of BAM file paths.
    n_workers : int, optional
        Number of worker processes (default: os.cpu_count()).

    Returns
    -------
    dict
        Dictionary mapping sample names to total base pair counts.
    """
    if n_workers is None:
        n_workers = os.cpu_count()

    logging.info(f"Starting base pair calculation for {len(bam_files)} BAM files using {n_workers} workers.")

    sample_base_pairs = {}

    # Use ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker_logging
    ) as executor:
        # Submit all BAM files for processing
        future_to_bam = {executor.submit(_process_bam_file, bam_path): bam_path for bam_path in bam_files}

        for future in concurrent.futures.as_completed(future_to_bam):
            bam_path = future_to_bam[future]
            try:
                sample_name, total_base_pairs = future.result()
                sample_base_pairs[sample_name] = total_base_pairs
            except Exception as exc:
                logging.error(f"Error processing '{bam_path}': {exc}")

    logging.info(f"Completed base pair calculation for {len(sample_base_pairs)} samples.")
    return sample_base_pairs

# import os
# import logging
# import pysam
# import concurrent.futures
# #from concurrent.futures import ProcessPoolExecutor, as_completed

# def _init_worker_logging():
#     """
#     Initializer for worker processes to configure logging.
#     This is necessary for logging to work in child processes on systems
#     that use 'spawn' instead of 'fork' (like Windows).
#     """
#     logging.basicConfig(
#         level=logging.DEBUG,
#         format='%(asctime)s - %(levelname)s - %(message)s',
#         datefmt='%Y-%m-%d %H:%M:%S'
#     )

# def _fast_count_base_pairs(filepath):
#     """
#     Count base pairs in a BAM file, handling compression.
#     """
#     bamfile = pysam.AlignmentFile(filepath, "rb")
#     total_bp = 0

#     for read in bamfile.fetch(until_eof=True):
#         # Flag-explained:
#         # not secondary (0x100) AND not supplementary (0x800)
#         if not read.is_secondary and not read.is_supplementary:
#             total_bp += read.query_length

#     bamfile.close()
                    
#     return total_bp

# def _process_file_for_bp_count(filepath):
#     """Helper for parallel base pair counting."""
#     filename = os.path.basename(filepath)
    
#     logging.info(f"Counting base pairs in: {filename}")
    
#     sample_name = "_".join(filename.split('_')[0:2])
#     total_bp = _fast_count_base_pairs(filepath)
    
#     logging.debug(f"Counted {total_bp} bp in {filename} for sample {sample_name}")
    
#     return sample_name, total_bp

# def calculate_total_base_pairs(bam_files, n_workers=None):
#     """
#     Calculate total base pairs from BAM files for each sample in parallel.
    
#     Parameters:
#     -----------
#     bam_files : list
#         List of BAM file paths 
#     n_workers : int, optional
#         Number of worker processes. Defaults to os.cpu_count().
    
#     Returns:
#     --------
#     dict: Dictionary mapping sample names to total base pairs
#     """
#     if n_workers is None:
#         n_workers = os.cpu_count()

#     sample_base_pairs = {}
#     logging.info(f"Calculating total base pairs from {len(bam_files)} samples using up to {n_workers} workers...")
    
#     with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker_logging) as executor:
#         future_to_file = {executor.submit(_process_file_for_bp_count, fp): fp for fp in bam_files}
        
#         for future in concurrent.futures.as_completed(future_to_file):
#             filepath = future_to_file[future]
#             try:
#                 sample_name, total_bp = future.result()
#                 sample_base_pairs[sample_name] = total_bp
#             except Exception as exc:
#                 logging.error(f'{filepath} generated an exception: {exc}')

#     logging.info(f"Finished calculating total base pairs for {len(sample_base_pairs)} samples.")

#     return sample_base_pairs
