import os
import gzip
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

def _init_worker_logging():
    """
    Initializer for worker processes to configure logging.
    This is necessary for logging to work in child processes on systems
    that use 'spawn' instead of 'fork' (like Windows).
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def fast_count_base_pairs(filepath):
    """
    Count base pairs in a FASTQ or FASTA file, handling compression.
    """
    total_bp = 0
    is_fastq = any(filepath.endswith(ext) for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fastq.gzip', '.fq.gzip'])
    
    open_func = gzip.open if filepath.endswith((".gz", ".gzip")) else open
    
    with open_func(filepath, "rt") as handle:
        if is_fastq:
            # For FASTQ files, the sequence is every 4th line, starting from the 2nd
            for i, line in enumerate(handle):
                if i % 4 == 1:
                    total_bp += len(line.strip())
        else: # Assumes FASTA
            # For FASTA files, count lines that don't start with '>'
            for line in handle:
                if not line.startswith('>'):
                    total_bp += len(line.strip())
                    
    return total_bp

def _process_file_for_bp_count(filepath):
    """Helper for parallel base pair counting."""
    filename = os.path.basename(filepath)
    # The user wants to see which file is being processed.
    logging.info(f"Counting base pairs in: {filename}")
    sample_name = filename.split('_')[0]
    total_bp = fast_count_base_pairs(filepath)
    logging.debug(f"Counted {total_bp} bp in {filename} for sample {sample_name}")
    return sample_name, total_bp

def calculate_total_base_pairs(sequence_files, n_workers=None):
    """
    Calculate total base pairs from sequence files for each sample in parallel.
    
    Parameters:
    -----------
    sequence_files : list
        List of sequence file paths (FASTQ or FASTA)
    n_workers : int, optional
        Number of worker processes. Defaults to os.cpu_count().
    
    Returns:
    --------
    dict: Dictionary mapping sample names to total base pairs
    """
    if n_workers is None:
        n_workers = os.cpu_count()

    sample_base_pairs = {}
    logging.info(f"Calculating total base pairs from {len(sequence_files)} files using up to {n_workers} workers...")
    
    with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker_logging) as executor:
        future_to_file = {executor.submit(_process_file_for_bp_count, fp): fp for fp in sequence_files}
        
        for future in as_completed(future_to_file):
            filepath = future_to_file[future]
            try:
                sample_name, total_bp = future.result()
                if sample_name in sample_base_pairs:
                    sample_base_pairs[sample_name] += total_bp
                else:
                    sample_base_pairs[sample_name] = total_bp
            except Exception as exc:
                logging.error(f'{filepath} generated an exception: {exc}')

    logging.info(f"Finished calculating total base pairs for {len(sample_base_pairs)} samples.")
    return sample_base_pairs
