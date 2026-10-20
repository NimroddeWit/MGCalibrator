import os
import gzip
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


# ----------------------------------------------------------------------
# FASTQ base pair counting
# ----------------------------------------------------------------------
def _open_fastq(fastq_path):
    """Open a (optionally gzip-compressed) FASTQ file in binary mode."""
    if fastq_path.endswith((".gz", ".gzip")):
        return gzip.open(fastq_path, "rb")
    return open(fastq_path, "rb")


def count_bases_in_fastq(fastq_path):
    """
    Count the total number of bases in a FASTQ file.

    Assumes the standard 4-line FASTQ format (header, sequence, '+', quality),
    where the sequence is on a single line. Gzip-compressed files are supported.

    Parameters
    ----------
    fastq_path : str
        Path to the FASTQ file (plain or gzip-compressed).

    Returns
    -------
    int
        Total number of bases (sum of all read lengths).
    """
    total_bases = 0
    n_lines = 0
    with _open_fastq(fastq_path) as handle:
        for line in handle:
            if n_lines == 0 and not line.startswith(b"@"):
                raise ValueError(
                    f"'{fastq_path}' does not look like a FASTQ file (first line does not start with '@')."
                )
            if n_lines % 4 == 1:
                total_bases += len(line.rstrip(b"\r\n"))
            n_lines += 1

    if n_lines % 4 != 0:
        logging.warning(
            f"'{fastq_path}' contains {n_lines} lines, which is not a multiple of 4. "
            "The file may be truncated or use multi-line sequences."
        )
    return total_bases


def _process_fastq_file(sample_name, fastq_path):
    """
    Helper function to count the bases of one FASTQ file.
    Used for parallel execution.

    Returns
    -------
    tuple
        (sample_name, fastq_path, total_bases)
    """
    logging.info(f"Counting bases in FASTQ file: {os.path.basename(fastq_path)}")
    total_bases = count_bases_in_fastq(fastq_path)
    logging.debug(f"Sample '{sample_name}': {total_bases} bases counted in '{fastq_path}'")
    return sample_name, fastq_path, total_bases


def calculate_total_base_pairs_fastq(fastq_files_by_sample, n_workers=None):
    """
    Calculate the total number of bases per sample from FASTQ files in parallel.
    All FASTQ files belonging to a sample (e.g. R1 + R2 for paired-end Illumina
    reads, or a single file for Nanopore reads) are summed.

    Parameters
    ----------
    fastq_files_by_sample : dict
        Mapping of sample name -> list of FASTQ file paths.
    n_workers : int, optional
        Number of worker processes (default: os.cpu_count()).

    Returns
    -------
    dict
        Dictionary mapping sample names to total base counts.
    """
    if n_workers is None:
        n_workers = os.cpu_count()

    n_files = sum(len(paths) for paths in fastq_files_by_sample.values())
    logging.info(
        f"Starting base counting for {n_files} FASTQ files "
        f"({len(fastq_files_by_sample)} samples) using {n_workers} workers."
    )

    sample_base_pairs = {sample: 0 for sample in fastq_files_by_sample}

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker_logging
    ) as executor:
        future_to_job = {
            executor.submit(_process_fastq_file, sample, path): (sample, path)
            for sample, paths in fastq_files_by_sample.items()
            for path in paths
        }

        for future in concurrent.futures.as_completed(future_to_job):
            sample, path = future_to_job[future]
            try:
                _, _, total_bases = future.result()
            except Exception as exc:
                # Without all files the scaling factor would be wrong, so fail loudly.
                raise RuntimeError(f"Error counting bases in '{path}': {exc}") from exc
            sample_base_pairs[sample] += total_bases

    logging.info(f"Completed base counting for {len(sample_base_pairs)} samples.")
    return sample_base_pairs