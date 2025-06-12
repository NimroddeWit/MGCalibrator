from Bio import SeqIO
import os

def read_fastq(filepath):
    return SeqIO.parse(filepath, "fastq")

def calculate_total_base_pairs(fastq_files):
    """
    Calculate total base pairs from FASTQ files for each sample.
    
    Parameters:
    -----------
    fastq_files : list
        List of FASTQ file paths
    
    Returns:
    --------
    dict: Dictionary mapping sample names to total base pairs
    """
    sample_base_pairs = {}
    
    for filepath in fastq_files:
        # Extract sample name from filename (assumes format like sample_X_R1.fastq)
        filename = os.path.basename(filepath)
        sample_name = filename.split('_')[0] + '_' + filename.split('_')[1]  # e.g., sample_1
        
        # Count base pairs in this file
        total_bp = 0
        with open(filepath, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_bp += len(record.seq)
        
        # Add to sample total (handles multiple files per sample like R1, R2)
        if sample_name in sample_base_pairs:
            sample_base_pairs[sample_name] += total_bp
        else:
            sample_base_pairs[sample_name] = total_bp
    
    return sample_base_pairs
