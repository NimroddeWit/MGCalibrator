from Bio import SeqIO

def read_fastq(filepath):
    return SeqIO.parse(filepath, "fastq")
