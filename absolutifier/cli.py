import argparse
import pandas as pd
from .parser import read_fastq
from .processor import compute_absolute_abundance
from .fileutils import list_fastq_files

def main():
    parser = argparse.ArgumentParser(description="Absolutifier - converti abbondanze relative in assolute")
    parser.add_argument("--counts", required=True, help="File CSV con i conteggi")
    parser.add_argument("--meta", required=True, help="File CSV con le concentrazioni")
    parser.add_argument("--output", required=True, help="File CSV di output")
    parser.add_argument("--volume", type=float, required=True, help="Volume DNA (µl), usato per tutti i campioni")
    parser.add_argument("--extension", default=".fastq", help="Estensione dei file FASTQ (default: .fastq)")
    parser.add_argument("--suffixes", nargs="*", help="Lista di suffissi da filtrare nei file")
    parser.add_argument("--singleton", nargs="*", help="Lista di file singleton da includere")
    parser.add_argument("--fastq_folder", help="Cartella contenente i file FASTQ")
    args = parser.parse_args()

    if args.fastq_folder:
        fastq_files = list_fastq_files(
            folder=args.fastq_folder,
            extension=args.extension,
            suffixes=args.suffixes,
            singleton_files=args.singleton
        )
        print("File FASTQ trovati:", fastq_files)

    counts = pd.read_csv(args.counts, index_col=0).T
    meta = pd.read_csv(args.meta)
    dna_conc = dict(zip(meta.sample_id, meta.DNA_conc))
    volume = {sample: float(args.volume) for sample in meta.sample_id}

    absolute = compute_absolute_abundance(counts, dna_conc, volume)
    absolute.to_csv(args.output)
