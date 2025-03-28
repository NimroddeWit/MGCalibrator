# Absolutifier

Absolutifier is a command-line tool designed to compute absolute quantification from metagenomic sequencing data.  
It processes large-scale paired-end or single-end FASTQ files, combines them with relative abundance matrices and sample metadata (e.g., DNA concentration), and returns an **absolute abundance matrix** for downstream analyses.

---

## 📦 Installation

To install Absolutifier, first clone the repository **with submodules**, then run the installation script:

```bash
git clone --recurse-submodules https://github.com/Fuschi/Absolutifier.git
cd Absolutifier
./install
```

> 🔧 The `install` script will automatically build the project and its dependencies.  
> Make sure you have a C++20-compatible compiler and CMake installed.

---

## ▶️ Example Usage

To run Absolutifier on a dataset of paired-end metagenomes:

```bash
./absolutifier \
  --metagenomes data/PE_fastq/ \
  --suffixes _R1,_R2 \
  --extension .fastq \
  --counts data/PE_fastq/counts.csv \
  --meta data/PE_fastq/metadata.csv \
  --volume 5000000 \
  --output test \
  --debug
```

### Explanation of parameters

| **Flag**          | **Description**                                                               |
|------------------|-------------------------------------------------------------------------------|
| `--metagenomes`  | Path to the directory containing FASTQ/FASTA files                            |
| `--suffixes`     | Comma-separated suffixes for paired-end reads (e.g., `_R1,_R2`)               |
| `--extension`    | File extension used for input files (e.g., `.fastq`, `.fq`, `.fasta`)         |
| `--counts`       | CSV file containing the abundance matrix per sample              |
| `--meta`         | CSV file with metadata including sample IDs and DNA concentrations            |
| `--volume`       | Processed volume per sample (in μL)                                           |
| `--output`       | Output file path where absolute counts will be saved                          |
| `--debug`        | *(Optional)* Enables verbose logging and saves intermediate debug files       |

---

## 🐞 Debug Mode

When `--debug` is enabled, the program generates intermediate files in a `debug/` folder, including:

- `abundance_matrix_debug.csv`
- `metadata_debug.csv`
- `base_counts_debug.csv`
- `fraction_sequenced_debug.csv`

These help with transparency and validation of each step in the quantification process.

---

Feel free to open issues or contribute to development! 
