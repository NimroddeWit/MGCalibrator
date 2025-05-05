# Absolutifier (Python Version)

Absolutifier is a Python-based command-line tool that computes **absolute quantification** from metagenomic data.
It combines relative abundance matrices (from reads or contigs), metadata (e.g., DNA concentration), and a fixed volume to return an **absolute abundance matrix** for downstream analysis.

---

## 📦 Installation

Clone or download the repository and install it using pip:

```bash
pip install .
```

Optionally, use a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
pip install .
```

---

## ▶️ Example Usage

Run Absolutifier with:

```bash
absolutifier \
  --counts data/PE_fastq/counts.csv \
  --meta data/PE_fastq/metadata.csv \
  --output output_pe_absolute.csv \
  --volume 500 \
  --fastq_folder data/PE_fastq \
  --suffixes _R1 _R2 \
  --extension .fastq
```

### Explanation of Parameters

| **Flag**         | **Description**                                                                 |
|------------------|---------------------------------------------------------------------------------|
| `--counts`       | CSV file with **sample rows** and **taxa columns** (e.g., `sample_id,taxa_1,...`) |
| `--meta`         | CSV file with `sample_id` and `DNA_conc` columns                               |
| `--volume`       | Fixed DNA volume (µL) to apply to all samples                                  |
| `--output`       | Output CSV file with absolute abundances                                       |
| `--fastq_folder` | *(Optional)* Folder with FASTQ/FASTA files                                     |
| `--suffixes`     | *(Optional)* List of suffixes for filtering file names (e.g., `_R1 _R2`)       |
| `--extension`    | *(Optional)* File extension (default: `.fastq`)                                |
| `--singleton`    | *(Optional)* Additional single FASTQ file(s) to include                        |

---

## 🧪 Input Format Example

**counts.csv**:

```csv
sample_id,taxa_1,taxa_2,taxa_3
sample_1,0,5,2
sample_2,2,0,10
sample_3,100,4,0
```

**metadata.csv**:

```csv
sample_id,DNA_conc
sample_1,22.5
sample_2,31.1
sample_3,29.3
```

---

## 🛠 Features in Development

- Base count summaries from FASTQ files
- Sample-wise normalization reports
- Optional debug outputs

---

Feel free to open issues or contribute to development!
