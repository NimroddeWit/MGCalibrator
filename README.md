# Absolutifier

Raw feature counts from metagenomic sequencing are not directly comparable across samples due to variations in sequencing depth and initial DNA input. `Absolutifier` addresses this by converting raw counts into **absolute abundances**.

The tool works by calculating a sample-specific scaling factor based on the total DNA weight (concentration × volume) present *before* sequencing. This scaling factor is used to *absolutify* the raw counts, allowing for valid simulation of scale across samples.`Absolutifier` estimates uncertainty in these absolute abundances by adding a pseudocount to the observed counts and then drawing Monte Carlo samples from a Dirichlet distribution and provides a measure of error (95% confidence intervals for each feature).

`Absolutifier` works with any type of biological features: **taxa, genes, functional categories, etc.**

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
| `--counts`       | CSV file with **sample rows** and **feature columns** (e.g., `sample_id,feature_1,...`) |
| `--meta`         | CSV file with `sample_id` and `DNA_conc` columns                               |
| `--volume`       | DNA volume (microL) before sequencing (fixed for all samples)                              |
| `--output`       | Output CSV file with absolute abundances                                       |
| `--fastq_folder` | **(Required)** Folder with FASTQ/FASTA files                                     |
| `--extension`    | *(Optional)* File extension (default: `.fastq`)                                |
| `--suffixes`     | *(Optional)* List of suffixes for filtering file names (e.g., `_R1 _R2`)       |
| `--singleton`    | *(Optional)* Additional single FASTQ file(s) to include                        |
| `--error_bars`   | *(Optional)* Flag to calculate 95% confidence intervals with Monte Carlo sampling |
| `--mc_samples`   | *(Optional)* Number of Monte Carlo samples for error bars (default: 1000)      |
| `--alpha`        | *(Optional)* Dirichlet prior for the Bayesian error model (default: 0.5)        |
| `--plot`         | *(Optional)* Flag to generate visualization plots of the results               |
| `--top_features` | *(Optional)* Number of top features to include in plots (default: 20)          |
| `--plot_format`  | *(Optional)* Image format for plots (png, pdf, svg; default: png)               |
| `--figsize`      | *(Optional)* Figure size for plots (width height; default: 12 8)                |

---

## 🧪 Input Format Example

**counts.csv**:

```csv
sample_id,feature_1,feature_2,feature_3
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

## 🛠️ Features

- **Scaling factor**: Calculates scaling factors by comparing DNA input to total sequenced base pairs.
- **Robust Error Propagation**: Computes 95% confidence intervals using a Bayesian model with Monte Carlo sampling.
- **Comprehensive Outputs**: Generates consolidated tables with counts, absolute abundances, confidence intervals, and scaling factors.
- **Publication-Ready Plots**: Creates high-quality bar plots, heatmaps, and confidence interval comparisons.
- **Flexible Input**: Supports paired-end, single-end, and mixed FASTQ/FASTA files.

---

Feel free to open issues or contribute to development!
