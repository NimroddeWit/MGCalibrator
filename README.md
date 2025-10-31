# 🧬 MGCalibrator


**MGCalibrator** is a command-line tool for **metagenomic absolute abundance calibration**.
It converts sequencing-derived coverage information from BAM files into **absolute abundances (number of copies)** with **statistical confidence intervals**.

This enables **cross-sample comparison** by correcting for:

* differences in sequencing depth, and
* differences in the total DNA input before sequencing.

The tool uses per-sample scaling factors and Monte Carlo simulation to estimate measurement uncertainty.

---

## 🚀 Overview

MGCalibrator performs the following steps:

1. **(Optional)** Filters BAM files by minimum read percent identity using [`CoverM`](https://github.com/wwood/CoverM).
2. **Extracts read depths** per reference sequence from BAM files (`samtools depth`).
3. **(Optional)** Clusters or bins references according to user-provided CSVs.
4. **Runs Monte Carlo simulations** to estimate confidence intervals for each reference’s mean depth (M98).
5. **Computes sample-specific scaling factors** from input DNA masses and sequenced base pairs.
6. **Calibrates raw depths** to absolute abundances (number of copies), with confidence intervals.

---

## 📦 Installation

Clone and install the repository with:

```bash
git clone https://github.com/NimroddeWit/MGCalibrator.git
cd MGCalibrator
pip install .
```

**Requirements:**

* Python ≥ 3.9
* [`samtools`](http://www.htslib.org/)
* [`coverm`](https://github.com/wwood/CoverM) (optional but recommended)
* Other dependencies are installed automatically (`pandas`, `numpy`, `pysam`, etc.)

---

## ▶️ Example Usage

```bash
MGCalibrator \
  --bam_folder data/bams/ \
  --dna_mass data/dna_mass.csv \
  --scaling_factors data/scaling_factors.txt \
  --output results/calibrated_abundances.csv \
  --reference_clusters data/reference_clusters.csv \
  --reference_bins data/reference_bins.csv \
  --perc_ident 97 \
  --threads 8 \
  --mc_samples 500 \
  --batch_size 500 \
  --pseudocount 1.0
```

---

## ⚙️ Command-Line Arguments

| **Argument**           | **Required** | **Description**                                                                    |
| ---------------------- | ------------ | ---------------------------------------------------------------------------------- |
| `--bam_folder`         | ✅            | Folder containing input `.bam` files.                                              |
| `--extensions`         | ❌            | File extensions to search for (default: `.sorted.bam`, `.sort.bam`).               |
| `--suffixes`           | ❌            | List of suffixes to filter file names.                                             |
| `--reference_clusters` | ❌            | CSV mapping reference sequences to cluster names.                                  |
| `--reference_bins`     | ❌            | CSV mapping reference sequences to bins.                                           |
| `--dna_mass`           | ✅            | CSV file with measured total DNA mass (ng) per sample (`sample_id, DNA_mass`).     |
| `--scaling_factors`    | ✅            | Path to a text file containing scaling factors (created or updated automatically). |
| `--output`             | ✅            | Path for the output CSV containing absolute abundances.                            |
| `--threads`            | ❌            | Number of CPU cores to use (default: all available).                               |
| `--skip_filtering`     | ❌            | Skip the CoverM filtering step (requires already filtered BAMs).                   |
| `--perc_ident`         | ❌            | Minimum read percent identity for filtering (default: 97).                         |
| `--mc_samples`         | ❌            | Number of Monte Carlo samples for confidence intervals (default: 100).             |
| `--batch_size`         | ❌            | Number of simulations per batch (default: 500).                                    |
| `--pseudocount`        | ❌            | Pseudocount for read-mapping probability (default: 1.0).                           |

---

## 🧪 Input File Requirements

### 1. **BAM files**

* In case of filtering they must be sorted (`.sorted.bam`).
* Filenames should follow this format:

  ```
  sample_1_(...).sorted.bam
  sample_2_(...).sorted.bam
  ```

### 2. **dna_mass.csv**

| sample_id | DNA_mass |
| --------- | -------- |
| sample_1  | 22.5     |
| sample_2  | 31.1     |

* The `sample_id` column must match the prefixes of the BAM files.

### 3. **reference_clusters.csv (optional)**

| sequence   | cluster   |
| ---------- | --------- |
| contig_001 | cluster_1 |
| contig_002 | cluster_1 |

### 4. **reference_bins.csv (optional)**

| sequence   | bin   |
| ---------- | ----- |
| contig_001 | bin_A |
| contig_002 | bin_A |

> ⚠️ Clustered references cannot also be binned.

---

## 📈 Output

The main output is a CSV with the following columns:

| Sample | Reference | M98_mean | M98_lower_ci | M98_upper_ci | scaling_factor | calibrated_depth | calibrated_lower_ci | calibrated_upper_ci |
| ------ | --------- | -------- | ------------ | ------------ | -------------- | ---------------- | ------------------- | ------------------- |

---

## 🧮 How It Works

### Step 1. Filter BAMs

Reads below the percent identity threshold (e.g., 97%) are filtered using **CoverM**.

### Step 2. Compute Depths

Depths per position are calculated using:

```bash
samtools depth -a file.bam
```

and stored as NumPy arrays per reference.

### Step 3. Monte Carlo Simulation

Each reference is resampled (`n_simulations` times) to estimate mean coverage (`M98`) and 95% confidence intervals.

### Step 4. Scaling

Scaling factors are computed as:
[
\text{Scaling factor} = \frac{\text{Initial DNA mass}}{\text{Final DNA mass after sequencing}}
]

### Step 5. Calibration

Final absolute abundances (in number of copies) are derived by multiplying scaled depths with the computed scaling factors.

---

## 💡 Tips

* For reproducibility, use consistent naming between BAM files and DNA mass table.
* Reuse previously computed scaling factors (they are cached in `--scaling_factors`).
* Increase `--mc_samples` for more precise confidence intervals (at the cost of runtime).

---

## 🧰 Dependencies

| Package                 | Purpose             |
| ----------------------- | ------------------- |
| `pandas`, `numpy`       | Data handling       |
| `pysam`                 | BAM file parsing    |
| `matplotlib`, `seaborn` | (Optional) plotting |
| `coverm`, `samtools`    | External tools      |

---

## 📄 Citation

If you use **MGCalibrator** in your research, please cite:

> de Wit, N. et al. (2026). [insert publication title here...]
> GitHub: [https://github.com/NimroddeWit/MGCalibrator](https://github.com/NimroddeWit/MGCalibrator)

---

## 🙌 Acknowledgements

Developed by **Nimrod de Wit**
📍 RIVM, 2025

---
