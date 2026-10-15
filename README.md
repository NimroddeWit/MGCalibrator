# 🧬 MGCalibrator


**MGCalibrator** is a command-line tool for **metagenomic absolute abundance calibration**.
It converts sequencing-derived coverage information into **absolute abundances (number of copies)**. Two input modes are supported (use one of them, not both):

* **BAM mode**: BAM files as input. Results include **statistical confidence intervals**.
* **SingleM mode**: SingleM output files + the raw FASTQ files as input. The SingleM coverage values are used as depths.

Both modes produce an output file with the same columns.

This enables **cross-sample comparison** by correcting for:

* differences in sequencing depth, and
* differences in the total DNA input before sequencing.

The tool uses per-sample scaling factors and, in BAM mode, Monte Carlo simulation to estimate measurement uncertainty.

---

## 🚀 Overview

### BAM mode

MGCalibrator performs the following steps:

1. **(Optional)** Filters BAM files by minimum read percent identity.
2. **Extracts read depths** per reference sequence from BAM files.
3. **(Optional)** Clusters or bins references according to user-provided CSVs.
4. **Runs Monte Carlo simulations** to estimate confidence intervals for each reference’s mean depth.
5. **Computes sample-specific scaling factors** from input DNA masses and sequenced base pairs.
6. **Calibrates raw depths** to absolute abundances (number of copies), with confidence intervals.

### SingleM mode

1. Reads the SingleM output files (`sample`, `full_coverage` and `taxonomy` columns).
2. Links every SingleM sample name to the sample name derived from the FASTQ file names.
3. Counts the sequenced bases in the raw FASTQ file(s) of each sample.
4. **Computes sample-specific scaling factors** from input DNA masses and sequenced base pairs.
5. **Calibrates the SingleM coverage values** (`full_coverage`, used as `depth`) to absolute abundances (number of copies). No filtering or Monte Carlo simulation is performed; the error columns are `NA`.

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
* [`samtools`](http://www.htslib.org/) (BAM mode only)
* [`coverm`](https://github.com/wwood/CoverM) (BAM mode only)
* Other dependencies are installed automatically (`pandas`, `numpy`, `pysam`, etc.)

---

## ▶️ Example Usage

**BAM mode:**

```bash
mgcalibrator \
    --bam_folder data/input/ \
    --dna_mass data/input/dna_mass.csv \
    --cluster data/input/reference_clusters.csv \
    --bin data/input/reference_bins.csv \
    --output data/output/output.csv \
    --scaling_factors data/output/scaling_factors.txt \
```

**SingleM mode:**

```bash
mgcalibrator \
    --singlem_folder data/input/singlem/ \
    --fastq_folder data/input/fastq/ \
    --dna_mass data/input/dna_mass.csv \
    --output data/output/output.csv \
    --scaling_factors data/output/scaling_factors.txt
```

---

## ⚙️ Command-Line Arguments

| **Argument**               | **Required** | **Description**                                                                                                                      |
|----------------------------|:------------:|--------------------------------------------------------------------------------------------------------------------------------------|
| `--bam_folder`             | ✅ (BAM)      | BAM mode: folder containing input `.bam` files. |
| `--singlem_folder`         | ✅ (SingleM)  | SingleM mode: folder containing the SingleM output files (TSV). Requires `--fastq_folder`. |
| `--fastq_folder`           | ✅ (SingleM)  | SingleM mode: folder containing the raw FASTQ files (plain or gzip-compressed): 1 file per sample for Nanopore reads, 2 files per sample for paired-end Illumina reads. Requires `--singlem_folder`. |
| `--extensions`             | ❌            | BAM mode: file extensions to search for (default: `.sorted.bam`, `.sort.bam`). |
| `--fastq_extensions`       | ❌            | SingleM mode: FASTQ file extensions to search for (default: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`). |
| `--singlem_extensions`     | ❌            | SingleM mode: SingleM output file extensions to search for (default: `.tsv`, `.tsv.gz`). |
| `--suffixes`               | ❌            | Optional list of suffixes to further filter BAM files (BAM mode) or FASTQ files (SingleM mode) by filename, e.g. `_1 _2` to skip singleton FASTQ files. |
| `--sample_map`             | ❌            | SingleM mode: CSV file linking SingleM sample names to FASTQ-derived sample names (columns `singlem_sample` and `sample_id`), for names that cannot be matched automatically. |
| `--cluster`                | ❌           | CSV file mapping reference sequences to clusters (e.g., to stack alleles into one gene).                                             |
| `--bin`                    | ❌           | CSV file mapping reference sequences to bins (e.g., to combine contigs from the same MAG).                                           |
| `--dna_mass`               | ✅           | CSV file containing measured DNA mass (in ng), with columns: `sample_id` and `DNA_mass`.                                             |
| `--scaling_factors`        | ✅           | Path to a file storing precomputed scaling factors (if it exists, it will be reused, otherwise created/updated).                     |
| `--output`                 | ✅           | Path for the output CSV containing absolute abundances (and errors).                                                                 |
| `--threads`                | ❌           | Number of parallel threads/CPU cores to use (default: all available).                                                                |
| `--skip_filtering`         | ❌           | Skip filtering based on mapping identity (assumes BAMs are already filtered).                                                        |
| `--filter_percent`         | ❌           | Minimum read percent identity for filtering (default: `97.0`).                                                                       |
| `--depth_percent`          | ❌           | Percentage of middle depth values to use for depth calculation (default: `50.0`).                                                    |
| `--qubit_error`            | ❌           | Qubit measurement error (as relative error, default: `0.0`).                                                                         |
| `--mc_samples`             | ❌           | Number of Monte Carlo samples for error/confidence interval estimation (default: `100`).                                             |
| `--batch_size`             | ❌           | Number of simulations per batch for Monte Carlo calculation (default: `500`).                                                        |
| `--pseudocount`            | ❌           | Pseudocount for read-mapping probability at zero-depth positions (default: `1.0`).                                                   |
| `--skip_error_calculation` | ❌           | Skip error calculation by Monte Carlo simulation (returns only point estimates).                                                     |

**Note:**
* Use **either** `--bam_folder` **or** `--singlem_folder` + `--fastq_folder`; providing both (or neither) is an error.
* The filtering, clustering/binning, depth, Qubit error and Monte Carlo options only apply to BAM mode. In SingleM mode they are ignored (with a warning).
* The `--qubit_error` option allows you to propagate laboratory measurement error into the uncertainty analysis.
* Error estimation is performed using Monte Carlo simulations, but can be disabled with `--skip_error_calculation` if you only want point estimates.


---

## 🧪 Input File Requirements

### 1. **BAM files** (BAM mode)

* In case of filtering they must be sorted (`.sorted.bam`).
* Filenames should follow this format:

  ```
  sample_1_(...).sorted.bam
  sample_2_(...).sorted.bam
  ```

### 2. **FASTQ files** (SingleM mode)

* Plain or gzip-compressed, standard 4-line FASTQ format.
* One file per sample for Nanopore reads, two files (R1 and R2) per sample for paired-end Illumina reads. The bases of all FASTQ files of a sample are summed.
* The sample name is derived from the file name in the same way as for BAM files (first two underscore-separated parts, ignoring everything after the first dot):

  ```
  sample_1_(...)_1.fastq.gz   ->  sample_1
  sample_1_(...)_2.fastq.gz   ->  sample_1
  sample_2_(...).fastq.gz     ->  sample_2
  ```

### 3. **SingleM output files** (SingleM mode)

* Tab-separated files with at least the columns `sample`, `full_coverage` and `taxonomy`. Other columns (such as `error_up` and `error_down`) are currently ignored.
* A file may contain one or more samples.
* The values in the `sample` column do not need to be identical to the sample names derived from the FASTQ files. MGCalibrator matches them automatically by (in this order): the manual `--sample_map`, an exact match, applying the same naming rule as for FASTQ files (first two underscore-separated parts), or a FASTQ sample name that the SingleM sample name starts with. If a SingleM sample cannot be matched, MGCalibrator stops and lists the unmatched names.
* If multiple SingleM samples are linked to the same FASTQ sample (for example when R1 and R2 were analysed separately), their `full_coverage` values are summed per taxonomy.

### 4. **dna_mass.csv**

| sample_id | DNA_mass |
| --------- | -------- |
| sample_1  | 22.5     |
| sample_2  | 31.1     |

* The `sample_id` column must match the prefixes of the BAM files (BAM mode) or the sample names derived from the FASTQ files (SingleM mode).

### 5. **reference_clusters.csv (optional, BAM mode)**

| sequence   | cluster   |
| ---------- | --------- |
| contig_001 | cluster_1 |
| contig_002 | cluster_1 |

### 6. **reference_bins.csv (optional, BAM mode)**

| sequence   | bin   |
| ---------- | ----- |
| contig_001 | bin_A |
| contig_002 | bin_A |

> ⚠️ Clustered references cannot also be binned.

### 7. **sample_map.csv (optional, SingleM mode)**

| singlem_sample | sample_id |
| -------------- | --------- |
| my_odd_name    | sample_1  |

---

## 📈 Output

The main output is a CSV with the following columns (identical in both modes):

| sample | reference | depth | error_down | error_up | calibrated_depth | calibrated_error_down | calibrated_error_up |
| ------ | --------- | ----- | ---------- | -------- | ---------------- | --------------------- | ------------------- |

In SingleM mode, `reference` is the SingleM `taxonomy`, `depth` is the SingleM `full_coverage`, and all error columns are `NA`.

---

## 🧮 How It Works

### BAM mode

#### Step 1. Filter BAMs

Reads below the percent identity threshold (e.g., 97%) are filtered using **CoverM**.

#### Step 2. Compute Depths

Depths per position are calculated using:

```bash
samtools depth -a file.bam
```

and stored as NumPy arrays per reference.

#### Step 3. Monte Carlo Simulation

Each reference is resampled (`n_simulations` times) to estimate mean coverage (`M50`) and 95% confidence intervals.

#### Step 4. Scaling

Scaling factors are computed as:

<p align="center">
  <strong>Scaling factor</strong> = <sup>Initial DNA mass</sup> / <sub>Final DNA mass after sequencing</sub>
</p>

#### Step 5. Calibration

Final absolute abundances (in number of copies) are derived by multiplying scaled depths with the computed scaling factors.

### SingleM mode

1. The columns `sample`, `full_coverage` and `taxonomy` are read from all SingleM files in `--singlem_folder`.
2. SingleM sample names are linked to the sample names derived from the FASTQ file names.
3. The bases in the FASTQ file(s) of each sample are counted (all files of a sample are summed) and the scaling factor is computed with the same formula as in BAM mode.
4. The SingleM `full_coverage` values (`depth`) are multiplied by the scaling factor.

---

## 💡 Tips

* For reproducibility, use consistent naming between BAM/FASTQ files and DNA mass table.
* In SingleM mode, use FASTQ files that are identical to the reads that were given to SingleM, and use `--suffixes _1 _2` if the FASTQ folder also contains files that should not be counted (such as singleton files).
* Reuse previously computed scaling factors (they are cached in `--scaling_factors`).
* Increase `--mc_samples` for more precise confidence intervals (at the cost of runtime).

---

## 🧰 Dependencies

| Package                 | Purpose             |
| ----------------------- | ------------------- |
| `pandas`, `numpy`       | Data handling       |
| `pysam`                 | BAM file parsing    |
| `coverm`, `samtools`    | External tools (BAM mode only) |

---

## 📄 Citation

If you use **MGCalibrator** in your research, please cite our preprint:

> de Wit, N. et al., 2026. Calibrating for absolute microbiome abundances without spike-ins. BioRxiv. https://doi.org/10.64898/2026.02.26.708180

---