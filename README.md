# 🧬 MGCalibrator


**MGCalibrator** is a command-line tool for **metagenomic absolute abundance calibration**.
It converts sequencing-derived coverage information into **absolute abundances (number of copies)**. Three input modes are supported (use exactly one of them):

* **BAM mode**: BAM files as input. Results include **statistical confidence intervals**.
* **SingleM mode**: SingleM output files + the raw FASTQ files as input. The SingleM species-level coverage values are used as depths.
* **Sylph mode**: Sylph results (`.sylphmpa` files from `sylph-tax taxprof`, optionally with the `sylph profile` TSV files) + the raw FASTQ files as input. The Sylph strain-level coverage values are used as depths. Sylph's confidence intervals can be included.

All modes produce an output file with the same columns.

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

1. Reads the SingleM output files (`sample`, `coverage` and `taxonomy` columns) and keeps the species-level rows.
2. Links every SingleM sample name to the sample name derived from the FASTQ file names.
3. Counts the sequenced bases in the raw FASTQ file(s) of each sample.
4. **Computes sample-specific scaling factors** from input DNA masses and sequenced base pairs.
5. **Calibrates the SingleM coverage values** (used as `depth`) to absolute abundances (number of copies). No filtering or Monte Carlo simulation is performed; the error columns are `NA`.

### Sylph mode

1. Reads the `.sylphmpa` files and keeps the strain-level rows (`t__...`) with their coverage.
2. **(Optional)** Reads the `sylph profile` TSV files and adds the 5th-95th percentile interval of the coverage as `error_down` / `error_up`.
3. Links every Sylph sample to the sample name derived from the FASTQ file names.
4. Counts the sequenced bases in the raw FASTQ file(s), **computes sample-specific scaling factors** and **calibrates the Sylph coverage values** (used as `depth`) to absolute abundances (number of copies).

---

## 📦 Installation

### With mamba / conda (Bioconda)

```bash
mamba install -c conda-forge -c bioconda mgcalibrator
```

This installs MGCalibrator together with all its dependencies, including `samtools` and `coverm` for BAM mode. (Available once the recipe has been added to Bioconda; see [`packaging/bioconda`](packaging/bioconda).)

### From source

Clone and install the repository with:

```bash
git clone https://github.com/NimroddeWit/MGCalibrator.git
cd MGCalibrator
mamba env create -f environment.yml   # optional: creates an environment with all dependencies
pip install .
```

**Requirements:**

* Python ≥ 3.9
* [`samtools`](http://www.htslib.org/) (BAM mode only)
* [`coverm`](https://github.com/wwood/CoverM) (BAM mode only)
* Other dependencies are installed automatically (`pandas`, `numpy`, `pysam`)

SingleM, Sylph and sylph-tax are not part of MGCalibrator: run them on your FASTQ files first (see [Preparing Sylph input](#-preparing-sylph-input)).

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
    --scaling_factors data/output/scaling_factors.txt
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

**Sylph mode** (here with two Sylph databases and confidence intervals):

```bash
mgcalibrator \
    --sylph_folder data/input/sylph/gtdb data/input/sylph/virus \
    --sylph_profile_folder data/input/sylph_profiles/ \
    --fastq_folder data/input/fastq/ \
    --dna_mass data/input/dna_mass.csv \
    --output data/output/output.csv \
    --scaling_factors data/output/scaling_factors.txt
```

---

## ⚙️ Command-Line Arguments

| **Argument**                 | **Required**   | **Description**                                                                                                                      |
|------------------------------|:--------------:|--------------------------------------------------------------------------------------------------------------------------------------|
| `--bam_folder`               | ✅ (BAM)        | BAM mode: folder containing input `.bam` files.                                                                                      |
| `--singlem_folder`           | ✅ (SingleM)    | SingleM mode: folder containing the SingleM output files (standard taxonomic profile TSV). Requires `--fastq_folder`.                |
| `--sylph_folder`             | ✅ (Sylph)      | Sylph mode: one or more folders (or files) with `.sylphmpa` files. Requires `--fastq_folder`.                                        |
| `--sylph_profile_folder`     | ❌             | Sylph mode: one or more folders (or files) with the `sylph profile` output TSV files, to add confidence intervals.                    |
| `--fastq_folder`             | ✅ (SingleM/Sylph) | SingleM/Sylph mode: folder containing the raw FASTQ files (plain or gzip-compressed): 1 file per sample for Nanopore reads, 2 files per sample for paired-end Illumina reads. |
| `--extensions`               | ❌             | BAM mode: file extensions to search for (default: `.sorted.bam`, `.sort.bam`).                                                       |
| `--fastq_extensions`         | ❌             | FASTQ file extensions to search for (default: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`).                                               |
| `--singlem_extensions`       | ❌             | SingleM output file extensions to search for (default: `.tsv`, `.tsv.gz`).                                                           |
| `--sylphmpa_extensions`      | ❌             | `.sylphmpa` file extensions to search for in folders (default: `.sylphmpa`, `.sylphmpa.gz`).                                         |
| `--sylph_profile_extensions` | ❌             | Sylph profile file extensions to search for in folders (default: `.tsv`, `.tsv.gz`).                                                 |
| `--suffixes`                 | ❌             | Optional list of suffixes to further filter BAM files (BAM mode) or FASTQ files (SingleM/Sylph mode) by filename, e.g. `_1 _2` to skip singleton FASTQ files. |
| `--sample_map`               | ❌             | SingleM/Sylph mode: CSV file linking sample names in the SingleM/Sylph output to FASTQ-derived sample names, for names that cannot be matched automatically. |
| `--cluster`                  | ❌             | BAM mode: CSV file mapping reference sequences to clusters (e.g., to stack alleles into one gene).                                   |
| `--bin`                      | ❌             | BAM mode: CSV file mapping reference sequences to bins (e.g., to combine contigs from the same MAG).                                 |
| `--dna_mass`                 | ✅             | CSV file containing measured DNA mass (in ng), with columns: `sample_id` and `DNA_mass`.                                             |
| `--scaling_factors`          | ✅             | Path to a file storing precomputed scaling factors (if it exists, it will be reused, otherwise created/updated).                     |
| `--output`                   | ✅             | Path for the output CSV containing absolute abundances (and errors).                                                                 |
| `--threads`                  | ❌             | Number of parallel threads/CPU cores to use (default: all available).                                                                |
| `--skip_filtering`           | ❌             | BAM mode: skip filtering based on mapping identity (assumes BAMs are already filtered).                                              |
| `--filter_percent`           | ❌             | BAM mode: minimum read percent identity for filtering (default: `97.0`).                                                             |
| `--depth_percent`            | ❌             | BAM mode: percentage of middle depth values to use for depth calculation (default: `50.0`).                                          |
| `--qubit_error`              | ❌             | Qubit measurement error (as relative error, default: `0.0`). BAM mode and Sylph mode with `--sylph_profile_folder`.                  |
| `--mc_samples`               | ❌             | BAM mode: number of Monte Carlo samples for error/confidence interval estimation (default: `100`).                                   |
| `--batch_size`               | ❌             | BAM mode: number of simulations per batch for Monte Carlo calculation (default: `500`).                                              |
| `--pseudocount`              | ❌             | BAM mode: pseudocount for read-mapping probability at zero-depth positions (default: `1.0`).                                         |
| `--skip_error_calculation`   | ❌             | BAM mode: skip error calculation by Monte Carlo simulation (returns only point estimates).                                           |
| `--version`                  | ❌             | Print the version and exit.                                                                                                          |

**Note:**
* Use exactly one input mode: `--bam_folder`, **or** `--singlem_folder` + `--fastq_folder`, **or** `--sylph_folder` + `--fastq_folder`. Providing several modes (or none) is an error.
* Options that do not apply to the selected mode are ignored, with a warning.
* The `--qubit_error` option allows you to propagate laboratory measurement error into the uncertainty analysis. It needs confidence intervals, so it is only used in BAM mode and in Sylph mode with `--sylph_profile_folder`.
* In BAM mode, error estimation is performed using Monte Carlo simulations, but can be disabled with `--skip_error_calculation` if you only want point estimates.


---

## 🧪 Input File Requirements

### 1. **BAM files** (BAM mode)

* In case of filtering they must be sorted (`.sorted.bam`).
* Filenames should follow this format:

  ```
  sample_1_(...).sorted.bam
  sample_2_(...).sorted.bam
  ```

### 2. **FASTQ files** (SingleM and Sylph mode)

* Plain or gzip-compressed, standard 4-line FASTQ format.
* One file per sample for Nanopore reads, two files (R1 and R2) per sample for paired-end Illumina reads. The bases of all FASTQ files of a sample are summed.
* The sample name is derived from the file name in the same way as for BAM files (first two underscore-separated parts, ignoring everything after the first dot):

  ```
  sample_1_(...)_1.fastq.gz   ->  sample_1
  sample_1_(...)_2.fastq.gz   ->  sample_1
  sample_2_(...).fastq.gz     ->  sample_2
  ```

### 3. **SingleM output files** (SingleM mode)

* The standard SingleM taxonomic profile: a tab-separated file with (at least) the columns `sample`, `coverage` and `taxonomy`, as written by `singlem pipe` (`--taxonomic-profile`) or `singlem condense`. Other columns are ignored.
* Only **species-level** rows are used (taxonomy ending in `s__<name>`): at species level the coverage is equal to the full coverage, because a species has no descendants. Rows at higher taxonomic levels are dropped.
* The complete taxonomy string is stored in the `reference` column of the output.
* A file may contain one or more samples.
* If multiple SingleM samples are linked to the same FASTQ sample (for example when R1 and R2 were analysed separately), their coverage values are summed per taxonomy.

### 4. **Sylph results** (Sylph mode)

* **`.sylphmpa` files** (required), the output of `sylph-tax taxprof`. `sylph-tax taxprof` writes one `.sylphmpa` file per sample; the sample name (the `Sample_file` of Sylph) is stored in the `#SampleID` header line of each file, so the file names do not matter. There is no merged `.sylphmpa` table with coverage (`sylph-tax merge` only supports ANI and abundances), so give the per-sample files.
* Only **strain-level** rows (`t__...`) are used, because Sylph only reports coverage at that level. The coverage is taken from the column `Coverage (if strain-level)` (`Eff_cov`, or `True_cov` if Sylph was run with `-u`). The complete taxonomic string (with `|` separators, as written by sylph-tax) is stored in the `reference` column of the output.
* `--sylph_folder` accepts **multiple folders and/or files**. Files for the same sample (for example from separate Sylph runs with different databases) are combined, based on the sample name in the file header. If the same genome occurs twice for a sample, both rows are kept and a warning is given.
* **`sylph profile` TSV files** (optional, `--sylph_profile_folder`): the output of `sylph profile` that the `.sylphmpa` files were made from. These can be single- or multi-sample files (one TSV with a `Sample_file` column is the usual output when profiling several samples at once). The `Lambda_5-95_percentile` column is used as `error_down` / `error_up`. Each `.sylphmpa` strain row is linked to its row in the profile by sample name and genome accession (the same way sylph-tax links them), and the coverage values of both files are compared as a consistency check. Rows for which Sylph could not estimate an interval (coverage too low or too high: `Eff_lambda` is `LOW` / `HIGH`) have `NA` errors. Files in the folder that are not `sylph profile` output are skipped with a warning.
* Sample names are matched to the FASTQ files like for SingleM (see below). The directory part of the Sylph sample name is ignored.

> ℹ️ Sylph reports a **5th–95th percentile** interval (a 90% interval), whereas the Monte Carlo intervals of BAM mode are 95% intervals. The intervals of the two modes are therefore not directly comparable.

### 5. **Matching sample names** (SingleM and Sylph mode)

The sample names in the SingleM/Sylph output do not need to be identical to the sample names derived from the FASTQ files. MGCalibrator matches them automatically by (in this order): the manual `--sample_map`, an exact match, applying the same naming rule as for FASTQ files (first two underscore-separated parts), or a FASTQ sample name that the SingleM/Sylph sample name starts with. If a sample cannot be matched, MGCalibrator stops and lists the unmatched names.

### 6. **dna_mass.csv**

| sample_id | DNA_mass |
| --------- | -------- |
| sample_1  | 22.5     |
| sample_2  | 31.1     |

* The `sample_id` column must match the prefixes of the BAM files (BAM mode) or the sample names derived from the FASTQ files (SingleM and Sylph mode).

### 7. **reference_clusters.csv (optional, BAM mode)**

| sequence   | cluster   |
| ---------- | --------- |
| contig_001 | cluster_1 |
| contig_002 | cluster_1 |

### 8. **reference_bins.csv (optional, BAM mode)**

| sequence   | bin   |
| ---------- | ----- |
| contig_001 | bin_A |
| contig_002 | bin_A |

> ⚠️ Clustered references cannot also be binned.

### 9. **sample_map.csv (optional, SingleM and Sylph mode)**

| singlem_sample | sample_id |
| -------------- | --------- |
| my_odd_name    | sample_1  |

The first column may also be called `sylph_sample`.

---

## 🐟 Preparing Sylph input

1. **Run Sylph.** Profiling several samples in one run is normal: all samples end up in one TSV file with a `Sample_file` column.

   ```bash
   sylph profile db1.syldb db2.syldb -1 A_1.fastq.gz B_1.fastq.gz -2 A_2.fastq.gz B_2.fastq.gz -o profiles.tsv
   ```

   If you use several databases, give them to the **same** `sylph profile` run where possible, so that Sylph can handle all genomes together. Separate runs per database are also supported, but then the same reads can contribute to genomes of different databases.

   Sylph's `Eff_cov` slightly underestimates the true sequencing depth. With `sylph profile -u`, Sylph reports `True_cov` instead; MGCalibrator supports both (the confidence intervals are rescaled automatically).

2. **Add taxonomy with sylph-tax.** The `-t` option must match the databases used for Sylph; use multiple identifiers if the run used multiple databases.

   ```bash
   sylph-tax taxprof profiles.tsv -t GTDB_r220 IMGVR_4.1 -o taxprof/
   ```

   This writes one `.sylphmpa` file per sample (named after the sample). If you profiled each database separately, put the `.sylphmpa` files of each database in their own folder (or use a different `-o` prefix): `sylph-tax` refuses to write files with the same name twice.

3. **Run MGCalibrator** on the folder(s) with `.sylphmpa` files, the Sylph profile TSV file(s) (optional, for confidence intervals) and the FASTQ files:

   ```bash
   mgcalibrator \
       --sylph_folder taxprof_db1/ taxprof_db2/ \
       --sylph_profile_folder profiles_db1.tsv profiles_db2.tsv \
       --fastq_folder fastq/ \
       --dna_mass dna_mass.csv \
       --scaling_factors scaling_factors.txt \
       --output calibrated.csv
   ```

   Use the FASTQ files that were given to Sylph (and `--suffixes _1 _2` if the folder also contains files that were not profiled, such as singleton files).

---

## 📈 Output

The main output is a CSV with the following columns (identical in all modes):

| sample | reference | depth | error_down | error_up | calibrated_depth | calibrated_error_down | calibrated_error_up |
| ------ | --------- | ----- | ---------- | -------- | ---------------- | --------------------- | ------------------- |

* **BAM mode**: `reference` is the reference sequence (or cluster/bin), `depth` the mean depth with Monte Carlo errors.
* **SingleM mode**: `reference` is the full SingleM taxonomy string of the species, `depth` is the SingleM `coverage`, and all error columns are `NA`.
* **Sylph mode**: `reference` is the full taxonomic string of the strain (including the `t__` genome identifier), `depth` is the Sylph coverage. The error columns are filled from the Sylph confidence intervals if `--sylph_profile_folder` is given (`NA` otherwise or where no interval is available).

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

### SingleM and Sylph mode

1. The coverage values (SingleM: species-level `coverage`; Sylph: strain-level `Coverage (if strain-level)`) are read from the SingleM/Sylph files and used as depths. For Sylph, the confidence intervals are read from the profile files if supplied.
2. SingleM/Sylph sample names are linked to the sample names derived from the FASTQ file names.
3. The bases in the FASTQ file(s) of each sample are counted (all files of a sample are summed) and the scaling factor is computed with the same formula as in BAM mode.
4. The depths (and errors) are multiplied by the scaling factor.

---

## 💡 Tips

* For reproducibility, use consistent naming between BAM/FASTQ files and DNA mass table.
* In SingleM and Sylph mode, use FASTQ files that are identical to the reads that were given to SingleM/Sylph, and use `--suffixes _1 _2` if the FASTQ folder also contains files that should not be counted (such as singleton files).
* Reuse previously computed scaling factors (they are cached in `--scaling_factors`).
* Increase `--mc_samples` for more precise confidence intervals in BAM mode (at the cost of runtime).

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
