# 🧬 MGCalibrator


**MGCalibrator** is a command-line tool for **metagenomic absolute abundance calibration**.
It converts sequencing-derived coverage information into **absolute abundances (number of copies)**. Three input modes are supported (use exactly one of them):

* **BAM mode**: BAM files as input. Results include **statistical confidence intervals**.
* **SingleM mode**: SingleM output files + the raw FASTQ files as input. The SingleM species-level coverage values are used as depths.
* **Sylph mode**: Sylph results (`.sylphmpa` files from `sylph-tax taxprof`, optionally with the `sylph profile` TSV files) + the raw FASTQ files as input. The Sylph strain-level coverage values (`True_cov`) are used as depths. Confidence intervals can be calculated from the estimated read counts of Sylph. **Sylph must be run with `--estimate-unknown`, and with `--estimate-read-counts` if you want confidence intervals** (see [Preparing Sylph input](#-preparing-sylph-input)).

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
2. Links every SingleM sample name to a `sample_id` of the DNA mass file.
3. Counts the sequenced bases in the raw FASTQ file(s) of each sample.
4. **Computes sample-specific scaling factors** from input DNA masses and sequenced base pairs.
5. **Calibrates the SingleM coverage values** (used as `depth`) to absolute abundances (number of copies). No filtering or Monte Carlo simulation is performed; the error columns are `NA`.

### Sylph mode

1. Reads the `.sylphmpa` files and keeps the strain-level rows (`t__...`) with their coverage.
2. **(Optional)** Reads the `sylph profile` TSV files (run with `--estimate-read-counts`) and takes the estimated number of reads per genome. From this, a **95% confidence interval around the coverage** is calculated (`error_down` / `error_up`): the exact Poisson (Garwood) 95% interval of the read count is calculated and the resulting relative interval is applied to the coverage. The confidence intervals that Sylph reports itself are **not** used. The interval is combined with `--qubit_error` after the calibration.
3. Links every Sylph sample to a `sample_id` of the DNA mass file.
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
* Other dependencies are installed automatically (`pandas`, `numpy`, `scipy`, `pysam`)

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

**Sylph mode** (here with two Sylph databases and confidence intervals; Sylph was run with `--estimate-read-counts`):

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
| `--sylph_folder`             | ✅ (Sylph)      | Sylph mode: one or more folders (or files) with `.sylphmpa` files. Sylph must be run with `--estimate-unknown`. Requires `--fastq_folder`. |
| `--sylph_profile_folder`     | ❌             | Sylph mode: one or more folders (or files) with the `sylph profile` output TSV files, to calculate confidence intervals from the estimated read counts. Requires Sylph to be run with `--estimate-read-counts`. |
| `--fastq_folder`             | ✅ (SingleM/Sylph) | SingleM/Sylph mode: folder containing the raw FASTQ files (plain or gzip-compressed): 1 file per sample for Nanopore reads, 2 files per sample for paired-end Illumina reads. |
| `--extensions`               | ❌             | BAM mode: file extensions to search for (default: `.sorted.bam`, `.sort.bam`).                                                       |
| `--fastq_extensions`         | ❌             | FASTQ file extensions to search for (default: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`).                                               |
| `--singlem_extensions`       | ❌             | SingleM output file extensions to search for (default: `.tsv`, `.tsv.gz`).                                                           |
| `--sylphmpa_extensions`      | ❌             | `.sylphmpa` file extensions to search for in folders (default: `.sylphmpa`, `.sylphmpa.gz`).                                         |
| `--sylph_profile_extensions` | ❌             | Sylph profile file extensions to search for in folders (default: `.tsv`, `.tsv.gz`).                                                 |
| `--suffixes`                 | ❌             | Optional list of suffixes to further filter BAM files (BAM mode) or FASTQ files (SingleM/Sylph mode) by filename, e.g. `_1 _2` to skip singleton FASTQ files. |
| `--sample_map`               | ❌             | CSV file linking input names (BAM/FASTQ file names or SingleM/Sylph sample names) to `sample_id`s of the DNA mass file, for names that cannot be matched automatically. |
| `--cluster`                  | ❌             | BAM mode: CSV file mapping reference sequences to clusters (e.g., to stack alleles into one gene).                                   |
| `--bin`                      | ❌             | BAM mode: CSV file mapping reference sequences to bins (e.g., to combine contigs from the same MAG).                                 |
| `--dna_mass`                 | ✅             | CSV file containing measured DNA mass (in ng), with columns: `sample_id` and `DNA_mass`.                                             |
| `--scaling_factors`          | ✅             | Path to a file storing scaling factors. If it exists, the factors in it are reused (with a warning; the sequence files are not counted again), otherwise it is created/updated. Use `--overwrite_scaling_factors` after changing sequence files or DNA masses. |
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
| `--seed`                     | ❌             | BAM mode: random seed, to make the Monte Carlo simulations reproducible (default: not set). The results are also independent of `--threads`. |
| `--min_depth`                | ❌             | All modes: remove references/taxa with a depth below this value before the calibration, e.g. to leave out noise at very low coverage (default: `0`, keep everything). |
| `--overwrite_scaling_factors`| ❌             | All modes: ignore an existing `--scaling_factors` file and calculate all scaling factors again.                                      |
| `--version`                  | ❌             | Print the version and exit.                                                                                                          |

**Note:**
* Use exactly one input mode: `--bam_folder`, **or** `--singlem_folder` + `--fastq_folder`, **or** `--sylph_folder` + `--fastq_folder`. Providing several modes (or none) is an error.
* Options that do not apply to the selected mode are ignored, with a warning.
* The `--qubit_error` option allows you to propagate laboratory measurement error into the uncertainty analysis. It needs confidence intervals, so it is only used in BAM mode and in Sylph mode with `--sylph_profile_folder` (Sylph run with `--estimate-read-counts`). The relative error of the interval and the Qubit error are combined in quadrature: `sqrt(relative_error² + qubit_error²)`, separately for the lower and upper side.
* In BAM mode, error estimation is performed using Monte Carlo simulations, but can be disabled with `--skip_error_calculation` if you only want point estimates.
* **Reproducibility:** every run writes `<output>_run_info.json` next to the output file, with the command line, all parameters (including defaults and the seed), the versions of Python, numpy, pandas, scipy and pysam (and samtools/CoverM in BAM mode), and a few numbers about the data (number of samples and rows, which scaling factors were reused). Sylph and SingleM mode contain no random steps; in BAM mode, use `--seed` to get identical results.
* **Scaling factors:** the scaling factors depend on the DNA mass and on the number of sequenced bases. If you reuse a `--scaling_factors` file, MGCalibrator warns that the factors are not recalculated, and also warns if the file contains samples that are not in the current run. Use `--overwrite_scaling_factors` (or a new file name) when the sequence files or `dna_mass.csv` have changed.


---

## 🧪 Input File Requirements

### 1. **BAM files** (BAM mode)

* In case of filtering they must be sorted (`.sorted.bam`).
* One BAM file per sample. Each file name must start with a `sample_id` of the DNA mass file, followed by `_`, `.` or `-` (see *Matching sample names* below):

  ```
  sample_1_(...).sorted.bam   ->  sample_1
  ERR1543552.sorted.bam       ->  ERR1543552
  ```

### 2. **FASTQ files** (SingleM and Sylph mode)

* Plain or gzip-compressed, standard 4-line FASTQ format.
* One file per sample for Nanopore reads, two files (R1 and R2) per sample for paired-end Illumina reads. The bases of all FASTQ files of a sample are summed.
* Each file is assigned to a sample by its file name, in the same way as for BAM files (see *Matching sample names* below). Mate numbers (`_1`, `_2`, `_R1`, `_R2`) never have to be interpreted:

  ```
  ERR1543552_1.fastq.gz       ->  ERR1543552   (sample_id: ERR1543552)
  ERR1543552_2.fastq.gz       ->  ERR1543552
  sample_1_(...)_1.fastq.gz   ->  sample_1     (sample_id: sample_1)
  Patient_A_001_R1.fastq.gz   ->  Patient_A_001
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
* **`sylph profile` TSV files** (optional, `--sylph_profile_folder`): the output of `sylph profile` that the `.sylphmpa` files were made from. They are only needed for the confidence intervals, and **only work if Sylph was run with `--estimate-read-counts`** (then the `Sequence_abundance` column holds estimated read counts instead of percentages). These can be single- or multi-sample files (one TSV with a `Sample_file` column is the usual output when profiling several samples at once). Each `.sylphmpa` strain row is linked to its row in the profile by sample name and genome accession (the same way sylph-tax links them), and the coverage values of both files are compared as a consistency check. Rows that cannot be linked have `NA` errors. If the `Sequence_abundance` of a sample looks like percentages (sums to at most 100), MGCalibrator gives a warning and no interval is calculated for that sample; if this is the case for all samples, MGCalibrator stops with an error. Files in the folder that are not `sylph profile` output are skipped with a warning.
* Sample names are matched to the `sample_id`s like for SingleM (see below). The directory part of the Sylph sample name is ignored.

**How the confidence interval is calculated.** The number of reads of a genome is Poisson distributed, and the coverage is directly proportional to the number of reads. For every genome, MGCalibrator calculates the exact (Garwood) 95% confidence interval of the estimated read count `N` and applies the bounds, relative to `N`, to the coverage. Example: a genome with coverage 0.5 and 1,000 reads gets a coverage interval of about 0.47–0.53; with 100,000 reads the interval is about 0.497–0.503. The interval is narrow for genomes with many reads (about ±1.96/√N) and wide for genomes with few reads (for example, 0.21–2.9 times the coverage with 3 reads).

* The interval is calculated directly (with `scipy`), not simulated, so it is deterministic. For roughly 100 reads or more it is practically identical to the percentile interval of simulated Poisson draws; with few reads it is wider and does not collapse to zero.
* The interval is a 95% interval, like the Monte Carlo intervals in BAM mode.

> ⚠️ **What the interval does and does not cover.** The interval only describes the *counting noise* of the reads. For abundant genomes it becomes very small (about 0.6% for 100,000 reads), and it should not be read as the total uncertainty of the coverage. It does not include:
>
> * the uncertainty of the read counts themselves: Sylph does not count reads but estimates them (reads = depth × genome size / average read length) and describes the values as rough;
> * systematic effects in Sylph's coverage estimate, such as differences between the reference genome and the genome in your sample, the reassignment of shared k-mers, the estimated read accuracy (`--read-seq-id`), and the fact that Sylph works with a subsample of the k-mers;
> * the uncertainty of the DNA mass and the base count, except for what you give in `--qubit_error`.
>
> The relative interval and `--qubit_error` are combined in quadrature, which assumes that both describe the same confidence level (the interval is a 95% interval). Define your Qubit error accordingly. If accuracy matters, validate the calibrated abundances with a mock community or spike-in with known copy numbers, or with an independent method such as ddPCR.

### 5. **Matching sample names** (all modes)

The `sample_id` column of `dna_mass.csv` defines the sample names. MGCalibrator does not guess a sample name from the *structure* of a file name; it matches every input name (BAM file, FASTQ file, or sample name in the SingleM/Sylph output) against the known `sample_id`s, in this order:

1. The manual `--sample_map` (see below).
2. The name equals a `sample_id` (directory and extension ignored).
3. The name starts with a `sample_id` followed by `_`, `.` or `-`. If several `sample_id`s match, the longest one wins (`S1_rep_R1.fastq.gz` belongs to `S1_rep`, not `S1`; `S10_R1.fastq.gz` never belongs to `S1`).

`sample_id`s can therefore contain any number of underscores, dots or hyphens. Sample IDs must be unique.

If a SingleM/Sylph sample or a BAM file cannot be matched, MGCalibrator stops and lists the unmatched names. FASTQ files that match no `sample_id` are ignored with a warning, as are `sample_id`s without any results. In BAM mode, several BAM files for one sample give an error (merge them first). Every assignment is logged (`Assigned 'x' -> sample 'y'`), so you can check it.

### 6. **dna_mass.csv**

| sample_id | DNA_mass |
| --------- | -------- |
| sample_1  | 22.5     |
| sample_2  | 31.1     |

* The `sample_id` column must be unique and must be the prefix of the corresponding BAM/FASTQ file names (see *Matching sample names*).

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

### 9. **sample_map.csv (optional, all modes)**

| name         | sample_id |
| ------------ | --------- |
| my_odd_name  | sample_1  |

* `name` can be a BAM/FASTQ file name (with or without extension) or a sample name as written in the SingleM/Sylph output. Only needed for names that do not start with their `sample_id`.
* The first column may also be called `singlem_sample` or `sylph_sample`.

-------------- | --------- |
| my_odd_name    | sample_1  |

The first column may also be called `sylph_sample`.

---

## 🐟 Preparing Sylph input

1. **Run Sylph.** Profiling several samples in one run is normal: all samples end up in one TSV file with a `Sample_file` column.

   ```bash
   sylph profile db1.syldb db2.syldb -1 A_1.fastq.gz B_1.fastq.gz -2 A_2.fastq.gz B_2.fastq.gz \
       --estimate-unknown --estimate-read-counts -o profiles.tsv
   ```

   If you use several databases, give them to the **same** `sylph profile` run where possible, so that Sylph can handle all genomes together. Separate runs per database are also supported, but then the same reads can contribute to genomes of different databases.

   > ⚠️ **Always run Sylph with `--estimate-unknown` (`-u`).** Without it, Sylph reports the effective coverage (`Eff_cov`), which underestimates the true sequencing depth. With `--estimate-unknown`, Sylph reports the true coverage (`True_cov`), which is much more accurate (see the Sylph paper). MGCalibrator assumes that the coverage in the `.sylphmpa` files is the true coverage.
   >
   > ⚠️ **Use the same `--read-seq-id` for all samples of a data set.** With `--estimate-unknown`, Sylph estimates the read accuracy for every sample separately unless you give `--read-seq-id` (for example `99` for Illumina reads). This estimate works poorly for complex metagenomes (soil, ocean) and at low sequencing depth, and it influences the coverage (`True_cov`) and therefore the calibration. A constant `--read-seq-id` makes the samples comparable.
   >
   > ⚠️ **If you want confidence intervals, also run Sylph with `--estimate-read-counts`** (Sylph ≥ 0.9.0; this implies `--estimate-unknown`). Sylph then writes the estimated number of reads per genome to the `Sequence_abundance` column, which MGCalibrator uses to calculate the confidence interval of the coverage (see above). Without this flag, no confidence intervals (and no Qubit error propagation) are possible.

2. **Add taxonomy with sylph-tax.** The `-t` option must match the databases used for Sylph; use multiple identifiers if the run used multiple databases.

   ```bash
   sylph-tax taxprof profiles.tsv -t GTDB_r220 IMGVR_4.1 -o taxprof/
   ```

   This writes one `.sylphmpa` file per sample (named after the sample). If you profiled each database separately, put the `.sylphmpa` files of each database in their own folder (or use a different `-o` prefix): `sylph-tax` refuses to write files with the same name twice.

3. **Run MGCalibrator** on the folder(s) with `.sylphmpa` files, the Sylph profile TSV file(s) (optional, for confidence intervals; requires `--estimate-read-counts`) and the FASTQ files:

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

| sample | reference | depth | error_down | error_up | calibrated_depth | calibrated_error_down | calibrated_error_up | scaling_factor | read_counts | error_note |
| ------ | --------- | ----- | ---------- | -------- | ---------------- | --------------------- | ------------------- | -------------- | ----------- | ---------- |

* `scaling_factor`: the factor that was used to convert `depth` to `calibrated_depth` (the same for all rows of a sample).
* `read_counts`: the number of reads behind the depth. BAM mode: the number of mapped reads of the reference. Sylph mode: the number of reads that Sylph estimated for the genome (with `--sylph_profile_folder`). SingleM mode: `NA`. A small count means a wide and unreliable interval.
* `error_note`: empty if the error columns are filled; otherwise the reason why they are `NA` (or zero): `no_error_estimate_for_singlem`, `error_calculation_skipped` (`--skip_error_calculation`), `no_sylph_profile` (no `--sylph_profile_folder`), `not_in_sylph_profile` (genome or sample not found in the profile), `coverage_differs_from_sylph_profile`, `read_counts_are_percentages` (Sylph was not run with `--estimate-read-counts`), `no_read_count`.
* The columns are the same in all modes. The new columns `scaling_factor`, `read_counts` and `error_note` are added at the end, so the position of the earlier columns is unchanged.

* **BAM mode**: `reference` is the reference sequence (or cluster/bin), `depth` the mean depth with Monte Carlo errors.
* **SingleM mode**: `reference` is the full SingleM taxonomy string of the species, `depth` is the SingleM `coverage`, and all error columns are `NA`.
* **Sylph mode**: `reference` is the full taxonomic string of the strain (including the `t__` genome identifier), `depth` is the Sylph coverage (`True_cov`). The error columns are calculated from the Poisson confidence interval of the coverage if `--sylph_profile_folder` is given (`NA` otherwise or where no read count is available).

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

1. The coverage values (SingleM: species-level `coverage`; Sylph: strain-level `Coverage (if strain-level)`) are read from the SingleM/Sylph files and used as depths. For Sylph, the estimated read counts are read from the profile files if supplied, and a Poisson confidence interval around the coverage is calculated from them.
2. SingleM/Sylph sample names are linked to the `sample_id`s of the DNA mass file.
3. The bases in the FASTQ file(s) of each sample are counted (all files of a sample are summed) and the scaling factor is computed with the same formula as in BAM mode.
4. The depths (and errors) are multiplied by the scaling factor.

---

## 💡 Tips

* For reproducibility, use consistent naming between BAM/FASTQ files and DNA mass table.
* In SingleM and Sylph mode, use FASTQ files that are identical to the reads that were given to SingleM/Sylph, and use `--suffixes _1 _2` if the FASTQ folder also contains files that should not be counted (such as singleton files).
* Scaling factors are cached in `--scaling_factors` and reused in later runs; use `--overwrite_scaling_factors` when sequence files or DNA masses have changed.
* Use `--min_depth` to drop references/taxa with a very low coverage, where the estimates are mostly noise.
* Use `--seed` in BAM mode for reproducible Monte Carlo intervals, and keep the `_run_info.json` file together with your results.
* Increase `--mc_samples` for more precise confidence intervals in BAM mode (at the cost of runtime).

---

## 🛠️ Development

```bash
pip install -e ".[test]"
pytest
```

The tests run automatically on GitHub (Python 3.9–3.12) for every push and pull request, see `.github/workflows/tests.yml`.

---

## 🧰 Dependencies

| Package                 | Purpose             |
| ----------------------- | ------------------- |
| `pandas`, `numpy`       | Data handling       |
| `scipy`                 | Exact Poisson confidence intervals (Sylph mode) |
| `pysam`                 | BAM file parsing    |
| `coverm`, `samtools`    | External tools (BAM mode only) |

---

## 📄 Citation

If you use **MGCalibrator** in your research, please cite our preprint:

> de Wit, N. et al., 2026. Calibrating for absolute microbiome abundances without spike-ins. BioRxiv. https://doi.org/10.64898/2026.02.26.708180

---
