# MGCalibrator

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![bioRxiv](https://img.shields.io/badge/preprint-bioRxiv-b31b1b.svg)](https://doi.org/10.64898/2026.02.26.708180)

MGCalibrator converts metagenomic coverage into **absolute abundances (copies)**, without spike-ins. It uses the measured DNA input mass and the number of sequenced bases to correct for differences in sequencing depth and DNA input between samples, enabling cross-sample and time-series comparisons. Where possible, confidence intervals are reported.

## Input modes

Use exactly one mode per run. All modes write the same output format.

| Mode | Input | Depth used | Confidence intervals |
|---|---|---|---|
| **Sylph** | `.sylphmpa` files (`sylph-tax taxprof`) + FASTQ; optionally `sylph profile` TSVs | Strain-level `True_cov` | Exact Poisson (Garwood) 95% interval from Sylph's estimated read counts |
| **BAM** | Coordinate-sorted BAM files | Mean per-reference depth | Monte Carlo resampling |
| **SingleM** | SingleM profile + FASTQ | Species-level `coverage` | None (`NA`) |

## Installation

```bash
mamba install -c conda-forge -c bioconda mgcalibrator
```

or from source:

```bash
git clone https://github.com/NimroddeWit/MGCalibrator.git
cd MGCalibrator
mamba env create -f environment.yml   # optional
pip install .
```

Requires Python ≥ 3.9. `samtools` and `coverm` are needed for BAM mode only. Sylph, sylph-tax and SingleM are **not** bundled: run them on your reads first.

## Quick start

All modes need `--dna_mass`, `--scaling_factors` and `--output`.

### Sylph

```bash
# 1. Profile (see "Preparing Sylph input" for required flags)
sylph profile db.syldb -1 *_1.fastq.gz -2 *_2.fastq.gz \
    --estimate-unknown --estimate-read-counts --read-seq-id 99.5 -o profiles.tsv
sylph-tax taxprof profiles.tsv -t GTDB_r220 -o taxprof/

# 2. Calibrate
mgcalibrator \
    --sylph_folder taxprof/ \
    --sylph_profile_folder profiles.tsv \
    --fastq_folder fastq/ \
    --dna_mass dna_mass.csv \
    --scaling_factors scaling_factors.txt \
    --output calibrated.csv
```

`--sylph_folder` accepts several folders/files (e.g. one per database). `--sylph_profile_folder` is optional and only needed for confidence intervals.

### BAM

```bash
mgcalibrator \
    --bam_folder bams/ \
    --dna_mass dna_mass.csv \
    --cluster reference_clusters.csv \
    --scaling_factors scaling_factors.txt \
    --output calibrated.csv
```

### SingleM

```bash
mgcalibrator \
    --singlem_folder singlem/ \
    --fastq_folder fastq/ \
    --dna_mass dna_mass.csv \
    --scaling_factors scaling_factors.txt \
    --output calibrated.csv
```

## Preparing Sylph input

- Always run `sylph profile` with **`--estimate-unknown`** (`-u`). This reports `True_cov`, which is more accurate than `Eff_cov`; MGCalibrator assumes `True_cov`.
- For confidence intervals, also use **`--estimate-read-counts`** (Sylph ≥ 0.9.0; implies `-u`).
- Use the **same `--read-seq-id`** for all samples of a dataset (e.g. `99` for Illumina). Otherwise Sylph estimates read accuracy per sample, which is unreliable for complex metagenomes and low depth, and affects the calibration.
- Give all databases to the same `sylph profile` run where possible. Separate runs per database are supported, but the same reads can then contribute to genomes of different databases.
- Run `sylph-tax taxprof` with `-t` matching the databases used. It writes one `.sylphmpa` file per sample; keep outputs of separate runs in separate folders.

## Input files

| File | Mode | Description |
|---|---|---|
| `dna_mass.csv` | all | Columns `sample_id` (unique) and `DNA_mass` (ng) |
| `.sylphmpa` files | Sylph | One per sample; strain-level (`t__`) rows only; coverage from `Coverage (if strain-level)` |
| `sylph profile` TSVs | Sylph (optional) | Single- or multi-sample; requires `--estimate-read-counts` |
| BAM files | BAM | One per sample; `.sorted.bam` if filtering |
| SingleM profile | SingleM | Standard TSV with `sample`, `coverage`, `taxonomy`; species-level rows only |
| FASTQ files | Sylph, SingleM | Plain or gzipped; 1 file/sample (Nanopore) or 2 (paired-end Illumina); bases are summed per sample |
| `reference_clusters.csv`, `reference_bins.csv` | BAM (optional) | Columns `sequence`, `cluster` / `bin`; a reference cannot be both clustered and binned |
| `sample_map.csv` | all (optional) | Columns `name`, `sample_id`; only for names that do not start with their `sample_id` |

**Sample matching.** The `sample_id` column of `dna_mass.csv` defines sample names. Each input name (file name, or sample name in SingleM/Sylph output) is matched in this order: (1) `--sample_map`; (2) exact match; (3) the name starts with a `sample_id` followed by `_`, `.` or `-` (longest match wins, so `S10_R1.fastq.gz` never matches `S1`). Mate suffixes (`_1`, `_R2`) need no interpretation. Unmatched names cause an error that lists them; every assignment is logged.

## Output

CSV with one row per sample and reference (identical columns in all modes):

| Column | Description |
|---|---|
| `sample`, `reference` | Sample ID; reference/cluster/bin, or full taxonomy string (SingleM, Sylph) |
| `depth`, `error_down`, `error_up` | Uncalibrated depth with 95% interval (`NA` if unavailable) |
| `calibrated_depth`, `calibrated_error_down`, `calibrated_error_up` | Absolute abundance (copies) with interval |
| `scaling_factor` | Factor used for the calibration (per sample) |
| `read_counts` | Reads behind the depth (BAM: mapped reads; Sylph: estimated reads; SingleM: `NA`) |
| `error_note` | Reason if error columns are `NA` (e.g. `no_error_estimate_for_singlem`, `error_calculation_skipped`, `no_sylph_profile`, `not_in_sylph_profile`, `read_counts_are_percentages`) |

Each run also writes `<output>_run_info.json` (command line, parameters including defaults and seed, tool versions, data summary).

## Command-line options

**Input** (choose one mode)

| Option | Description |
|---|---|
| `--sylph_folder` | Sylph: folder(s)/file(s) with `.sylphmpa` files. Requires `--fastq_folder` |
| `--sylph_profile_folder` | Sylph: folder(s)/file(s) with `sylph profile` TSVs, for confidence intervals |
| `--bam_folder` | BAM: folder with input BAM files |
| `--singlem_folder` | SingleM: folder with SingleM profiles. Requires `--fastq_folder` |
| `--fastq_folder` | Sylph/SingleM: folder with raw FASTQ files |

**Common**

| Option | Description |
|---|---|
| `--dna_mass` (required) | DNA mass table |
| `--scaling_factors` (required) | Scaling-factor file; reused with a warning if it exists (sequence files are not recounted) |
| `--output` (required) | Output CSV |
| `--overwrite_scaling_factors` | Recompute all scaling factors; use after changing sequence files or DNA masses |
| `--sample_map` | Manual name → `sample_id` map |
| `--suffixes` | Filter BAM/FASTQ files by suffix (e.g. `_1 _2` to skip singletons) |
| `--min_depth` | Drop references/taxa below this depth (default: `0`) |
| `--qubit_error` | Relative DNA-quantification error, combined in quadrature with the interval; used in BAM mode and Sylph mode with `--sylph_profile_folder` (default: `0.0`) |
| `--threads` | CPU cores (default: all) |
| `--extensions`, `--fastq_extensions`, `--singlem_extensions`, `--sylphmpa_extensions`, `--sylph_profile_extensions` | File extensions to search for |
| `--version` | Print version |

**BAM only**

| Option | Description |
|---|---|
| `--cluster`, `--bin` | Cluster / bin tables |
| `--skip_filtering` | Skip identity filtering (BAMs already filtered) |
| `--filter_percent` | Minimum read identity in % (default: `97.0`) |
| `--depth_percent` | Percentage of middle depth values used for mean depth (default: `50.0`) |
| `--mc_samples`, `--batch_size` | Monte Carlo samples (default: `100`) and batch size (default: `500`) |
| `--pseudocount` | Pseudocount for zero-depth positions (default: `1.0`) |
| `--skip_error_calculation` | Point estimates only |
| `--seed` | Seed for reproducible Monte Carlo runs |

Options that do not apply to the selected mode are ignored with a warning.

## Method

1. **Depth.** Sylph/SingleM: coverage from the profile. BAM: per-position depth (`samtools depth -a`) after optional CoverM identity filtering, summarised over the middle `--depth_percent` of positions.
2. **Uncertainty.** BAM: Monte Carlo resampling of positions. Sylph: exact Poisson (Garwood) 95% interval of the estimated read count `N`, applied relative to `N` to the coverage. SingleM: none.
3. **Scaling factor.** Per sample, from the DNA mass and the number of sequenced bases (initial DNA mass / DNA mass after sequencing), as described in the preprint.
4. **Calibration.** Depth (and interval) × scaling factor = absolute abundance in copies.

## Caveats

- Intervals describe **counting noise only** (plus `--qubit_error`, if set). They exclude systematic effects: Sylph's estimated (not counted) reads, reference-genome divergence, shared k-mer reassignment, `--read-seq-id`, and uncertainty in DNA mass and base counts. Qubit error and the interval are assumed to be at the same (95%) confidence level.
- BAM mode works best with non-redundant references with minimal sequence overlap.
- Use FASTQ files identical to those given to Sylph/SingleM.
- Where accuracy matters, validate against a mock community, spike-in, or ddPCR/qPCR.

## Development

```bash
pip install -e ".[test]"
pytest
```

Tests run on GitHub Actions (Python 3.9–3.12).

## Citation

> de Wit, N. et al. (2026). Calibrating for absolute microbiome abundances without spike-ins. *bioRxiv*. https://doi.org/10.64898/2026.02.26.708180

## License and acknowledgements

MIT License (see [LICENSE](LICENSE)). MGCalibrator is a fork of [Fuschi/Absolutifier](https://github.com/Fuschi/Absolutifier). Bugs and feature requests: [GitHub Issues](https://github.com/NimroddeWit/MGCalibrator/issues).