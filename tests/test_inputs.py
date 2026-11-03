"""Tests for reading SingleM / Sylph input and matching sample names."""

import numpy as np
import pandas as pd
import pytest

from MGCalibrator.fileutils import collect_files, strip_extension
from MGCalibrator.processor import filtered_bam_path
from MGCalibrator.samples import (
    assign_files_to_samples,
    build_sample_map,
    check_sample_ids,
    group_by_sample,
    match_sample_id,
)
from MGCalibrator.singlem import read_singlem_tables, singlem_to_depth_table
from MGCalibrator.sylph import (
    add_poisson_confidence_intervals,
    add_read_counts,
    poisson_coverage_ci,
    read_sylph_profiles,
    read_sylphmpa_files,
    sylph_to_depth_table,
)

PROFILE_HEADER = (
    "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t{cov}\t"
    "ANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\t"
    "Containment_ind\tNaive_ANI\tkmers_reassigned\tContig_name\n"
)


def write_sylphmpa(path, sample, rows):
    lines = [f"#SampleID\t{sample}\tTaxonomies_used:['x']\n",
             "clade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)\n",
             "d__Bacteria\t100.0\t100.0\tNA\tNA\n",
             "d__Bacteria|g__Escherichia|s__Escherichia coli\t100.0\t100.0\tNA\tNA\n"]
    lines += [f"{clade}\t100.0\t100.0\t99.0\t{cov}\n" for clade, cov in rows]
    path.write_text("".join(lines))


def profile_row(sample, genome, cov, reads, contig="NC_1 Escherichia coli"):
    # The Sylph confidence interval columns (Eff_lambda, Lambda_5-95_percentile) are ignored by MGCalibrator
    return (f"{sample}\t{genome}\t50\t{reads}\t99.0\t{cov}\tNA-NA\t{cov}\t0.01-9.99\t1\t1.0\t10/20\t99.0\t0\t{contig}\n")


# ----------------------------------------------------------------------
# SingleM
# ----------------------------------------------------------------------
def test_singlem_keeps_species_only_and_full_taxonomy(tmp_path):
    f = tmp_path / "s.tsv"
    f.write_text(
        "sample\tcoverage\ttaxonomy\n"
        "a_1\t0.5\tRoot; d__Bacteria\n"
        "a_1\t0.8\tRoot; d__Bacteria; g__Bacillus\n"
        "a_1\t2.0\tRoot; d__Bacteria; g__Escherichia; s__Escherichia coli\n"
        "a_1\t1.0\tRoot; d__Bacteria; g__Escherichia; s__\n"
    )
    df = read_singlem_tables([str(f)])
    assert list(df["taxonomy"]) == ["Root; d__Bacteria; g__Escherichia; s__Escherichia coli"]
    depth = singlem_to_depth_table(df, {"a_1": "a_1"})
    assert depth.loc[0, "depth"] == 2.0
    assert depth.loc[0, "reference"] == "Root; d__Bacteria; g__Escherichia; s__Escherichia coli"


def test_singlem_requires_coverage_column(tmp_path):
    f = tmp_path / "s.tsv"
    f.write_text("sample\tfull_coverage\ttaxonomy\na\t1\tRoot; s__x y\n")
    with pytest.raises(ValueError, match="coverage"):
        read_singlem_tables([str(f)])


# ----------------------------------------------------------------------
# Sylph: .sylphmpa
# ----------------------------------------------------------------------
def test_sylphmpa_strain_rows_only_and_sample_from_header(tmp_path):
    f = tmp_path / "whatever_name.sylphmpa"
    write_sylphmpa(f, "/data/reads/s1_x.fastq.gz",
                   [("d__Bacteria|g__Escherichia|s__Escherichia coli|t__GCF_1.1", "0.5")])
    df = read_sylphmpa_files([str(f)])
    assert len(df) == 1
    assert df.loc[0, "sample"] == "/data/reads/s1_x.fastq.gz"
    assert df.loc[0, "reference"] == "d__Bacteria|g__Escherichia|s__Escherichia coli|t__GCF_1.1"
    assert df.loc[0, "depth"] == 0.5
    assert df.loc[0, "accession"] == "GCF_1.1"


def test_sylphmpa_files_of_same_sample_are_combined(tmp_path):
    (tmp_path / "db1").mkdir()
    (tmp_path / "db2").mkdir()
    write_sylphmpa(tmp_path / "db1" / "s.sylphmpa", "s1_x.fq.gz", [("d__A|t__G1", "1.0")])
    write_sylphmpa(tmp_path / "db2" / "s.sylphmpa", "s1_x.fq.gz", [("r__V|t__V1", "2.0")])
    files = collect_files([str(tmp_path / "db1"), str(tmp_path / "db2")], [".sylphmpa"])
    df = read_sylphmpa_files(files)
    assert set(df["sample"]) == {"s1_x.fq.gz"}
    assert sorted(df["depth"]) == [1.0, 2.0]


def test_pavian_sylphmpa_is_rejected(tmp_path):
    f = tmp_path / "p.sylphmpa"
    f.write_text("#mpa_v3_sylphmock#SampleID\ts\tTaxonomies_used:['x']\n#clade_name\tplaceholder\trelative_abundance\tplaceholder2\nd__A\t0\t1\t\n")
    with pytest.raises(ValueError, match="coverage"):
        read_sylphmpa_files([str(f)])


# ----------------------------------------------------------------------
# Sylph: profile TSV -> read counts -> Poisson confidence intervals
# ----------------------------------------------------------------------
def make_strains(tmp_path, sample="s1_x.fq.gz"):
    f = tmp_path / "s.sylphmpa"
    write_sylphmpa(f, sample, [
        ("d__A|s__E coli|t__GCF_000005845.2", "0.412"),
        ("d__A|s__B subtilis|t__GCF_000009045.1", "12.35"),
    ])
    return read_sylphmpa_files([str(f)])


ECOLI = "GCF_000005845.2_ASM584v2_genomic.fna.gz"
BSUB = "GCF_000009045.1_ASM904v1_genomic.fna.gz"


def test_read_counts_joined_on_sample_and_accession(tmp_path):
    p = tmp_path / "profile.tsv"
    p.write_text(
        PROFILE_HEADER.format(cov="True_cov")
        + profile_row("s1_x.fq.gz", ECOLI, "0.412", "12700")
        + profile_row("s1_x.fq.gz", BSUB, "12.350", "346000")
        + profile_row("other.fq.gz", ECOLI, "0.412", "999")
    )
    result = add_read_counts(make_strains(tmp_path), read_sylph_profiles([str(p)]))
    assert list(result["read_counts"]) == [12700, 346000]


def test_sylph_confidence_intervals_are_not_used(tmp_path):
    # The profile has a very wide Sylph interval (0.01-9.99); the result must only depend on the read count
    p = tmp_path / "profile.tsv"
    p.write_text(PROFILE_HEADER.format(cov="True_cov") + profile_row("s1_x.fq.gz", ECOLI, "0.412", "1000000"))
    strains = add_read_counts(make_strains(tmp_path), read_sylph_profiles([str(p)]))
    result = add_poisson_confidence_intervals(strains)
    assert result.loc[0, "lower_ci"] == pytest.approx(0.412 * 0.9980, abs=0.0005)
    assert result.loc[0, "upper_ci"] == pytest.approx(0.412 * 1.0020, abs=0.0005)
    assert result.loc[0, "error_note"] == ""
    # The second genome is not in the profile: no interval, with the reason
    assert np.isnan(result.loc[1, "lower_ci"]) and result.loc[1, "error_note"] == "not_in_sylph_profile"


def test_inconsistent_coverage_gives_no_read_counts(tmp_path):
    p = tmp_path / "profile.tsv"
    p.write_text(PROFILE_HEADER.format(cov="True_cov") + profile_row("s1_x.fq.gz", ECOLI, "0.900", "12700"))
    result = add_read_counts(make_strains(tmp_path), read_sylph_profiles([str(p)]))
    assert result["read_counts"].isna().all()


def test_percentages_instead_of_read_counts_are_rejected(tmp_path):
    # Sylph run without --estimate-read-counts: Sequence_abundance is a percentage (sums to 100)
    p = tmp_path / "profile.tsv"
    p.write_text(
        PROFILE_HEADER.format(cov="True_cov")
        + profile_row("s1_x.fq.gz", ECOLI, "0.412", "60.0")
        + profile_row("s1_x.fq.gz", BSUB, "12.350", "40.0")
    )
    with pytest.raises(ValueError, match="--estimate-read-counts"):
        read_sylph_profiles([str(p)])


def test_percentage_samples_among_read_count_samples_get_no_interval(tmp_path, caplog):
    p = tmp_path / "profile.tsv"
    p.write_text(
        PROFILE_HEADER.format(cov="True_cov")
        + profile_row("s1_x.fq.gz", ECOLI, "0.412", "12700")
        + profile_row("pct.fq.gz", ECOLI, "0.412", "100.0")
    )
    profiles = read_sylph_profiles([str(p)])
    assert profiles.set_index("sample")["read_counts"].isna().to_dict() == {"s1_x.fq.gz": False, "pct.fq.gz": True}
    assert "pct.fq.gz" in caplog.text


def test_effective_coverage_gives_warning(tmp_path, caplog):
    p = tmp_path / "profile.tsv"
    p.write_text(PROFILE_HEADER.format(cov="Eff_cov") + profile_row("s", "GCF_1.1_ASM1v1_genomic.fna.gz", "1.000", "5000"))
    read_sylph_profiles([str(p)])
    assert "--estimate-unknown" in caplog.text


def test_non_profile_tsv_is_skipped(tmp_path):
    good = tmp_path / "profile.tsv"
    good.write_text(PROFILE_HEADER.format(cov="True_cov")
                    + profile_row("s", "GCF_1.1_ASM1v1_genomic.fna.gz", "1.000", "5000"))
    other = tmp_path / "merged.tsv"
    other.write_text("clade_name\ts1\nd__A\t1\n")
    assert len(read_sylph_profiles([str(good), str(other)])) == 1
    with pytest.raises(ValueError):
        read_sylph_profiles([str(other)])


def test_poisson_ci_matches_the_documented_examples():
    # Same coverage, 1,000 vs 100,000 reads: relative interval about +-6% vs +-0.6%
    lower, upper = poisson_coverage_ci([0.5, 0.5], [1000, 100000])
    assert (lower[0], upper[0]) == pytest.approx((0.47, 0.53), abs=0.005)
    assert (lower[1], upper[1]) == pytest.approx((0.497, 0.503), abs=0.001)


def test_poisson_ci_is_the_exact_garwood_interval():
    from scipy.stats import chi2

    reads = np.array([1.0, 3.0, 20.0, 1000.0, 123456.0, 12.5])  # also non-integer (estimated) counts
    lower, upper = poisson_coverage_ci(np.ones_like(reads), reads)
    np.testing.assert_allclose(lower, chi2.ppf(0.025, 2 * reads) / 2 / reads, rtol=1e-9)
    np.testing.assert_allclose(upper, chi2.ppf(0.975, 2 * reads + 2) / 2 / reads, rtol=1e-9)
    # Few reads: wide and asymmetric interval, but never below zero
    assert lower[1] == pytest.approx(0.2062, abs=1e-3) and upper[1] == pytest.approx(2.9224, abs=1e-3)


def test_poisson_ci_contains_coverage_and_handles_missing_values():
    coverage = [0.5, 0.5, 0.5, np.nan, 0.5]
    reads = [3, 0, np.nan, 100, -5]
    lower, upper = poisson_coverage_ci(coverage, reads)
    assert lower[0] <= 0.5 <= upper[0]
    assert np.isnan(lower[1:]).all() and np.isnan(upper[1:]).all()


def test_poisson_ci_is_deterministic_and_narrows_with_more_reads():
    coverage, reads = np.full(25, 1.0), np.linspace(10, 5000, 25)
    a = poisson_coverage_ci(coverage, reads)
    b = poisson_coverage_ci(coverage, reads)
    np.testing.assert_array_equal(a[0], b[0])
    width = a[1] - a[0]
    assert (np.diff(width) < 0).all()  # more reads -> narrower interval


def test_depth_table_without_profile_has_nan_errors(tmp_path):
    strains = make_strains(tmp_path)
    table = sylph_to_depth_table(strains, {"s1_x.fq.gz": "s1_x"})
    assert table["lower_ci"].isna().all() and table["upper_ci"].isna().all()
    assert list(table["sample"]) == ["s1_x", "s1_x"]


# ----------------------------------------------------------------------
# Sample matching
# ----------------------------------------------------------------------
@pytest.mark.parametrize("name", [
    "test_1_input_1",                       # SingleM: forward file name without extension
    "/data/reads/test_1_input_1.fastq.gz",  # Sylph: path of the reads
    "test_1_input_R1.fastq.gz",
    "test_1",
])
def test_sample_names_are_matched_to_fastq_samples(name):
    assert build_sample_map([name], ["test_1", "test_2"])[name] == "test_1"


def test_unmatched_sample_lists_names_and_mentions_sample_map():
    with pytest.raises(ValueError, match="sylph_sample"):
        build_sample_map(["unknown"], ["test_1"], source="Sylph")


def test_manual_map_overrides_and_accepts_base_name():
    mapping = build_sample_map(["/x/odd.fq.gz"], ["test_1"], manual_map={"odd.fq.gz": "test_1"}, source="Sylph")
    assert mapping == {"/x/odd.fq.gz": "test_1"}


# ----------------------------------------------------------------------
# Robust sample naming: the sample_id of the DNA mass file is the source of truth
# ----------------------------------------------------------------------
def test_paired_files_with_mate_numbers_belong_to_one_sample():
    # ENA-style names: the sample ID has no underscore, the mate number is separated by one
    files = ["/r/ERR1543552_1.fastq.gz", "/r/ERR1543552_2.fastq.gz", "/r/ERR1543553_1.fastq.gz"]
    sample_by_file, unmatched = assign_files_to_samples(files, ["ERR1543552", "ERR1543553"])
    assert unmatched == []
    assert group_by_sample(sample_by_file) == {
        "ERR1543552": ["/r/ERR1543552_1.fastq.gz", "/r/ERR1543552_2.fastq.gz"],
        "ERR1543553": ["/r/ERR1543553_1.fastq.gz"],
    }


def test_single_end_samples_ending_in_a_number_are_not_merged():
    # 'sample_1' and 'sample_2' are two Nanopore samples, not the mates of 'sample'
    files = ["sample_1.fastq.gz", "sample_2.fastq.gz"]
    sample_by_file, _ = assign_files_to_samples(files, ["sample_1", "sample_2"])
    assert sample_by_file == {"sample_1.fastq.gz": "sample_1", "sample_2.fastq.gz": "sample_2"}


@pytest.mark.parametrize("sample_id, name", [
    ("Patient_A_001", "Patient_A_001_R1.fastq.gz"),   # more than two underscore-separated parts
    ("S1.2", "S1.2_R1.fastq.gz"),                      # dot in the sample ID
    ("plate-1_B02", "plate-1_B02.sorted.bam"),         # hyphen in the sample ID
    ("A", "A.fastq.gz"),                               # single-character ID, no tags at all
    ("ERR1543552", "ERR1543552.sorted.bam"),
])
def test_sample_ids_can_have_any_structure(sample_id, name):
    assert match_sample_id(name, [sample_id, "other"]) == sample_id


def test_similar_sample_ids_are_not_confused():
    ids = ["S1", "S10", "S1_rep"]
    assert match_sample_id("S10_R1.fastq.gz", ids) == "S10"       # 'S1' is not followed by a separator
    assert match_sample_id("S1_rep_R1.fastq.gz", ids) == "S1_rep"  # the longest ID wins
    assert match_sample_id("S1_R1.fastq.gz", ids) == "S1"
    assert match_sample_id("S2_R1.fastq.gz", ids) is None


def test_unmatched_files_are_returned():
    sample_by_file, unmatched = assign_files_to_samples(["a_1.fq", "b_1.fq"], ["a"])
    assert sample_by_file == {"a_1.fq": "a"} and unmatched == ["b_1.fq"]


def test_manual_map_works_for_file_names_and_is_validated():
    assert match_sample_id("/x/odd.fq.gz", ["s1"], {"odd": "s1"}) == "s1"
    with pytest.raises(ValueError, match="not a sample_id"):
        match_sample_id("odd.fq.gz", ["s1"], {"odd": "typo"})


def test_duplicated_or_empty_sample_ids_are_rejected():
    assert check_sample_ids(["a", "b"]) == ["a", "b"]
    with pytest.raises(ValueError, match="Duplicated"):
        check_sample_ids(["a", "a"])
    with pytest.raises(ValueError, match="empty"):
        check_sample_ids(["a", float("nan")])


def test_strip_extension_keeps_dots_in_the_name():
    assert strip_extension("/d/S1.2_R1.fastq.gz") == "S1.2_R1"
    assert strip_extension("x.sorted.bam") == "x"
    assert strip_extension("x.sort.bam") == "x"
    assert strip_extension("x.sylphmpa") == "x"
    assert strip_extension("no_extension") == "no_extension"


def test_filtered_bam_path_handles_all_bam_extensions():
    assert filtered_bam_path("/in/s1.sorted.bam", "/out", 97.0) == "/out/s1_filtered97.bam"
    assert filtered_bam_path("/in/s1.sort.bam", "/out", 97.0) == "/out/s1_filtered97.bam"


def test_cli_fastq_mode_with_ena_style_names(tmp_path, monkeypatch):
    """End to end: ERR1543552_1/_2.fastq.gz (no underscore in the sample ID) in SingleM mode."""
    import gzip
    import sys

    from MGCalibrator.cli import main

    fastq = tmp_path / "fastq"
    singlem = tmp_path / "singlem"
    fastq.mkdir()
    singlem.mkdir()
    read = "@r\nACGTACGTAC\n+\nIIIIIIIIII\n"
    for name in ["ERR1543552_1", "ERR1543552_2", "Patient_A_001_R1", "Patient_A_002_R1"]:
        with gzip.open(fastq / f"{name}.fastq.gz", "wt") as handle:
            handle.write(read * 100)
    species = "Root; d__Bacteria; g__Escherichia; s__Escherichia coli"
    (singlem / "profile.tsv").write_text(
        "sample\tcoverage\ttaxonomy\n"
        f"ERR1543552_1\t2.0\t{species}\n"
        f"Patient_A_001_R1\t1.0\t{species}\n"
        f"Patient_A_002_R1\t3.0\t{species}\n"
    )
    mass = tmp_path / "dna_mass.csv"
    mass.write_text("sample_id,DNA_mass\nERR1543552,10\nPatient_A_001,20\nPatient_A_002,30\n")
    output = tmp_path / "out" / "result.csv"

    monkeypatch.setattr(sys, "argv", [
        "mgcalibrator", "--singlem_folder", str(singlem), "--fastq_folder", str(fastq),
        "--dna_mass", str(mass), "--scaling_factors", str(tmp_path / "out" / "sf.txt"),
        "--output", str(output),
    ])
    main()

    result = pd.read_csv(output).set_index("sample")
    assert set(result.index) == {"ERR1543552", "Patient_A_001", "Patient_A_002"}
    # Both mates are counted for ERR1543552 (2 x 1000 bases), one file for the others (1000 bases)
    scaling = dict(line.split(": ") for line in (tmp_path / "out" / "sf.txt").read_text().splitlines())
    assert float(scaling["ERR1543552"]) == pytest.approx(10 / (2000 * 1.079e-12))
    assert float(scaling["Patient_A_001"]) == pytest.approx(20 / (1000 * 1.079e-12))


def test_cli_sylph_mode_combines_poisson_error_with_qubit_error(tmp_path, monkeypatch):
    """End to end: Poisson CI from read counts, combined in quadrature with the Qubit error."""
    import sys

    from MGCalibrator.cli import main

    fastq = tmp_path / "fastq"
    sylph = tmp_path / "sylph"
    fastq.mkdir()
    sylph.mkdir()
    (fastq / "s1_R1.fastq").write_text("@r\nACGTACGTAC\n+\nIIIIIIIIII\n" * 100)
    write_sylphmpa(sylph / "s1.sylphmpa", "s1_R1.fastq", [("d__A|s__E coli|t__GCF_000005845.2", "0.5")])
    (sylph / "profile.tsv").write_text(
        PROFILE_HEADER.format(cov="True_cov") + profile_row("s1_R1.fastq", ECOLI, "0.5", "1000")
    )
    (tmp_path / "dna_mass.csv").write_text("sample_id,DNA_mass\ns1,10\n")
    output = tmp_path / "out.csv"
    monkeypatch.setattr(sys, "argv", [
        "mgcalibrator", "--sylph_folder", str(sylph), "--sylph_profile_folder", str(sylph / "profile.tsv"),
        "--fastq_folder", str(fastq), "--dna_mass", str(tmp_path / "dna_mass.csv"),
        "--scaling_factors", str(tmp_path / "sf.txt"), "--output", str(output), "--qubit_error", "0.1",
    ])
    main()

    row = pd.read_csv(output).iloc[0]
    assert row.read_counts == 1000 and row.error_note != row.error_note  # no note: NaN after reading the CSV
    poisson_down, poisson_up = row.error_down / row.depth, row.error_up / row.depth
    assert 0.03 < poisson_down < 0.09 and 0.03 < poisson_up < 0.09  # about 1.96 / sqrt(1000) = 6%
    assert row.calibrated_error_down / row.calibrated_depth == pytest.approx(np.hypot(poisson_down, 0.1))
    assert row.calibrated_error_up / row.calibrated_depth == pytest.approx(np.hypot(poisson_up, 0.1))


# ----------------------------------------------------------------------
# Reproducibility, output columns, --min_depth, scaling factor file, run info
# ----------------------------------------------------------------------
def test_make_rng_is_reproducible_and_independent_of_processing_order():
    from MGCalibrator.processor import make_rng

    assert make_rng(5, "s1").integers(1 << 30) == make_rng(5, "s1").integers(1 << 30)
    assert make_rng(5, "s1").integers(1 << 30) != make_rng(5, "s2").integers(1 << 30)
    assert make_rng(5, "s1").integers(1 << 30) != make_rng(6, "s1").integers(1 << 30)


def test_monte_carlo_depth_is_reproducible_with_seed():
    from MGCalibrator.processor import MC_simulation_for_depth, make_rng

    depths = np.full(1000, 5, dtype=np.int32)
    reads = np.array([[i * 10, 100] for i in range(50)])
    first = MC_simulation_for_depth(depths, reads, n_simulations=20, rng=make_rng(3, "s1"))
    second = MC_simulation_for_depth(depths, reads, n_simulations=20, rng=make_rng(3, "s1"))
    assert first == second


def test_apply_min_depth():
    from MGCalibrator.cli import apply_min_depth

    depths = pd.DataFrame({"depth": [0.01, 0.5, 2.0]})
    assert apply_min_depth(depths, 0)[1] == 0
    kept, n_removed = apply_min_depth(depths, 0.1)
    assert list(kept["depth"]) == [0.5, 2.0] and n_removed == 1
    with pytest.raises(ValueError, match="No rows left"):
        apply_min_depth(depths, 10)


def test_scaling_factors_file_is_reused_with_warning_and_can_be_overwritten(tmp_path, caplog):
    from MGCalibrator.processor import calculate_scaling_factors

    fastq = tmp_path / "s1.fastq"
    fastq.write_text("@r\nACGTACGTAC\n+\nIIIIIIIIII\n" * 100)  # 1000 bases
    factors_file = tmp_path / "sf.txt"
    factors_file.write_text("s1: 123.0\nold_sample: 5.0\n")
    kwargs = dict(samples={"s1"}, sample_by_bam=None, initial_dna_mass={"s1": 10.0},
                  scaling_factors_file=str(factors_file), fastq_files_by_sample={"s1": [str(fastq)]})

    factors, reused = calculate_scaling_factors(**kwargs)
    assert factors["s1"] == 123.0 and reused == ["s1"]
    assert "--overwrite_scaling_factors" in caplog.text and "old_sample" in caplog.text

    factors, reused = calculate_scaling_factors(**kwargs, overwrite=True)
    assert reused == [] and factors["s1"] == pytest.approx(10.0 / (1000 * 1.079e-12))
    assert "old_sample" not in factors_file.read_text()


def run_cli(monkeypatch, *arguments):
    import sys

    from MGCalibrator.cli import main

    monkeypatch.setattr(sys, "argv", ["mgcalibrator", *map(str, arguments)])
    main()


def make_singlem_dataset(tmp_path):
    import gzip

    fastq, singlem = tmp_path / "fastq", tmp_path / "singlem"
    fastq.mkdir()
    singlem.mkdir()
    with gzip.open(fastq / "S1_1.fastq.gz", "wt") as handle:
        handle.write("@r\nACGTACGTAC\n+\nIIIIIIIIII\n" * 100)
    (singlem / "profile.tsv").write_text(
        "sample\tcoverage\ttaxonomy\n"
        "S1_1\t2.0\tRoot; d__Bacteria; g__Escherichia; s__Escherichia coli\n"
        "S1_1\t0.01\tRoot; d__Bacteria; g__Bacillus; s__Bacillus subtilis\n"
    )
    (tmp_path / "dna_mass.csv").write_text("sample_id,DNA_mass\nS1,10\n")
    return fastq, singlem


def test_cli_output_columns_error_note_min_depth_and_run_info(tmp_path, monkeypatch):
    import json

    fastq, singlem = make_singlem_dataset(tmp_path)
    output = tmp_path / "out" / "result.csv"
    common = ["--singlem_folder", singlem, "--fastq_folder", fastq, "--dna_mass", tmp_path / "dna_mass.csv",
              "--scaling_factors", tmp_path / "out" / "sf.txt", "--output", output]
    run_cli(monkeypatch, *common, "--min_depth", "0.1")

    result = pd.read_csv(output)
    assert list(result.columns) == [
        "sample", "reference", "depth", "error_down", "error_up", "calibrated_depth",
        "calibrated_error_down", "calibrated_error_up", "scaling_factor", "read_counts", "error_note",
    ]
    assert len(result) == 1  # the row with depth 0.01 is removed
    assert result.loc[0, "scaling_factor"] == pytest.approx(10 / (1000 * 1.079e-12))
    assert result.loc[0, "error_note"] == "no_error_estimate_for_singlem"

    info = json.loads((tmp_path / "out" / "result_run_info.json").read_text())
    assert info["input_mode"] == "SingleM + FASTQ"
    assert info["parameters"]["min_depth"] == 0.1 and info["seed"] is None
    assert info["data"]["n_rows_removed_by_min_depth"] == 1
    assert info["scaling_factors"]["calculated_in_this_run"] == ["S1"]
    assert {"numpy", "pandas", "scipy", "pysam"} <= set(info["packages"])

    # Second run: the scaling factor is now taken from the file
    run_cli(monkeypatch, *common)
    info = json.loads((tmp_path / "out" / "result_run_info.json").read_text())
    assert info["scaling_factors"]["reused_from_file"] == ["S1"]


def test_cli_sylph_error_note_without_profile(tmp_path, monkeypatch):
    fastq = tmp_path / "fastq"
    sylph = tmp_path / "sylph"
    fastq.mkdir()
    sylph.mkdir()
    (fastq / "s1_R1.fastq").write_text("@r\nACGTACGTAC\n+\nIIIIIIIIII\n" * 100)
    write_sylphmpa(sylph / "s1.sylphmpa", "s1_R1.fastq", [("d__A|s__E coli|t__GCF_000005845.2", "0.5")])
    (tmp_path / "dna_mass.csv").write_text("sample_id,DNA_mass\ns1,10\n")
    output = tmp_path / "out.csv"
    run_cli(monkeypatch, "--sylph_folder", sylph, "--fastq_folder", fastq, "--dna_mass", tmp_path / "dna_mass.csv",
            "--scaling_factors", tmp_path / "sf.txt", "--output", output)
    result = pd.read_csv(output)
    assert result.loc[0, "error_note"] == "no_sylph_profile"
    assert pd.isna(result.loc[0, "read_counts"]) or result.loc[0, "read_counts"] == "NA"


def test_tool_version_handles_non_utf8_output_and_missing_tools(tmp_path):
    import os
    import stat

    from MGCalibrator.runinfo import _tool_version

    tool = tmp_path / "faketool"
    tool.write_bytes(b"#!/bin/sh\nprintf 'faketool 1.2 \\xa9 copyright\\\\n'\n")
    tool.chmod(tool.stat().st_mode | stat.S_IEXEC)
    assert _tool_version(str(tool)).startswith("faketool 1.2")
    assert _tool_version(str(tmp_path / "does_not_exist")) is None
