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
    add_confidence_intervals,
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


def profile_row(sample, genome, cov, lam, ci, contig="NC_1 Escherichia coli"):
    return (f"{sample}\t{genome}\t50\t50\t99.0\t{cov}\tNA-NA\t{lam}\t{ci}\t1\t1.0\t10/20\t99.0\t0\t{contig}\n")


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
# Sylph: profile TSV -> confidence intervals
# ----------------------------------------------------------------------
def make_strains(tmp_path, sample="s1_x.fq.gz"):
    f = tmp_path / "s.sylphmpa"
    write_sylphmpa(f, sample, [
        ("d__A|s__E coli|t__GCF_000005845.2", "0.412"),
        ("d__A|s__B subtilis|t__GCF_000009045.1", "12.35"),
    ])
    return read_sylphmpa_files([str(f)])


def test_confidence_intervals_joined_on_sample_and_accession(tmp_path):
    p = tmp_path / "profile.tsv"
    p.write_text(
        PROFILE_HEADER.format(cov="Eff_cov")
        + profile_row("s1_x.fq.gz", "GCF_000005845.2_ASM584v2_genomic.fna.gz", "0.412", "0.412", "0.25-0.61")
        + profile_row("s1_x.fq.gz", "GCF_000009045.1_ASM904v1_genomic.fna.gz", "12.350", "HIGH", "NA-NA")
        + profile_row("other.fq.gz", "GCF_000005845.2_ASM584v2_genomic.fna.gz", "0.412", "0.412", "0.01-0.99")
    )
    result = add_confidence_intervals(make_strains(tmp_path), read_sylph_profiles([str(p)]))
    ecoli, bsub = result.iloc[0], result.iloc[1]
    assert (ecoli["lower_ci"], ecoli["upper_ci"]) == (0.25, 0.61)
    assert np.isnan(bsub["lower_ci"]) and np.isnan(bsub["upper_ci"])  # HIGH coverage: no interval


def test_true_cov_rescales_interval(tmp_path):
    p = tmp_path / "profile.tsv"
    # True_cov (0.824) is twice Eff_lambda (0.412), so the interval is doubled
    p.write_text(
        PROFILE_HEADER.format(cov="True_cov")
        + profile_row("s1_x.fq.gz", "GCF_000005845.2_ASM584v2_genomic.fna.gz", "0.412", "0.206", "0.10-0.30")
    )
    profiles = read_sylph_profiles([str(p)])
    assert profiles.loc[0, "ci_low"] == pytest.approx(0.20)
    assert profiles.loc[0, "ci_high"] == pytest.approx(0.60)


def test_inconsistent_coverage_gives_no_interval(tmp_path):
    p = tmp_path / "profile.tsv"
    p.write_text(
        PROFILE_HEADER.format(cov="Eff_cov")
        + profile_row("s1_x.fq.gz", "GCF_000005845.2_ASM584v2_genomic.fna.gz", "0.900", "0.900", "0.5-1.2")
    )
    result = add_confidence_intervals(make_strains(tmp_path), read_sylph_profiles([str(p)]))
    assert result["lower_ci"].isna().all()


def test_non_profile_tsv_is_skipped(tmp_path):
    good = tmp_path / "profile.tsv"
    good.write_text(PROFILE_HEADER.format(cov="Eff_cov")
                    + profile_row("s", "GCF_1.1_ASM1v1_genomic.fna.gz", "1.000", "1.000", "0.5-1.5"))
    other = tmp_path / "merged.tsv"
    other.write_text("clade_name\ts1\nd__A\t1\n")
    assert len(read_sylph_profiles([str(good), str(other)])) == 1
    with pytest.raises(ValueError):
        read_sylph_profiles([str(other)])


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
