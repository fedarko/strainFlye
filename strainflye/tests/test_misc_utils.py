import os
import shutil
import tempfile
import pytest
import pysam
import pandas as pd
import strainflye.misc_utils as mu
from io import StringIO
from strainflye.bcf_utils import parse_sf_bcf
from strainflye.errors import ParameterError
from strainflye.config import DI_PREF


TEST_DIR = os.path.join(
    "strainflye",
    "tests",
    "inputs",
    "small",
)
BAM = os.path.join(TEST_DIR, "alignment.bam")
BCF = os.path.join(TEST_DIR, "call-r-min3-di12345", "naive-calls.bcf")


def test_make_output_dir_exists():
    with tempfile.TemporaryDirectory() as tdname:
        assert os.path.exists(tdname)
        child_file_name = os.path.join(tdname, "test_file.txt")
        with open(child_file_name, "w") as f:
            f.write("please don't overwrite me!")

        mu.make_output_dir(tdname)

        # Check that this didn't overwrite anything (the directory itself, or
        # any of its contents)
        assert os.path.exists(tdname)
        assert os.path.exists(child_file_name)
        with open(child_file_name, "r") as f:
            ctext = f.read()
        assert ctext == "please don't overwrite me!"


def test_make_output_dir_notexists_multilevel():
    tdname_top = os.path.join(tempfile.gettempdir(), "_strainflye_td1")
    tdname_bot = os.path.join(tdname_top, "_strainflye_td2")
    try:
        assert not os.path.exists(tdname_top)
        assert not os.path.exists(tdname_bot)
        mu.make_output_dir(tdname_bot)
        assert os.path.exists(tdname_top)
        assert os.path.exists(tdname_bot)
    finally:
        # https://stackoverflow.com/a/13118112
        shutil.rmtree(tdname_top)


def test_verify_contig_subset_raises_error():
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_subset(set("abcdef"), set("abdef"), "s1", "s2")
    assert str(ei.value) == "All contigs in s1 must also be contained in s2."


def test_verify_contig_subset_good():
    # Sets are identical
    mu.verify_contig_subset(set("abcdef"), set("abcdef"), "s1", "s2")
    mu.verify_contig_subset(set("a"), set("a"), "s1", "s2")

    # Proper subset
    mu.verify_contig_subset(set(""), set(""), "s1", "s2")
    mu.verify_contig_subset(set(""), set("ab"), "s1", "s2")
    mu.verify_contig_subset(set("a"), set("ab"), "s1", "s2")
    mu.verify_contig_subset(set("ab"), set("abc"), "s1", "s2")


def test_verify_contig_subset_exact():
    # Check that exact doesn't change mandate that child is subset of parent
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_subset(
            set("abcdef"), set("abdef"), "s1", "s2", exact=True
        )
    assert str(ei.value) == "All contigs in s1 must also be contained in s2."

    # Now, check that exact ensures that the two sets are equal, even if child
    # is a subset of parent
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_subset(
            set("abcd"), set("abcdef"), "s1", "s2", exact=True
        )
    assert str(ei.value) == "All contigs in s2 must also be contained in s1."

    # Now, check that exact succeeds
    mu.verify_contig_subset(set("abcd"), set("abcd"), "s1", "s2", exact=True)
    mu.verify_contig_subset(set(""), set(""), "s1", "s2", exact=True)


def test_verify_contig_lengths_bam_bcf_missing():
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_lengths({"c1": 12, "c2": 5})
    assert str(ei.value) == "Neither bam_obj nor bcf_obj is provided."


def test_verify_contig_lengths_good_both():
    bam_obj = pysam.AlignmentFile(BAM, "rb")
    bcf_obj, tt, tm = parse_sf_bcf(BCF)
    mu.verify_contig_lengths(
        {"c1": 23, "c2": 12, "c3": 16}, bam_obj=bam_obj, bcf_obj=bcf_obj
    )

    # check that subsets are ok
    mu.verify_contig_lengths({"c2": 12}, bam_obj=bam_obj, bcf_obj=bcf_obj)


def test_verify_contig_lengths_good_just_bam():
    bam_obj = pysam.AlignmentFile(BAM, "rb")
    mu.verify_contig_lengths({"c1": 23, "c2": 12, "c3": 16}, bam_obj=bam_obj)


def test_verify_contig_lengths_good_just_bcf():
    bcf_obj, tt, tm = parse_sf_bcf(BCF)
    mu.verify_contig_lengths({"c1": 23, "c2": 12, "c3": 16}, bcf_obj=bcf_obj)


def test_verify_contig_lengths_mismatch_with_both():
    # We make the bam check before the bcf check, so the bam error should get
    # thrown. but this is less a hard mandate on how this function should
    # perform and more a quirk of our implementation, yada yada
    bam_obj = pysam.AlignmentFile(BAM, "rb")
    bcf_obj, tt, tm = parse_sf_bcf(BCF)
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_lengths(
            {"c1": 22, "c2": 12, "c3": 16}, bam_obj=bam_obj, bcf_obj=bcf_obj
        )

    assert str(ei.value) == (
        "Contig c1 has length 22 in the FASTA file, but length 23 in the BAM "
        "file."
    )


def test_verify_contig_lengths_mismatch_with_just_bam():
    bam_obj = pysam.AlignmentFile(BAM, "rb")
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_lengths(
            {"c1": 22, "c2": 12, "c3": 16}, bam_obj=bam_obj
        )

    assert str(ei.value) == (
        "Contig c1 has length 22 in the FASTA file, but length 23 in the BAM "
        "file."
    )


def test_verify_contig_lengths_mismatch_with_just_bcf():
    bcf_obj, tt, tm = parse_sf_bcf(BCF)
    with pytest.raises(ParameterError) as ei:
        mu.verify_contig_lengths(
            {"c1": 23, "c2": 12, "c3": 1024}, bcf_obj=bcf_obj
        )

    assert str(ei.value) == (
        "Contig c3 has length 1,024 in the FASTA file, but length 16 in the "
        "BCF file."
    )


def test_load_and_sanity_check_diversity_indices_not_enough_di_cols():
    # this function is already indirectly tested by the FDR utils tests,
    # because it originated as code in fdr_utils.autoselect_decoy().
    # So we can just test some of the weird cases that aren't covered by that.

    tsv_text = (
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35.2\t100\t0.5\t0.2\n"
    )

    # Case 1: Just the DI cols error
    tsv = StringIO(tsv_text)
    with pytest.raises(ParameterError) as ei:
        mu.load_and_sanity_check_diversity_indices(tsv, min_num_di_columns=3)
    assert str(ei.value) == (
        "Diversity indices file describes < 3 column(s) of diversity indices."
    )

    # Case 2: Both the DI cols and contigs error -- in this case, the contigs
    # error is found first
    # (we gotta recreate the StringIO because at this point it's "empty." [IDK
    # if "empty" is the correct term to use in this case but whatever you get
    # the idea])
    tsv = StringIO(tsv_text)
    with pytest.raises(ParameterError) as ei:
        mu.load_and_sanity_check_diversity_indices(
            tsv, min_num_di_columns=3, min_num_contigs=2
        )
    assert str(ei.value) == ("Diversity indices file describes < 2 contig(s).")


def test_load_and_sanity_check_diversity_indices_one_of_each():
    tsv = StringIO(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\n"
        "edge_1\t35.2\t100\t0.5\n"
    )
    obs_df = mu.load_and_sanity_check_diversity_indices(tsv)
    exp_df = pd.DataFrame(
        {
            "AverageCoverage": [35.2],
            "Length": [100],
            f"{DI_PREF}1": [0.5],
        },
        index=pd.Index(["edge_1"], name="Contig"),
    )
    pd.testing.assert_frame_equal(obs_df, exp_df)
