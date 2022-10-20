import os
import stat
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
from .utils_for_testing import mock_log, write_indexed_bcf


TEST_DIR = os.path.join(
    "strainflye",
    "tests",
    "inputs",
    "small",
)
FASTA = os.path.join(TEST_DIR, "contigs.fasta")
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


def test_load_triplet_good(capsys):
    n2l, bam_obj, bcf_obj, tt, tm = mu.load_triplet(FASTA, BAM, BCF, mock_log)
    # The actual fasta / bam / bcf loading done in this file isn't super
    # important to test, since these are either already tested elsewhere (the
    # fasta and bcf loading) or are completely done by an external library (bam
    # loading).
    #
    # So we just cursorily check that these loaded files look ok. the main
    # thing we wanna check is that the logged output is correct.
    assert n2l == {"c1": 23, "c2": 12, "c3": 16}
    exp_contigs = set(["c1", "c2", "c3"])
    assert set(bam_obj.references) == exp_contigs
    assert set(bcf_obj.header.contigs) == exp_contigs

    assert tt is None
    assert tm is None

    captured = capsys.readouterr()
    assert captured.out == (
        "PREFIX\nMockLog: Loading and checking FASTA, BAM, and BCF files...\n"
        "MockLog: The FASTA file describes 3 contig(s).\n"
        "MockLog: All FASTA contig(s) are included in the BAM file (this BAM "
        "file has 3 reference(s)).\n"
        "MockLog: All FASTA contig(s) are included in the BCF file (the "
        "header of this BCF file describes 3 contig(s)).\n"
        "MockLog: The lengths of all contig(s) in the FASTA file match the "
        "corresponding lengths in the BAM and BCF files.\n"
        "MockLog: So far, these files seem good.\n"
    )


def test_load_triplet_fasta_not_in_bcf_all():
    # This is actually a valid VCF file (we use it in the parse_sf_bcf()
    # tests); the problem with it is that it doesn't describe c1, c2, or c3.
    with write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=1000>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    ) as fh:
        with pytest.raises(ParameterError) as ei:
            mu.load_triplet(FASTA, BAM, fh.name, mock_log)
        assert str(ei.value) == (
            "All contigs in the FASTA file must also be contained in the BCF "
            "file."
        )


def test_load_triplet_fasta_not_in_bcf_partial():
    with write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=c2,length=12>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "c2\t1\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "c2\t4\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "c2\t5\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "c2\t9\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    ) as fh:
        with pytest.raises(ParameterError) as ei:
            mu.load_triplet(FASTA, BAM, fh.name, mock_log)
        assert str(ei.value) == (
            "All contigs in the FASTA file must also be contained in the BCF "
            "file."
        )


def test_load_triplet_fasta_not_in_bam_all():
    # The BAM check should get performed before the BCF check, so we can test
    # the case where FASTA contigs aren't in the BAM by just using a different
    # FASTA file
    with pytest.raises(ParameterError) as ei:
        mu.load_triplet(StringIO(">missing_contig\nACGT"), BAM, BCF, mock_log)
    assert str(ei.value) == (
        "All contigs in the FASTA file must also be contained in the BAM "
        "file."
    )


def test_check_executable_nonexecutable_file(tmp_path):
    nonex_file = tmp_path / "nonex_file.sh"
    with open(nonex_file, "w") as f:
        f.write('#! /usr/bin/env bash\necho "hi"\n')
    with pytest.raises(ParameterError) as ei:
        mu.check_executable(nonex_file)
    assert str(ei.value) == f"{nonex_file} is not executable."


def test_check_executable_actually_executable_file(tmp_path):
    exf = tmp_path / "ex"
    with open(exf, "w") as f:
        f.write('#! /usr/bin/env bash\necho "hi"\n')
    # yoinked over from test_smooth_utils
    new_mode_bits = (
        os.stat(exf).st_mode | stat.S_IXGRP | stat.S_IXUSR | stat.S_IXOTH
    )
    os.chmod(exf, new_mode_bits)
    mu.check_executable(exf)
