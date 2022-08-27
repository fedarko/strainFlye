import os
import tempfile
import pytest
import numpy as np
import pandas as pd
from strainflye.errors import ParameterError, SequencingDataError
from strainflye.call_utils import (
    get_alt_pos_info,
    get_pos_info_str,
    call_r_mutation,
    call_p_mutation,
    parse_di_list,
    get_min_sufficient_coverages_p,
    get_min_sufficient_coverages_r,
    run,
)
from strainflye.bcf_utils import parse_sf_bcf
from .utils_for_testing import mock_log


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BAM = os.path.join(IN_DIR, "alignment.bam")


def test_get_alt_pos_info():
    # Test data reused from
    # https://github.com/fedarko/sheepgut/blob/main/notebooks/test_pileup.py
    assert get_alt_pos_info({"A": 10, "C": 0, "G": 3, "T": 50}) == (
        63,
        10,
        "A",
        50,
        "T",
    )

    # Ties are broken arbitrarily, so there are multiple correct answers here
    tie_api = get_alt_pos_info({"A": 10, "C": 10, "G": 10, "T": 10})
    assert tie_api[0] == 40
    assert tie_api[1] == tie_api[3] == 10
    assert tie_api[2] in "ACGT"
    assert tie_api[4] in "ACGT"
    assert tie_api[2] != tie_api[4]

    tie_api_2 = get_alt_pos_info({"A": 5, "C": 10, "G": 10, "T": 0})
    assert tie_api_2[0] == 25
    assert tie_api_2[1] == tie_api_2[3] == 10
    assert tie_api_2[2] in "CG"
    assert tie_api_2[4] in "CG"
    assert tie_api_2[2] != tie_api_2[4]

    # Ties are broken arbitrarily, so there are multiple correct answers here
    tie_api = get_alt_pos_info({"A": 0, "C": 0, "G": 0, "T": 0})
    assert tie_api[0] == 0
    assert tie_api[1] == tie_api[3] == 0
    assert tie_api[2] in "ACGT"
    assert tie_api[4] in "ACGT"
    assert tie_api[2] != tie_api[4]


def test_call_r_mutation():
    # paranoia
    for r in range(1, 6):
        assert call_r_mutation(5, r)
        assert call_r_mutation(100, r)
    for r in range(6, 20):
        assert not call_r_mutation(5, r)
        assert call_r_mutation(100, r)
    assert call_r_mutation(1, 1)


def test_call_p_mutation():
    # bear in mind that these values of p are at this point scaled up by
    # 10,000, relative to (0, 0.5]. i.e. 5,000 --> 5,000 / 10,000 = 0.5,
    # aka p = 50%.
    assert call_p_mutation(5, 10, 5000, 2)
    assert not call_p_mutation(2, 10, 5000, 2)

    # 2 / 10 = 20%. This is a p-mutation for p <= 2,000 (p <= 20%).
    for p in range(1, 2001):
        assert call_p_mutation(2, 10, p, 2)
    for p in range(2001, 3001):
        assert not call_p_mutation(2, 10, p, 2)

    # 5 / 1,000 = 0.5%. This is a p-mutation for p <= 50 (p <= 0.5%).
    for p in range(1, 51):
        assert call_p_mutation(5, 1000, p, 2)
    for p in (51, 200):
        assert not call_p_mutation(5, 1000, p, 2)

    # 1 / 10,000 = 0.01%. This is a p-mutation for p == 1 (p == 0.01%).
    # We could in theory go lower, but... probs not needed
    assert call_p_mutation(1, 10000, 1, 1)
    assert not call_p_mutation(1, 10000, 2, 1)

    assert call_p_mutation(2, 10, 1, 2)
    # failure due to min alt pos
    assert not call_p_mutation(2, 10, 1, 100)


def test_get_pos_info_str():
    assert get_pos_info_str(5, 1000) == "MDP=1000;AAD=5"
    assert get_pos_info_str(1000, 10000) == "MDP=10000;AAD=1000"
    assert get_pos_info_str(1, 1) == "MDP=1;AAD=1"


def test_parse_di_list_whitespace_ok():
    assert parse_di_list("1\n, 2,  \t  3,\n5", "p") == [1, 2, 3, 5]


def test_parse_di_list_normal_ok():
    # These are the defaults, at least of writing
    assert parse_di_list("50,100,200,500,1000,2500,5000", "p") == [
        50,
        100,
        200,
        500,
        1000,
        2500,
        5000,
    ]
    assert parse_di_list("5,10,20,50,100,250,500", "r") == [
        5,
        10,
        20,
        50,
        100,
        250,
        500,
    ]


def test_parse_di_list_non_int():
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("1,2,burger,3", "p")
    assert (
        "We couldn't parse \"burger\". Doesn't seem to be an integer?"
    ) == str(errorinfo.value)

    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("1,2,1.5,3", "p")
    assert (
        "We couldn't parse \"1.5\". Doesn't seem to be an integer?"
    ) == str(errorinfo.value)

    # check that surrounding whitespace is ignored in the error message
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("1,2, 5\n9 ,3", "p")
    assert (
        "We couldn't parse \"5\n9\". Doesn't seem to be an integer?"
    ) == str(errorinfo.value)

    # check that empty entries are flagged -- this also
    # implicitly handles empty lists
    for bad_str in ("1, 2,    , 3", "", "1, ", "1,"):
        for param in ("p", "r"):
            with pytest.raises(ParameterError) as errorinfo:
                parse_di_list(bad_str, param)
            assert (
                "We couldn't parse \"\". Doesn't seem to be an integer?"
            ) == str(errorinfo.value)

    # when the r-value has the drip :O
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("3, 4, \u1f4af, 5", "r")
    assert (
        "We couldn't parse \"\u1f4af\". Doesn't seem to be an integer?"
    ) == str(errorinfo.value)


def test_parse_di_list_p_out_of_range():
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("0,1,2,3", "p")
    assert (
        "0 is not in the range (0, 5000], and is thus not a valid value "
        "of p."
    ) == str(errorinfo.value)

    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("100, 200, 300, 400, -20", "p")
    assert (
        "-20 is not in the range (0, 5000], and is thus not a valid value "
        "of p."
    ) == str(errorinfo.value)

    # The first problematic value we see should trigger the error
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("100, 200, 300, 400, 666666, -20", "p")
    assert (
        "666666 is not in the range (0, 5000], and is thus not a valid value "
        "of p."
    ) == str(errorinfo.value)


def test_parse_di_list_r_out_of_range():
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("0,1,2,3", "r")
    assert "0 is not >= 1, and is thus not a valid value of r." == str(
        errorinfo.value
    )

    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("100, 200, 300, 400, -20", "r")
    assert ("-20 is not >= 1, and is thus not a valid value of r.") == str(
        errorinfo.value
    )


def test_parse_di_list_nonunique():
    with pytest.raises(ParameterError) as errorinfo:
        parse_di_list("10,1,2,3, 1", "r")
    assert "The list of diversity index threshold values isn't unique." == str(
        errorinfo.value
    )


def test_parse_di_list_badparam():
    for bad_param in ("R", "P", "asdf", "hamborgar", 1, "2", None, ""):
        with pytest.raises(ParameterError) as errorinfo:
            parse_di_list("10,1,2,3", bad_param)
        assert 'param must be either "p" or "r".' == str(errorinfo.value)


def test_parse_di_list_few_entries():
    assert parse_di_list("  1  ", "r") == [1]


def test_parse_di_list_sorts_results():
    assert parse_di_list("  10, 3  ", "r") == [3, 10]
    assert parse_di_list("1, 20, 5, 10, 100, 3, 15, 1000, 39", "r") == [
        1,
        3,
        5,
        10,
        15,
        20,
        39,
        100,
        1000,
    ]


def test_get_min_sufficient_coverages_p():
    assert get_min_sufficient_coverages_p(
        [100, 200, 300, 2000, 5000, 1, 50], 5
    ) == [500, 250, (500 / 3), 25, 10, 50000, 1000]


def test_get_min_sufficient_coverages_r():
    r_vals = [100, 200, 300, 2000, 5000, 1, 50]
    assert get_min_sufficient_coverages_r(r_vals, 2) == [r * 2 for r in r_vals]
    assert get_min_sufficient_coverages_r(r_vals, 1) == r_vals


def test_run_p_r_conflict():
    with pytest.raises(ParameterError) as errorinfo:
        run(
            "c",
            "b",
            "od",
            mock_log,
            True,
            min_p=5,
            min_r=6,
            min_alt_pos=2,
            div_index_p_list=[1, 2, 3],
            div_index_r_list=[1, 2, 3],
            min_read_number=5,
        )
    assert (
        "p and r can't be specified at the same time. Please choose one."
        == str(errorinfo.value)
    )

    with pytest.raises(ParameterError) as errorinfo:
        run("c", "b", "od", mock_log, True)
    assert "Either p or r needs to be specified." == str(errorinfo.value)


def test_run_no_muts_called():
    with tempfile.TemporaryDirectory() as td:
        with pytest.raises(ParameterError) as ei:
            run(
                FASTA,
                BAM,
                td,
                mock_log,
                True,
                min_r=200,
                div_index_r_list=[1, 2, 3, 4, 5, 6, 7],
                min_cov_factor=2,
            )
        assert str(ei.value) == (
            "Could not call any mutations across any of the 3 contigs. Try "
            "adjusting the minimum threshold value, or other parameters?"
        )


def test_run_small_dataset_r(capsys):
    # "Yeah integration tests are good I guess" -- Sun Tzu
    with tempfile.TemporaryDirectory() as td:
        run(
            FASTA,
            BAM,
            td,
            mock_log,
            True,
            min_r=2,
            div_index_r_list=[1, 2, 3, 4, 5, 6, 7],
            min_cov_factor=2,
        )

        # Check BCF
        bcf_obj, tt, tm = parse_sf_bcf(os.path.join(td, "naive-calls.bcf"))
        assert tt == "r"
        assert tm == 2
        # There should be exactly five mutations, three on c1 and two on c2.
        # See the logging checks later on in this test for details. (Note that
        # the positions here are 1-indexed, like the ones in the comments
        # below.)
        #
        # These lists of expected values are in order for
        # each of these mutations: so the first expected mutation is in c1, at
        # position 3, with ref (consensus) G and alt T, coverage (MDP) 12x, and
        # alt(pos) of 3
        exp_contigs = ["c1", "c1", "c1", "c3", "c3"]
        exp_pos = [4, 11, 13, 7, 8]
        exp_ref = ["G", "G", "G", "A", "T"]
        exp_mdp = [12, 12, 12, 13, 13]
        exp_alt = ["T", "A", "A", "T", "C"]
        exp_aad = [3, 5, 5, 6, 3]
        for ri, rec in enumerate(bcf_obj.fetch()):
            assert rec.contig == exp_contigs[ri]
            assert rec.pos == exp_pos[ri]
            assert rec.ref == exp_ref[ri]
            assert rec.info.get("MDP") == exp_mdp[ri]

            # Contingent upon the whole only-1-mut-at-a-position thing
            assert len(rec.alts) == 1
            assert rec.alts[0] == exp_alt[ri]
            rec_aad = rec.info.get("AAD")
            assert len(rec_aad) == 1
            assert rec_aad[0] == exp_aad[ri]

        # Check diversity index file
        obs_di_df = pd.read_csv(
            os.path.join(td, "diversity-indices.tsv"), sep="\t", index_col=0
        )

        # All of c1's positions are covered at 12x.
        c1_avg_cov = 12.0

        # Out of c2's 12 positions, 10 are covered at 11x, and the first two
        # are covered at just 10x (one of the reads doesn't cover them).
        c2_avg_cov = ((11 * 10) + (10 * 2)) / 12

        # Out of c3's 16 positions, 15 are covered at 13x; the last position is
        # covered at just 2x.
        c3_avg_cov = ((13 * 15) + 2) / 16

        exp_di_df = pd.DataFrame(
            {
                "AverageCoverage": [c1_avg_cov, c2_avg_cov, c3_avg_cov],
                "Length": [23, 12, 16],
                "DivIdx(r=1,minSuffCov=2)": [3 / 23, 0, 3 / 16],
                "DivIdx(r=2,minSuffCov=4)": [3 / 23, 0, 2 / 15],
                "DivIdx(r=3,minSuffCov=6)": [3 / 23, 0, 2 / 15],
                "DivIdx(r=4,minSuffCov=8)": [2 / 23, 0, 1 / 15],
                "DivIdx(r=5,minSuffCov=10)": [2 / 23, 0, 1 / 15],
                "DivIdx(r=6,minSuffCov=12)": [0, np.nan, 1 / 15],
                "DivIdx(r=7,minSuffCov=14)": [np.nan, np.nan, np.nan],
            },
            index=pd.Index(["c1", "c2", "c3"], name="Contig"),
        )
        pd.testing.assert_frame_equal(obs_di_df, exp_di_df)

    # Verify that the log messages written out match up with our expectations.
    # This is, of course, very brittle -- any changes to these messages will
    # break this test. But I don't imagine these messages will need to change
    # too much (knock on wood).

    # Basic intro stuff
    exp_out = (
        "PREFIX\nMockLog: Loading and checking contig information...\n"
        "MockLog: The FASTA file describes 3 contigs.\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "references), with the same lengths.\n"
        "PREFIX\nMockLog: Creating and writing diversity index and VCF file "
        "headers...\n"
        "MockLog: Wrote out diversity index and VCF file headers.\n"
        "MockLog: (We'll convert the VCF file to an indexed BCF afterwards.)\n"
    )
    # Messages from during calling.
    # For reference, here's some info about the mutated positions in these
    # contigs (using 1-indexing).
    #
    # -c1 has three mutated positions:
    #  -Position  4: 1 A, 2 C, 3 T, 6 G
    #  -Position 11: 5 A, 7 G
    #  -Position 13: 5 A, 7 G
    #
    # -c2 has no mutated positions. It does, however, have an "unreasonable"
    #  position at its position 1: this position is an A in the reference c2
    #  sequence, but it is completely a T in the alignment. strainFlye should
    #  categorize this as not a mutation.
    #
    # -c3 has two mutated positions:
    #  -Position 7: 6 T,  7 A
    #  -Position 8: 3 C, 10 T
    #  -Additionally, the final position in c3 (16) is skipped with a deletion
    #   in most of the alignments to c3; only two alignments cover this
    #   position (1 G, 1 T). This doesn't pass the min r threshold that we set.
    exp_out += (
        "PREFIX\nMockLog: Running r-mutation calling (--min-r = 2) and "
        "computing diversity indices...\n"
        "MockLog: On contig c1 (23 bp) (1 / 3 = 33.33% done).\n"
        "MockLog: Called 3 r-mutation(s) (using --min-r = 2) in contig c1.\n"
        "MockLog: 6 / 7 diversity indices were defined for contig c1.\n"
        "MockLog: On contig c2 (12 bp) (2 / 3 = 66.67% done).\n"
        "MockLog: Called 0 r-mutation(s) (using --min-r = 2) in contig c2.\n"
        "MockLog: 5 / 7 diversity indices were defined for contig c2.\n"
        "MockLog: On contig c3 (16 bp) (3 / 3 = 100.00% done).\n"
        "MockLog: Called 2 r-mutation(s) (using --min-r = 2) in contig c3.\n"
        "MockLog: 6 / 7 diversity indices were defined for contig c3.\n"
    )

    # Basic wrap-up stuff (making sure to include the logging statements from
    # bcf_utils.index_bcf())
    exp_out += (
        "MockLog: Done running r-mutation calling (--min-r = 2) and computing "
        "diversity indices.\n"
        "PREFIX\nMockLog: Converting the VCF file to a "
        "compressed BCF file...\n"
        "MockLog: Done.\n"
        "PREFIX\nMockLog: Indexing the BCF file...\n"
        "MockLog: Done indexing the BCF file.\n"
    )

    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_run_small_dataset_p(capsys):
    # "More integration tests are also nice too" -- Barack Obama
    #
    # A lot of test code here is reused from test_run_small_dataset_r() above.
    # Sorry this code is a bit gross.
    with tempfile.TemporaryDirectory() as td:
        run(
            FASTA,
            BAM,
            td,
            mock_log,
            True,
            min_p=20,
            div_index_p_list=[1, 50, 500, 5000],
            # This min_alt_pos value causes us to miss out on two of the
            # r-mutations found in the above integration test
            min_alt_pos=4,
            min_read_number=5,
        )

        # Check BCF
        bcf_obj, tt, tm = parse_sf_bcf(os.path.join(td, "naive-calls.bcf"))
        assert tt == "p"
        assert tm == 20
        exp_contigs = ["c1", "c1", "c3"]
        exp_pos = [11, 13, 7]
        exp_ref = ["G", "G", "A"]
        exp_mdp = [12, 12, 13]
        exp_alt = ["A", "A", "T"]
        exp_aad = [5, 5, 6]
        for ri, rec in enumerate(bcf_obj.fetch()):
            assert rec.contig == exp_contigs[ri]
            assert rec.pos == exp_pos[ri]
            assert rec.ref == exp_ref[ri]
            assert rec.info.get("MDP") == exp_mdp[ri]

            # Contingent upon the whole only-1-mut-at-a-position thing
            assert len(rec.alts) == 1
            assert rec.alts[0] == exp_alt[ri]
            rec_aad = rec.info.get("AAD")
            assert len(rec_aad) == 1
            assert rec_aad[0] == exp_aad[ri]

        # Check diversity index file
        obs_di_df = pd.read_csv(
            os.path.join(td, "diversity-indices.tsv"), sep="\t", index_col=0
        )
        c1_avg_cov = 12.0
        c2_avg_cov = ((11 * 10) + (10 * 2)) / 12
        c3_avg_cov = ((13 * 15) + 2) / 16
        exp_di_df = pd.DataFrame(
            {
                "AverageCoverage": [c1_avg_cov, c2_avg_cov, c3_avg_cov],
                "Length": [23, 12, 16],
                "DivIdx(p=1,minSuffCov=50000.0)": [np.nan, np.nan, np.nan],
                "DivIdx(p=50,minSuffCov=1000.0)": [np.nan, np.nan, np.nan],
                "DivIdx(p=500,minSuffCov=100.0)": [np.nan, np.nan, np.nan],
                "DivIdx(p=5000,minSuffCov=10.0)": [0.0, 0.0, 0.0],
            },
            index=pd.Index(["c1", "c2", "c3"], name="Contig"),
        )
        pd.testing.assert_frame_equal(obs_di_df, exp_di_df)

    exp_out = (
        "PREFIX\nMockLog: Loading and checking contig information...\n"
        "MockLog: The FASTA file describes 3 contigs.\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "references), with the same lengths.\n"
        "PREFIX\nMockLog: Creating and writing diversity index and VCF file "
        "headers...\n"
        "MockLog: Wrote out diversity index and VCF file headers.\n"
        "MockLog: (We'll convert the VCF file to an indexed BCF afterwards.)\n"
    )
    exp_out += (
        "PREFIX\nMockLog: Running p-mutation calling (--min-p = 0.20%) and "
        "computing diversity indices...\n"
        "MockLog: On contig c1 (23 bp) (1 / 3 = 33.33% done).\n"
        "MockLog: Called 2 p-mutation(s) (using --min-p = 0.20%) in contig c1."
        "\n"
        "MockLog: 1 / 4 diversity indices were defined for contig c1.\n"
        "MockLog: On contig c2 (12 bp) (2 / 3 = 66.67% done).\n"
        "MockLog: Called 0 p-mutation(s) (using --min-p = 0.20%) in contig c2."
        "\n"
        "MockLog: 1 / 4 diversity indices were defined for contig c2.\n"
        "MockLog: On contig c3 (16 bp) (3 / 3 = 100.00% done).\n"
        "MockLog: Called 1 p-mutation(s) (using --min-p = 0.20%) in contig c3."
        "\n"
        "MockLog: 1 / 4 diversity indices were defined for contig c3.\n"
    )
    exp_out += (
        "MockLog: Done running p-mutation calling (--min-p = 0.20%) and "
        "computing diversity indices.\n"
        "PREFIX\nMockLog: Converting the VCF file to a "
        "compressed BCF file...\n"
        "MockLog: Done.\n"
        "PREFIX\nMockLog: Indexing the BCF file...\n"
        "MockLog: Done indexing the BCF file.\n"
    )

    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_run_small_dataset_p_noverbose_diff_params(capsys):
    # Really, this is just me being nitpicky and trying to cover a few branches
    # that the other p-mutation integration test missed.
    with tempfile.TemporaryDirectory() as td:
        run(
            FASTA,
            BAM,
            td,
            mock_log,
            False,
            # Much larger min p
            min_p=3300,
            # Smaller min alt pos
            min_alt_pos=2,
            # Smaller min read number and a div idx p of 1000 (so now we
            # actually have some nonzero div indices)
            div_index_p_list=[1, 500, 1000, 5000],
            min_read_number=1,
        )

        # Check BCF
        bcf_obj, tt, tm = parse_sf_bcf(os.path.join(td, "naive-calls.bcf"))
        assert tt == "p"
        assert tm == 3300
        exp_contigs = ["c1", "c1", "c3"]
        exp_pos = [11, 13, 7]
        exp_ref = ["G", "G", "A"]
        exp_mdp = [12, 12, 13]
        exp_alt = ["A", "A", "T"]
        exp_aad = [5, 5, 6]
        for ri, rec in enumerate(bcf_obj.fetch()):
            assert rec.contig == exp_contigs[ri]
            assert rec.pos == exp_pos[ri]
            assert rec.ref == exp_ref[ri]
            assert rec.info.get("MDP") == exp_mdp[ri]
            assert len(rec.alts) == 1
            assert rec.alts[0] == exp_alt[ri]
            rec_aad = rec.info.get("AAD")
            assert len(rec_aad) == 1
            assert rec_aad[0] == exp_aad[ri]

        # Check diversity index file -- this'll be different due to us lowering
        # min_alt_pos and min_read_number
        obs_di_df = pd.read_csv(
            os.path.join(td, "diversity-indices.tsv"), sep="\t", index_col=0
        )
        c1_avg_cov = 12.0
        c2_avg_cov = ((11 * 10) + (10 * 2)) / 12
        c3_avg_cov = ((13 * 15) + 2) / 16
        exp_di_df = pd.DataFrame(
            {
                "AverageCoverage": [c1_avg_cov, c2_avg_cov, c3_avg_cov],
                "Length": [23, 12, 16],
                "DivIdx(p=1,minSuffCov=10000.0)": [np.nan, np.nan, np.nan],
                "DivIdx(p=500,minSuffCov=20.0)": [np.nan, np.nan, np.nan],
                # The denominator for c3 is 15 instead of 16, because one of
                # the positions in c3 has < 10x (mis)match coverage.
                "DivIdx(p=1000,minSuffCov=10.0)": [(3 / 23), 0.0, (2 / 15)],
                "DivIdx(p=5000,minSuffCov=2.0)": [0.0, 0.0, 0.0],
            },
            index=pd.Index(["c1", "c2", "c3"], name="Contig"),
        )
        pd.testing.assert_frame_equal(obs_di_df, exp_di_df)

    exp_out = (
        "PREFIX\nMockLog: Loading and checking contig information...\n"
        "MockLog: The FASTA file describes 3 contigs.\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "references), with the same lengths.\n"
        "PREFIX\nMockLog: Creating and writing diversity index and VCF file "
        "headers...\n"
        "MockLog: Wrote out diversity index and VCF file headers.\n"
        "MockLog: (We'll convert the VCF file to an indexed BCF afterwards.)\n"
    )
    # We passed verbose=False, so we skip a lot of the stuff in the middle
    exp_out += (
        "PREFIX\nMockLog: Running p-mutation calling (--min-p = 33.00%) and "
        "computing diversity indices...\n"
        "MockLog: Done running p-mutation calling (--min-p = 33.00%) and "
        "computing diversity indices.\n"
        "PREFIX\nMockLog: Converting the VCF file to a "
        "compressed BCF file...\n"
        "MockLog: Done.\n"
        "PREFIX\nMockLog: Indexing the BCF file...\n"
        "MockLog: Done indexing the BCF file.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_run_contig_with_no_alns_verbose_p(capsys):
    with tempfile.TemporaryDirectory() as td:
        run(
            os.path.join(IN_DIR, "contigs-with-c4.fasta"),
            os.path.join(IN_DIR, "c4-and-secondary.bam"),
            td,
            mock_log,
            True,
            min_p=3300,
            min_alt_pos=2,
            div_index_p_list=[1, 500, 1000, 5000],
            min_read_number=1,
        )

        # Check BCF
        bcf_obj, tt, tm = parse_sf_bcf(os.path.join(td, "naive-calls.bcf"))
        # Make sure that c4 is included in the BCF header, even though it has
        # no mutations called
        assert set(bcf_obj.header.contigs) == set(["c1", "c2", "c3", "c4"])
        # other stuff in the bcf -- same as before
        assert tt == "p"
        assert tm == 3300
        exp_contigs = ["c1", "c1", "c3"]
        exp_pos = [11, 13, 7]
        exp_ref = ["G", "G", "A"]
        exp_mdp = [12, 12, 13]
        exp_alt = ["A", "A", "T"]
        exp_aad = [5, 5, 6]
        for ri, rec in enumerate(bcf_obj.fetch()):
            assert rec.contig == exp_contigs[ri]
            assert rec.pos == exp_pos[ri]
            assert rec.ref == exp_ref[ri]
            assert rec.info.get("MDP") == exp_mdp[ri]
            assert len(rec.alts) == 1
            assert rec.alts[0] == exp_alt[ri]
            rec_aad = rec.info.get("AAD")
            assert len(rec_aad) == 1
            assert rec_aad[0] == exp_aad[ri]

        # Check diversity index file
        # c4 should be included here, but with completely undefined diversity
        # indices.
        obs_di_df = pd.read_csv(
            os.path.join(td, "diversity-indices.tsv"), sep="\t", index_col=0
        )
        c1_avg_cov = 12.0
        c2_avg_cov = ((11 * 10) + (10 * 2)) / 12
        c3_avg_cov = ((13 * 15) + 2) / 16
        exp_di_df = pd.DataFrame(
            {
                "AverageCoverage": [c1_avg_cov, c2_avg_cov, c3_avg_cov, 0],
                "Length": [23, 12, 16, 100],
                "DivIdx(p=1,minSuffCov=10000.0)": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "DivIdx(p=500,minSuffCov=20.0)": [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                "DivIdx(p=1000,minSuffCov=10.0)": [
                    (3 / 23),
                    0.0,
                    (2 / 15),
                    np.nan,
                ],
                "DivIdx(p=5000,minSuffCov=2.0)": [0.0, 0.0, 0.0, np.nan],
            },
            index=pd.Index(["c1", "c2", "c3", "c4"], name="Contig"),
        )
        pd.testing.assert_frame_equal(obs_di_df, exp_di_df)

    exp_out = (
        "PREFIX\nMockLog: Loading and checking contig information...\n"
        "MockLog: The FASTA file describes 4 contigs.\n"
        "MockLog: All of these are included in the BAM file (which has 4 "
        "references), with the same lengths.\n"
        "PREFIX\nMockLog: Creating and writing diversity index and VCF file "
        "headers...\n"
        "MockLog: Wrote out diversity index and VCF file headers.\n"
        "MockLog: (We'll convert the VCF file to an indexed BCF afterwards.)\n"
    )
    exp_out += (
        "PREFIX\nMockLog: Running p-mutation calling (--min-p = 33.00%) and "
        "computing diversity indices...\n"
        "MockLog: On contig c1 (23 bp) (1 / 4 = 25.00% done).\n"
        "MockLog: Called 2 p-mutation(s) (using --min-p = 33.00%) in contig "
        "c1.\n"
        "MockLog: 2 / 4 diversity indices were defined for contig c1.\n"
        "MockLog: On contig c2 (12 bp) (2 / 4 = 50.00% done).\n"
        "MockLog: Called 0 p-mutation(s) (using --min-p = 33.00%) in contig "
        "c2.\n"
        "MockLog: 2 / 4 diversity indices were defined for contig c2.\n"
        "MockLog: On contig c3 (16 bp) (3 / 4 = 75.00% done).\n"
        "MockLog: Called 1 p-mutation(s) (using --min-p = 33.00%) in contig "
        "c3.\n"
        "MockLog: 2 / 4 diversity indices were defined for contig c3.\n"
        "MockLog: On contig c4 (100 bp) (4 / 4 = 100.00% done).\n"
        "MockLog: No linear alignments to contig c4 exist in the BAM file. "
        "Not performing mutation calling within this contig.\n"
    )
    exp_out += (
        "MockLog: Done running p-mutation calling (--min-p = 33.00%) and "
        "computing diversity indices.\n"
        "PREFIX\nMockLog: Converting the VCF file to a "
        "compressed BCF file...\n"
        "MockLog: Done.\n"
        "PREFIX\nMockLog: Indexing the BCF file...\n"
        "MockLog: Done indexing the BCF file.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_run_contig_with_no_alns_nonverbose_r(capsys):
    with tempfile.TemporaryDirectory() as td:
        run(
            os.path.join(IN_DIR, "contigs-with-c4.fasta"),
            os.path.join(IN_DIR, "c4-and-secondary.bam"),
            td,
            mock_log,
            False,
            min_r=2,
            div_index_r_list=[1],
            min_cov_factor=2,
        )

        # Check BCF
        bcf_obj, tt, tm = parse_sf_bcf(os.path.join(td, "naive-calls.bcf"))
        assert set(bcf_obj.header.contigs) == set(["c1", "c2", "c3", "c4"])
        assert tt == "r"
        assert tm == 2

        exp_contigs = ["c1", "c1", "c1", "c3", "c3"]
        exp_pos = [4, 11, 13, 7, 8]
        exp_ref = ["G", "G", "G", "A", "T"]
        exp_mdp = [12, 12, 12, 13, 13]
        exp_alt = ["T", "A", "A", "T", "C"]
        exp_aad = [3, 5, 5, 6, 3]
        for ri, rec in enumerate(bcf_obj.fetch()):
            assert rec.contig == exp_contigs[ri]
            assert rec.pos == exp_pos[ri]
            assert rec.ref == exp_ref[ri]
            assert rec.info.get("MDP") == exp_mdp[ri]
            assert len(rec.alts) == 1
            assert rec.alts[0] == exp_alt[ri]
            rec_aad = rec.info.get("AAD")
            assert len(rec_aad) == 1
            assert rec_aad[0] == exp_aad[ri]

        # Check diversity index file
        obs_di_df = pd.read_csv(
            os.path.join(td, "diversity-indices.tsv"), sep="\t", index_col=0
        )
        c1_avg_cov = 12.0
        c2_avg_cov = ((11 * 10) + (10 * 2)) / 12
        c3_avg_cov = ((13 * 15) + 2) / 16
        exp_di_df = pd.DataFrame(
            {
                "AverageCoverage": [c1_avg_cov, c2_avg_cov, c3_avg_cov, 0],
                "Length": [23, 12, 16, 100],
                "DivIdx(r=1,minSuffCov=2)": [3 / 23, 0, 3 / 16, np.nan],
            },
            index=pd.Index(["c1", "c2", "c3", "c4"], name="Contig"),
        )
        pd.testing.assert_frame_equal(obs_di_df, exp_di_df)

    exp_out = (
        "PREFIX\nMockLog: Loading and checking contig information...\n"
        "MockLog: The FASTA file describes 4 contigs.\n"
        "MockLog: All of these are included in the BAM file (which has 4 "
        "references), with the same lengths.\n"
        "PREFIX\nMockLog: Creating and writing diversity index and VCF file "
        "headers...\n"
        "MockLog: Wrote out diversity index and VCF file headers.\n"
        "MockLog: (We'll convert the VCF file to an indexed BCF afterwards.)\n"
    )
    exp_out += (
        "PREFIX\nMockLog: Running r-mutation calling (--min-r = 2) and "
        "computing diversity indices...\n"
    )
    exp_out += (
        "MockLog: Done running r-mutation calling (--min-r = 2) and "
        "computing diversity indices.\n"
        "PREFIX\nMockLog: Converting the VCF file to a "
        "compressed BCF file...\n"
        "MockLog: Done.\n"
        "PREFIX\nMockLog: Indexing the BCF file...\n"
        "MockLog: Done indexing the BCF file.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_run_degen_nts():
    with tempfile.TemporaryDirectory() as td:
        with pytest.raises(SequencingDataError) as ei:
            run(
                FASTA,
                os.path.join(IN_DIR, "degen.bam"),
                td,
                mock_log,
                True,
                min_r=2,
                div_index_r_list=[1, 2, 3, 4, 5, 6, 7],
                min_cov_factor=2,
            )
        assert str(ei.value) == (
            "Found a degenerate nucleotide aligned to contig c1 at "
            "(1-indexed) position 9. Alignments including degenerate "
            "nucleotides (e.g. N) are not supported."
        )
