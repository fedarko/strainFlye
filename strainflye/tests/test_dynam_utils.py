import os
import pandas as pd
from pytest import approx
from strainflye import dynam_utils as du
from strainflye.tests.utils_for_testing import mock_log


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BAM = os.path.join(IN_DIR, "alignment.bam")


def test_skew():
    assert du.skew("ATATATATT") == 0
    assert du.skew("") == 0

    assert du.skew("CCCCC") == -1
    assert du.skew("GGGGG") == 1
    assert du.skew("GC") == 0
    assert du.skew("CG") == 0
    assert du.skew("GGGGC") == 3 / 5
    assert du.skew("CCCCG") == -3 / 5
    assert du.skew("CGCCC") == -3 / 5


def check_lists_approx_equal(exp, obs):
    assert len(exp) == len(obs)
    for (e, o) in zip(exp, obs):
        assert e == approx(o)


def test_update_cumulative_binned_skews():
    c = []
    du.update_cumulative_binned_skews(c, 1)
    assert c == [1]

    du.update_cumulative_binned_skews(c, -1)
    check_lists_approx_equal(c, [1, 0])

    du.update_cumulative_binned_skews(c, 0)
    check_lists_approx_equal(c, [1, 0, 0])

    du.update_cumulative_binned_skews(c, 0.5)
    check_lists_approx_equal(c, [1, 0, 0, 0.5])

    du.update_cumulative_binned_skews(c, -0.8)
    check_lists_approx_equal(c, [1, 0, 0, 0.5, -0.3])


def test_log_bin_ct_info(capsys):
    du.log_bin_ct_info("c1", 23, 5, mock_log)
    assert capsys.readouterr().out == (
        "MockLog: Creating 4 bins of length 5 bp and 1 smaller bin of length "
        "3 bp for contig c1...\n"
    )

    du.log_bin_ct_info("c1", 23, 500, mock_log)
    assert capsys.readouterr().out == (
        "MockLog: Creating 1 smaller bin of length 23 bp for contig c1...\n"
    )

    du.log_bin_ct_info("c1", 500, 10, mock_log)
    assert capsys.readouterr().out == (
        "MockLog: Creating 50 bins of length 10 bp for contig c1...\n"
    )

    # of course we make it say "bin" instead of "bins"
    du.log_bin_ct_info("c1", 1234, 1234, mock_log)
    assert capsys.readouterr().out == (
        "MockLog: Creating 1 bin of length 1,234 bp for contig c1...\n"
    )


def test_contig_covskew_c1_binlen10():
    nbcovs, cbskews, lp, cp = du.contig_covskew("c1", FASTA, BAM, 10, 0.7, 1.3)
    assert lp == [1, 11, 21]
    assert cp == [(1 + 10) / 2, (11 + 20) / 2, (21 + 23) / 2]

    # c1 has uniform coverage
    assert nbcovs == [1, 1, 1]

    # c1's sequence is ACTGACACCCAAACCAAACCTAC. Using bins of length 10:
    #
    # Bin 1: ACTGACACCC (5 C, 1 G): Skew = (1-5) / 6 = -4/6 = -0.6666...
    # Bin 2: AAACCAAACC (4 C, 0 G): Skew = (0-4) / 4 = -1
    # Bin 3: TAC        (1 C, 0 G): Skew = (0-1) / 1 = -1
    check_lists_approx_equal(
        cbskews, [-4 / 6, (-4 / 6) + (-1), (-4 / 6) + (-2)]
    )


def test_contig_covskew_c2_binlen2_lowclamp():
    nbcovs, cbskews, lp, cp = du.contig_covskew(
        "c2", FASTA, BAM, 2, 0.95, 1.05
    )
    assert lp == [1, 3, 5, 7, 9, 11]
    check_lists_approx_equal(cp, [1.5, 3.5, 5.5, 7.5, 9.5, 11.5])

    # First bin has median coverage 10; remaining bins have median coverage 11.
    # So, median of medians is 11. The first bin's normalized coverage is thus
    # 10 / 11 = 0.909090..., and -- since we've set the clamp range to [0.95,
    # 1.05] -- we should clamp this coverage to 0.95.
    check_lists_approx_equal(nbcovs, [0.95, 1, 1, 1, 1, 1])

    # c2's sequence is AAAAAAGGGGGG. Using bins of length 2:
    #
    # Bin 1: AA: Skew = 0
    # Bin 2: AA: Skew = 0
    # Bin 3: AA: Skew = 0
    # Bin 4: GG: Skew = 1
    # Bin 5: GG: Skew = 1
    # Bin 6: GG: Skew = 1
    check_lists_approx_equal(cbskews, [0, 0, 0, 1, 2, 3])


def test_contig_covskew_c3_binlen1_upperclamp():
    nbcovs, cbskews, lp, cp = du.contig_covskew("c3", FASTA, BAM, 1, 0.7, 0.8)

    just_1idxed_positions = list(range(1, 17))
    assert lp == just_1idxed_positions
    assert cp == just_1idxed_positions

    # all positions but the last have coverage 15; the last position has
    # coverage 2. so, the median coverage is 15, and the last position's
    # normalized coverage is clamped to the lower bound. All other positions
    # have a normalized coverage of 1, which gets clamped to the silly upper
    # bound of 0.8 I set here. (This is just for testing -- in practice, the
    # upper bound should always be > 1.)
    check_lists_approx_equal(nbcovs, ([0.8] * 15) + [0.7])

    # c3's sequence is AAAAAAGGGGGG. Using bins of length 2:
    #
    # Bin 1: AA: Skew = 0
    # Bin 2: AA: Skew = 0
    # Bin 3: AA: Skew = 0
    # Bin 4: GG: Skew = 1
    # Bin 5: GG: Skew = 1
    # Bin 6: GG: Skew = 1
    assert cbskews == [0] * 16


def test_run_covskew_good_noverbose_binlen10(capsys, tmp_path):
    du.run_covskew(FASTA, BAM, 10, 0.3, tmp_path, False, mock_log)

    assert sorted(os.listdir(tmp_path)) == [
        "c1_covskew.tsv",
        "c2_covskew.tsv",
        "c3_covskew.tsv",
    ]

    # this one should be the same as test_contig_covskew_c1_binlen10() above
    c1 = pd.read_csv(tmp_path / "c1_covskew.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        c1,
        pd.DataFrame(
            {
                "LeftPos_1IndexedInclusive": [1, 11, 21],
                "CenterPos": [(1 + 10) / 2, (11 + 20) / 2, (21 + 23) / 2],
                "NormalizedCoverage": [1.0, 1.0, 1.0],
                "CumulativeSkew": [-4 / 6, (-4 / 6) + (-1), (-4 / 6) + (-2)],
            }
        ),
    )

    c2 = pd.read_csv(tmp_path / "c2_covskew.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        c2,
        pd.DataFrame(
            {
                "LeftPos_1IndexedInclusive": [1, 11],
                "CenterPos": [(1 + 10) / 2, 11.5],
                "NormalizedCoverage": [1.0, 1.0],
                # both bins only include Gs, so they both have skews of 1.
                "CumulativeSkew": [1.0, 2.0],
            }
        ),
    )

    c3 = pd.read_csv(tmp_path / "c3_covskew.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        c3,
        pd.DataFrame(
            {
                "LeftPos_1IndexedInclusive": [1, 11],
                "CenterPos": [(1 + 10) / 2, (11 + 16) / 2],
                "NormalizedCoverage": [1.0, 1.0],
                "CumulativeSkew": [0, 0],
            }
        ),
    )

    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Loading and checking FASTA and BAM files...\n"
        "MockLog: The FASTA file describes 3 contig(s).\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "reference(s)), with the same lengths.\n"
        "PREFIX\nMockLog: Going through contigs and computing coverage/skew "
        "information...\n"
        "MockLog: Done.\n"
    )


def test_run_covskew_good_verbose_binlen1000(capsys, tmp_path):
    du.run_covskew(FASTA, BAM, 1000, 0.5, tmp_path, True, mock_log)

    assert sorted(os.listdir(tmp_path)) == [
        "c1_covskew.tsv",
        "c2_covskew.tsv",
        "c3_covskew.tsv",
    ]

    c1 = pd.read_csv(tmp_path / "c1_covskew.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        c1,
        pd.DataFrame(
            {
                "LeftPos_1IndexedInclusive": [1],
                "CenterPos": [(1 + 23) / 2],
                "NormalizedCoverage": [1.0],
                "CumulativeSkew": [(1 - 10) / 11],
            }
        ),
    )

    c2 = pd.read_csv(tmp_path / "c2_covskew.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        c2,
        pd.DataFrame(
            {
                "LeftPos_1IndexedInclusive": [1],
                "CenterPos": [(1 + 12) / 2],
                "NormalizedCoverage": [1.0],
                "CumulativeSkew": [1.0],
            }
        ),
    )

    c3 = pd.read_csv(tmp_path / "c3_covskew.tsv", sep="\t")
    pd.testing.assert_frame_equal(
        c3,
        pd.DataFrame(
            {
                "LeftPos_1IndexedInclusive": [1],
                "CenterPos": [(1 + 16) / 2],
                "NormalizedCoverage": [1.0],
                "CumulativeSkew": [0],
            }
        ),
    )

    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Loading and checking FASTA and BAM files...\n"
        "MockLog: The FASTA file describes 3 contig(s).\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "reference(s)), with the same lengths.\n"
        "PREFIX\nMockLog: Going through contigs and computing coverage/skew "
        "information...\n"
        "MockLog: On contig c1 (23 bp) (1 / 3 contigs = 33.33%).\n"
        "MockLog: Creating 1 smaller bin of length 23 bp for contig c1...\n"
        "MockLog: On contig c2 (12 bp) (2 / 3 contigs = 66.67%).\n"
        "MockLog: Creating 1 smaller bin of length 12 bp for contig c2...\n"
        "MockLog: On contig c3 (16 bp) (3 / 3 contigs = 100.00%).\n"
        "MockLog: Creating 1 smaller bin of length 16 bp for contig c3...\n"
        "MockLog: Done.\n"
    )
