from pytest import approx
from strainflye import dynam_utils as du
from strainflye.tests.utils_for_testing import mock_log


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
