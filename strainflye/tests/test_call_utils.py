import pytest
from strainflye.call_utils import (
    p2filter,
    r2filter,
    get_r_increments,
    get_p_increments,
    get_alt_pos_info,
    call_r_mutations,
    call_p_mutations,
)
from strainflye.errors import ParameterError
from pytest import approx


def check_approx_equal(obs, exp):
    for op, ep in zip(obs, exp):
        assert op == approx(ep)


def test_get_r_increments_d1():
    rv = get_r_increments(1, 5, 1, print)
    assert rv == [1, 2, 3, 4, 5]

    rv = get_r_increments(1, 2, 1, print)
    assert rv == [1, 2]

    rv = get_r_increments(5, 10, 1, print)
    assert rv == [5, 6, 7, 8, 9, 10]


def test_get_r_increments_d2():
    rv = get_r_increments(1, 5, 2, print)
    assert rv == [1, 3, 5]


def test_get_r_increments_one_r_due_to_large_delta(capsys):
    rv = get_r_increments(5, 50, 100, print)
    captured = capsys.readouterr()
    exp_out = (
        "Computing r-mutations for 1 value of r: 5.\n"
        "Warning: --max-r = 50 will not be included in the r-values used. "
        "This is due to --max-r minus --min-r not being divisible by "
        "--delta-r: 50 - 5 = 45, and 45 mod 100 = 45.\n"
    )
    assert captured.out == exp_out
    assert rv == [5]


def test_get_r_increments_nomax(capsys):
    rv = get_r_increments(5, 50, 2, print)
    captured = capsys.readouterr()
    exp_out = (
        "Computing r-mutations for 23 values of r.\n"
        "Warning: --max-r = 50 will not be included in the r-values used. "
        "This is due to --max-r minus --min-r not being divisible by "
        "--delta-r: 50 - 5 = 45, and 45 mod 2 = 1.\n"
    )
    assert captured.out == exp_out
    assert rv == [
        5,
        7,
        9,
        11,
        13,
        15,
        17,
        19,
        21,
        23,
        25,
        27,
        29,
        31,
        33,
        35,
        37,
        39,
        41,
        43,
        45,
        47,
        49,
    ]


def test_get_r_increments_bad():
    with pytest.raises(ParameterError) as errorinfo:
        get_r_increments(5, 2, 1, print)
    assert "Minimum r must be less than maximum r." == str(errorinfo.value)


def test_get_p_increments_bad():
    with pytest.raises(ParameterError) as errorinfo:
        get_p_increments(5, 2, 1, print)
    assert "Minimum p must be less than maximum p." == str(errorinfo.value)


def test_get_p_increments(capsys):
    pv = get_p_increments(1, 5, 1, print)
    check_approx_equal(pv, [0.01, 0.02, 0.03, 0.04, 0.05])
    captured = capsys.readouterr()
    exp_out = "Computing p-mutations for 5 values of p.\n"
    assert captured.out == exp_out


def test_get_p_increments_one_p_due_to_large_delta(capsys):
    pv = get_p_increments(2, 5, 45, print)
    check_approx_equal(pv, [0.02])
    captured = capsys.readouterr()
    exp_out = (
        "Computing p-mutations for 1 value of p: 0.02%.\n"
        "Warning: --max-p = 5 (aka 0.05%) will not be included in the values "
        "of p used. This is due to --max-p minus --min-p not being divisible "
        "by --delta-p: 5 - 2 = 3, and 3 mod 45 = 3.\n"
    )
    assert captured.out == exp_out


def test_get_p_increments_sf_paper_fig2(capsys):
    obs_pv = get_p_increments(15, 200, 1, print)
    exp_pv = [p / 100 for p in range(15, 201)]
    check_approx_equal(obs_pv, exp_pv)
    captured = capsys.readouterr()
    exp_out = "Computing p-mutations for 186 values of p.\n"
    assert captured.out == exp_out


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


def test_call_r_mutations():
    # paranoia
    assert call_r_mutations(5, [1, 2, 3, 4, 5, 6, 7, 8]) == (
        [True, True, True, True, True, False, False, False],
        True,
    )
    assert call_r_mutations(5, [4]) == ([True], True)
    assert call_r_mutations(5, [6]) == ([False], False)
    assert call_r_mutations(5, [100, 101]) == ([False, False], False)
    assert call_r_mutations(100, [5, 6]) == ([True, True], True)
    assert call_r_mutations(1, [1]) == ([True], True)


def test_call_p_mutations():
    # TODO add more stuff here -- e.g. more detailed floating-point tests
    assert call_p_mutations(5, 10, [50], 2) == ([True], True)
    assert call_p_mutations(2, 10, [50], 2) == ([False], False)
    assert call_p_mutations(2, 10, [1, 2, 3, 5, 10, 15, 18, 20, 25], 2) == (
        [True, True, True, True, True, True, True, True, False],
        True,
    )
    assert call_p_mutations(
        5, 1000, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1, 2, 3], 2
    ) == (
        [True, True, True, True, True, False, False, False, False],
        True,
    )


def test_p2filter():
    # should always use exactly two points of precision
    # and since we have the min, max, and delta p parameters as integers with
    # min 0.01, this works
    assert p2filter(0.01) == "p0.01"
    assert p2filter(0.010000003) == "p0.01"
    assert p2filter(3) == "p3.00"
    assert p2filter(1.5) == "p1.50"
    assert p2filter(1.53219) == "p1.53"
    assert p2filter(1.53819) == "p1.54"


def test_r2filter():
    assert r2filter(1) == "r1"
    assert r2filter(1000) == "r1000"
    assert r2filter(10000) == "r10000"
    assert r2filter(35) == "r35"
