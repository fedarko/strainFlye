import pytest
from strainflye.errors import ParameterError
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
from .utils_for_testing import mock_log


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
