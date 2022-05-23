import pytest
from strainflye.call_utils import (
    get_alt_pos_info,
    get_pos_info_str,
    call_r_mutation,
    call_p_mutation,
)
from strainflye.errors import ParameterError
from pytest import approx


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
