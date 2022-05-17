import pytest
from strainflye.call_utils import (
    get_alt_pos_info,
    call_r_mutation,
    call_p_mutation,
)
from strainflye.errors import ParameterError


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
    assert call_r_mutation(5, 5)
    assert call_r_mutation(5, 4)
    assert not call_r_mutation(5, 6)
    assert not call_r_mutation(5, 100)
    assert call_r_mutation(100, 5)
    assert call_r_mutation(1, 1)


def test_call_r_mutation_rangecheck():
    with pytest.raises(ParameterError) as errorinfo:
        call_r_mutation(1, 0)
    assert "r must be > 0" == str(errorinfo.value)


def test_call_p_mutation():
    # TODO add more stuff here
    assert call_p_mutation(5, 10, 50, 2)
    assert not call_p_mutation(2, 10, 50, 2)


def test_call_p_mutation_rangecheck():
    for bad_p in (-140, -100, -51, -1, -0.5, 0, 50.1, 51, 100, 140):
        with pytest.raises(ParameterError) as errorinfo:
            call_p_mutation(5, 10, bad_p, 2)
        assert "p must be in the range (0, 50]" == str(errorinfo.value)

    for bad_map in (-150, -100, -50, -1, -0.5):
        with pytest.raises(ParameterError) as errorinfo:
            call_p_mutation(5, 10, 2, bad_map)
        assert "--min-alt-pos must be >= 0" == str(errorinfo.value)
