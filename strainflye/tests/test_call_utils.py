import pytest
from strainflye.call_utils import call_r_mutation
from strainflye.errors import ParameterError


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
