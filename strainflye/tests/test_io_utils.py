import re
import pytest
import strainflye.io_utils as iu
from io import StringIO
from strainflye.errors import SequencingDataError


# ... In practice, the user will see "... in asdf.fasta ..." in error messages
# instead of this nonsense representing the "name" of a StringIO
SIO_REPR = "<_io.StringIO object at [a-zA-Z0-9]*>"


def test_get_fasta_name2len_good():
    sio = StringIO(">asdf\nAACCCAAATGA\n>ghij\nC\n")
    n2l = iu.get_fasta_name2len(sio)
    assert n2l == {"asdf": 11, "ghij": 1}

    sio = StringIO(">lonely\nGTAC\n")
    n2l = iu.get_fasta_name2len(sio)
    assert n2l == {"lonely": 4}


def test_get_fasta_name2len_dups():
    sio = StringIO(">seq1\nACGT\n>seq1\nTA\n")
    with pytest.raises(SequencingDataError) as ei:
        iu.get_fasta_name2len(sio)

    exp_pattern = f"Duplicate sequence name in {SIO_REPR}: seq1"
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_fasta_name2len_blankheader():
    # gotta escape the parens due to regex stuff. Using \(, \) made flake8 (?)
    # yell at me so I'm using \\(, \\) instead.
    exp_pattern = (
        "There exists a blank \\(or just whitespace\\) sequence name in "
        f"{SIO_REPR}"
    )

    # Test 1: has some whitepsace in header (skbio strips it)
    sio = StringIO(">seq1\nACGT\n>  \t  \nTA\n")
    with pytest.raises(SequencingDataError) as ei:
        iu.get_fasta_name2len(sio)
    assert re.match(exp_pattern, str(ei.value)) is not None

    # Test 2: completely blank
    sio = StringIO(">seq1\nACGT\n>\nTA\n")
    with pytest.raises(SequencingDataError) as ei:
        iu.get_fasta_name2len(sio)
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_fasta_name2len_degen():
    sio = StringIO(">seq0\nGTAC\n>seq1\nACMGT\n")
    with pytest.raises(SequencingDataError) as ei:
        iu.get_fasta_name2len(sio)
    exp_pattern = (
        f"Sequence seq1 in {SIO_REPR} has at least one degenerate nucleotide. "
        "This isn't supported at the moment, sorry."
    )
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_fasta_name2len_gap():
    sio = StringIO(">seq0\nGTA-C\n>seq1\nACGT\n")
    with pytest.raises(SequencingDataError) as ei:
        iu.get_fasta_name2len(sio)
    exp_pattern = (
        f"Sequence seq0 in {SIO_REPR} has at least one gap. "
        "This isn't supported at the moment, sorry."
    )
    assert re.match(exp_pattern, str(ei.value)) is not None
