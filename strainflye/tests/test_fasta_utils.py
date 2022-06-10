import re
import pytest
import strainflye.fasta_utils as utils
from io import StringIO
from strainflye.errors import SequencingDataError, ParameterError


# ... In practice, the user will see "... in asdf.fasta ..." in error messages
# instead of this nonsense representing the "name" of a StringIO
SIO_REPR = "<_io.StringIO object at [a-zA-Z0-9]*>"


def test_get_name2len_basic():
    sio = StringIO(">asdf\nAACCCAAATGA\n>ghij\nC\n")
    n2l = utils.get_name2len(sio)
    assert n2l == {"asdf": 11, "ghij": 1}


def test_get_name2len_oneseq():
    sio = StringIO(">lonely\nGTAC\n")
    n2l = utils.get_name2len(sio, min_num_contigs=1)
    assert n2l == {"lonely": 4}


def test_get_name2len_multiline():
    sio = StringIO(">longboi\nGTAC\nCCT\nTA\n>hi\nCAAT")
    n2l = utils.get_name2len(sio)
    assert n2l == {"longboi": 9, "hi": 4}


def test_get_name2len_dups():
    sio = StringIO(">seq1\nACGT\n>seq1\nTA\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)

    exp_pattern = f"Duplicate sequence name in {SIO_REPR}: seq1"
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_name2len_blankheader():
    # gotta escape the parens due to regex stuff. Using \(, \) made flake8 (?)
    # yell at me so I'm using \\(, \\) instead.
    exp_pattern = (
        "There exists a blank \\(or just whitespace\\) sequence name in "
        f"{SIO_REPR}"
    )

    # Test 1: has some whitepsace in header (skbio strips it)
    sio = StringIO(">seq1\nACGT\n>  \t  \nTA\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)
    assert re.match(exp_pattern, str(ei.value)) is not None

    # Test 2: completely blank
    sio = StringIO(">seq1\nACGT\n>\nTA\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_name2len_degen():
    sio = StringIO(">seq0\nGTAC\n>seq1\nACMGT\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)
    exp_pattern = (
        f"Sequence seq1 in {SIO_REPR} has at least one degenerate nucleotide. "
        "This isn't supported at the moment, sorry."
    )
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_name2len_gap():
    sio = StringIO(">seq0\nGTA-C\n>seq1\nACGT\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)
    exp_pattern = (
        f"Sequence seq0 in {SIO_REPR} has at least one gap. "
        "This isn't supported at the moment, sorry."
    )
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_name2len_min_num_contigs():
    sio = StringIO(">lonely\nGTAC\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)
    exp_pattern = f"Less than 2 contigs are given in {SIO_REPR}."
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_get_name2len_min_num_contigs_zero():
    # I assumed that this case would cause skbio to raise an error, but it just
    # raises a warning. Fortunately, the min_num_contigs check automatically
    # accounts for this case.
    sio = StringIO("")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_name2len(sio)
    exp_pattern = f"Less than 2 contigs are given in {SIO_REPR}."
    assert re.match(exp_pattern, str(ei.value)) is not None


def test_verify_contigs_subset_raises_error():
    with pytest.raises(ParameterError) as ei:
        utils.verify_contigs_subset(set("abcdef"), set("abdef"), "s1", "s2")
    assert str(ei.value) == "All contigs in s1 must also be contained in s2."


def test_verify_contigs_subset_good():
    # Sets are identical
    utils.verify_contigs_subset(set("abcdef"), set("abcdef"), "s1", "s2")
    utils.verify_contigs_subset(set("a"), set("a"), "s1", "s2")

    # Proper subset
    utils.verify_contigs_subset(set(""), set(""), "s1", "s2")
    utils.verify_contigs_subset(set(""), set("ab"), "s1", "s2")
    utils.verify_contigs_subset(set("a"), set("ab"), "s1", "s2")
    utils.verify_contigs_subset(set("ab"), set("abc"), "s1", "s2")


def test_verify_contigs_subset_exact():
    # Check that exact doesn't change mandate that child is subset of parent
    with pytest.raises(ParameterError) as ei:
        utils.verify_contigs_subset(
            set("abcdef"), set("abdef"), "s1", "s2", exact=True
        )
    assert str(ei.value) == "All contigs in s1 must also be contained in s2."

    # Now, check that exact ensures that the two sets are equal, even if child
    # is a subset of parent
    with pytest.raises(ParameterError) as ei:
        utils.verify_contigs_subset(
            set("abcd"), set("abcdef"), "s1", "s2", exact=True
        )
    assert str(ei.value) == "All contigs in s2 must also be contained in s1."

    # Now, check that exact succeeds
    utils.verify_contigs_subset(
        set("abcd"), set("abcd"), "s1", "s2", exact=True
    )
    utils.verify_contigs_subset(set(""), set(""), "s1", "s2", exact=True)


def test_get_single_seq_basic():
    sio = StringIO(">asdf\nAACCCAAATGA\n>ghij\nC\n")
    seq = utils.get_single_seq(sio, "ghij")
    assert seq.metadata["id"] == "ghij"
    assert len(seq) == 1
    assert str(seq) == "C"


def test_get_single_seq_lonely():
    sio = StringIO(">asdf\nAACCCAAATGA\n")
    seq = utils.get_single_seq(sio, "asdf")
    assert seq.metadata["id"] == "asdf"
    assert len(seq) == 11
    assert str(seq) == "AACCCAAATGA"


def test_get_single_seq_missing():
    sio = StringIO(">asdf\nAACCCAAATGA\n")
    with pytest.raises(SequencingDataError) as ei:
        utils.get_single_seq(sio, "ghij")
    exp_pattern = f"Contig ghij is not in {SIO_REPR}."
    assert re.match(exp_pattern, str(ei.value)) is not None
