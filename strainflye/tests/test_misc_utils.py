import os
import shutil
import tempfile
import pytest
import strainflye.misc_utils as mu
from strainflye.errors import ParameterError


def test_make_output_dir_exists():
    with tempfile.TemporaryDirectory() as tdname:
        assert os.path.exists(tdname)
        child_file_name = os.path.join(tdname, "test_file.txt")
        with open(child_file_name, "w") as f:
            f.write("please don't overwrite me!")

        mu.make_output_dir(tdname)

        # Check that this didn't overwrite anything (the directory itself, or
        # any of its contents)
        assert os.path.exists(tdname)
        assert os.path.exists(child_file_name)
        with open(child_file_name, "r") as f:
            ctext = f.read()
        assert ctext == "please don't overwrite me!"


def test_make_output_dir_notexists_multilevel():
    tdname_top = os.path.join(tempfile.gettempdir(), "_strainflye_td1")
    tdname_bot = os.path.join(tdname_top, "_strainflye_td2")
    try:
        assert not os.path.exists(tdname_top)
        assert not os.path.exists(tdname_bot)
        mu.make_output_dir(tdname_bot)
        assert os.path.exists(tdname_top)
        assert os.path.exists(tdname_bot)
    finally:
        # https://stackoverflow.com/a/13118112
        shutil.rmtree(tdname_top)


def test_verify_contigs_subset_raises_error():
    with pytest.raises(ParameterError) as ei:
        mu.verify_contigs_subset(set("abcdef"), set("abdef"), "s1", "s2")
    assert str(ei.value) == "All contigs in s1 must also be contained in s2."


def test_verify_contigs_subset_good():
    # Sets are identical
    mu.verify_contigs_subset(set("abcdef"), set("abcdef"), "s1", "s2")
    mu.verify_contigs_subset(set("a"), set("a"), "s1", "s2")

    # Proper subset
    mu.verify_contigs_subset(set(""), set(""), "s1", "s2")
    mu.verify_contigs_subset(set(""), set("ab"), "s1", "s2")
    mu.verify_contigs_subset(set("a"), set("ab"), "s1", "s2")
    mu.verify_contigs_subset(set("ab"), set("abc"), "s1", "s2")


def test_verify_contigs_subset_exact():
    # Check that exact doesn't change mandate that child is subset of parent
    with pytest.raises(ParameterError) as ei:
        mu.verify_contigs_subset(
            set("abcdef"), set("abdef"), "s1", "s2", exact=True
        )
    assert str(ei.value) == "All contigs in s1 must also be contained in s2."

    # Now, check that exact ensures that the two sets are equal, even if child
    # is a subset of parent
    with pytest.raises(ParameterError) as ei:
        mu.verify_contigs_subset(
            set("abcd"), set("abcdef"), "s1", "s2", exact=True
        )
    assert str(ei.value) == "All contigs in s2 must also be contained in s1."

    # Now, check that exact succeeds
    mu.verify_contigs_subset(set("abcd"), set("abcd"), "s1", "s2", exact=True)
    mu.verify_contigs_subset(set(""), set(""), "s1", "s2", exact=True)
