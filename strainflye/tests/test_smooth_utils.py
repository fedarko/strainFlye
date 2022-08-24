import os
import stat
import tempfile
import pytest
import strainflye.smooth_utils as su
from strainflye.errors import ParameterError
from strainflye.tests.utils_for_testing import mock_log


def test_convert_to_runs():
    assert su.convert_to_runs([]) == []
    assert su.convert_to_runs([3]) == [(3, 3)]
    assert su.convert_to_runs([3, 4, 9, 10, 12, 22]) == [
        (3, 4),
        (9, 10),
        (12, 12),
        (22, 22),
    ]
    assert su.convert_to_runs(
        [1, 2, 3, 5, 6, 7, 10, 20, 31, 32, 33, 34, 35, 40]
    ) == [(1, 3), (5, 7), (10, 10), (20, 20), (31, 35), (40, 40)]


def test_find_lja_bin_direct(capsys):
    # We don't validate the bin, we just assume that it exists and is
    # executable (in practice, click should verify these things, although ofc
    # we're vulnerable to race conditions so if the user wants to yank the
    # binary out from under us they're welcome to do that -- it'll just crash
    # stuff)
    assert su.find_lja_bin("asdf/asdf", mock_log) == "asdf/asdf"
    assert capsys.readouterr().out == ""

    # TODO: test directly by adding to $PATH temporarily? hm


def test_find_lja_bin_pathsearch_notfound(capsys):
    # This test assumes that the system on which this is running doesn't have
    # "lja" added to the PATH. If your computer has LJA in the PATH, then this
    # will cause this test to fail. (It shouldn't be a problem unless GitHub
    # Actions' machines start having LJA pre-installed by default, in which
    # case we've got a problem...)
    with pytest.raises(ParameterError) as ei:
        su.find_lja_bin(None, mock_log)
    assert str(ei.value) == (
        '--lja-bin was not specified, and we couldn\'t find "lja" in any of '
        "the system-wide executable locations given by $PATH."
    )
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Since --lja-bin wasn't specified, looking in $PATH "
        'for "lja"...\n'
    )


def test_find_lja_bin_pathsearch_found(capsys):
    # Rather than actually install LJA, we just make a fake LJA binary file and
    # modify $PATH to point to the location of this file.
    with tempfile.TemporaryDirectory() as fake_lja_dir:

        # Make the (cough cough) very legitimate real binary file
        fake_lja_bin_loc = os.path.join(fake_lja_dir, "lja")
        with open(fake_lja_bin_loc, "w") as f:
            f.write(
                "#! /usr/bin/env bash\n"
                'echo "hemblo i am a real assembler blease run me on your '
                'compuper"\n'
            )

        # Make the file executable, so that it turns up when we run
        # shutil.which(): see
        # https://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python#comment26692909_12792002
        # I guess this process is analogous to running "chmod +x" on the file.
        #
        # In the ~9 years since that Stack Overflow post (oofa doofa), it seems
        # like it is no longer kosher to OR the result of os.stat() with other
        # "mode bits." However, it looks like we can just use os.stat().st_mode
        # to get the current mode bits for the file.
        new_mode_bits = (
            os.stat(fake_lja_bin_loc).st_mode
            | stat.S_IXGRP
            | stat.S_IXUSR
            | stat.S_IXOTH
        )
        os.chmod(fake_lja_bin_loc, new_mode_bits)
        # Add fake_lja_dir to the PATH environment variable. os.pathsep is
        # probably a ":" character on most systems; see
        # https://stackoverflow.com/a/1681244.
        os.environ["PATH"] += os.pathsep + fake_lja_dir

        # Note that this change shouldn't persist after these changes,
        # thankfully; see https://stackoverflow.com/a/66390638.

        assert su.find_lja_bin(None, mock_log) == fake_lja_bin_loc
        assert capsys.readouterr().out == (
            "PREFIX\nMockLog: Since --lja-bin wasn't specified, looking in "
            '$PATH for "lja"...\n'
            f"MockLog: Found it at {fake_lja_bin_loc}!\n"
        )


def test_verify_vrf2_good(capsys):
    su.verify_vrf2({"D": 101, "E": 10000}, 50, mock_log)
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: All contigs must be > (2 \u00d7 "
        "--virtual-read-flank) = 100 bp long. Checking this...\n"
        "MockLog: All contigs meet this minimum length.\n"
    )


def test_verify_vrf2_bad():
    with pytest.raises(ParameterError) as ei:
        su.verify_vrf2({"C": 100, "D": 101, "E": 10000}, 50, mock_log)
    assert str(ei.value) == (
        "Contig C is 100 bp, which is \u2264 100 bp. Depending on your "
        "goals, you may want to remove short contigs from this FASTA file, "
        "lower --virtual-read-flank, or set --no-virtual-reads."
    )
