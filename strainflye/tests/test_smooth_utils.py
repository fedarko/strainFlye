import os
import stat
import tempfile
import pytest
import pysam
import strainflye.smooth_utils as su
from strainflye.errors import ParameterError, WeirdError
from strainflye.tests.utils_for_testing import mock_log


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BAM = os.path.join(IN_DIR, "alignment.bam")


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


def fetch_specific_aln(contig, aln_seq, alnfile="alignment.bam"):
    # or, more accurately, an alignment with a particular sequence.
    bf = pysam.AlignmentFile(os.path.join(IN_DIR, alnfile))
    found_aln = False
    for aln in bf.fetch(contig):
        if aln.query_sequence == aln_seq:
            found_aln = True
            break
    # This should never happen, but let's check
    if not found_aln:
        raise EnvironmentError(
            f"There aren't any linear alns to {contig} in our test BAM file "
            f"that match the exact sequence we're looking for ({aln_seq})."
        )
    return aln


def test_get_smooth_aln_replacements_good():
    aln = fetch_specific_aln("c1", "ACTGACACCCAAACCAAACCTAC")
    mp = [3, 10, 12]
    mp2ra = {3: ("G", "T"), 10: ("G", "A"), 12: ("G", "A")}
    repls = su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert repls == {3: "G", 10: "A", 12: "A"}


def test_get_smooth_aln_replacements_all_alts():
    # still works, even if aln has the "alt" at all mutated positions.
    aln = fetch_specific_aln("c1", "ACTGACACCCAAACCAAACCTAC")
    mp = [3, 10, 12]
    mp2ra = {3: ("T", "G"), 10: ("G", "A"), 12: ("G", "A")}
    repls = su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert repls == {3: "G", 10: "A", 12: "A"}


def test_get_smooth_aln_replacements_not_ref_or_alt():
    # Although aln has a G at 0-indexed position 3, we change mp2ra so that the
    # ref here is now C. The fact that G is no longer the ref or the alt here
    # means that we should ignore aln, and not generate a smoothed read from
    # it.
    aln = fetch_specific_aln("c1", "ACTGACACCCAAACCAAACCTAC")
    mp = [3, 10, 12]
    mp2ra = {3: ("C", "T"), 10: ("G", "A"), 12: ("G", "A")}
    repls = su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert repls is None


def test_get_smooth_aln_replacements_deletion_in_aln():
    aln = fetch_specific_aln("c3", "TTTTTTTTTTTTTTT")

    # If the deletion isn't aligned to a mutated position, things are fine
    for second_mp_pos in range(7, 15):
        mp = [6, second_mp_pos]
        mp2ra = {6: ("A", "T"), second_mp_pos: ("T", "C")}
        repls = su.get_smooth_aln_replacements(aln, mp, mp2ra)
        assert repls == {6: "T", second_mp_pos: "T"}

    # However, the deletion is aligned to a mutated position, we have to ignore
    # this aln
    mp = [6, 15]
    mp2ra = {6: ("A", "T"), 15: ("T", "C")}
    repls = su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert repls is None


def test_get_smooth_aln_replacements_insertion_in_aln():
    aln = fetch_specific_aln(
        "c1", "ACTTACACCCAAACCTTTAAACCTAC", alnfile="insertion.bam"
    )
    # ensure that insertions are ignored -- the TTT in the middle of this
    # linear alignment is an insertion, but we can handle mutations before and
    # after it
    mp = [3, 10, 12, 22]
    mp2ra = {3: ("T", "G"), 10: ("G", "A"), 12: ("G", "A"), 22: ("C", "G")}
    repls = su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert repls == {3: "T", 10: "A", 12: "A", 22: "C"}


def test_get_smooth_aln_replacements_mp_and_mp2ra_differ_subtly():
    aln = fetch_specific_aln("c3", "TTTTTTTTTTTTTTT")
    mp = [6, 8, 11]
    mp2ra = {6: ("A", "T"), 8: ("G", "T")}
    with pytest.raises(WeirdError) as ei:
        su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert str(ei.value) == (
        "The mutated positions given in mutated_positions and mp2ra differ.\n"
        "mutated_positions: [6, 8, 11];\n"
        "mp2ra: [6, 8]"
    )


def test_get_smooth_aln_replacements_mp_and_mp2ra_differ_a_lot():
    aln = fetch_specific_aln("c3", "TTTTTTTTTTTTTTT")
    mp = [1, 2, 3]
    mp2ra = {6: ("A", "T"), 8: ("G", "T"), 9: ("T", "C")}
    with pytest.raises(WeirdError) as ei:
        su.get_smooth_aln_replacements(aln, mp, mp2ra)
    assert str(ei.value) == (
        "The mutated positions given in mutated_positions and mp2ra differ.\n"
        "mutated_positions: [1, 2, 3];\n"
        "mp2ra: [6, 8, 9]"
    )


def test_compute_average_coverages_verbose(capsys):
    bf = pysam.AlignmentFile(BAM)
    assert su.compute_average_coverages(
        FASTA, {"c1": 23}, bf, True, mock_log
    ) == {"c1": 12.0}
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Computing average coverages in each contig, for "
        "use with virtual reads...\n"
        "MockLog: On contig c1 (23 bp) (1 / 1 = 100.00% done).\n"
        "MockLog: c1 has average coverage 12.00x.\n"
        "MockLog: Done.\n"
    )


def test_compute_average_coverages_noverbose(capsys):
    bf = pysam.AlignmentFile(BAM)
    assert su.compute_average_coverages(
        FASTA, {"c1": 23}, bf, False, mock_log
    ) == {"c1": 12.0}
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Computing average coverages in each contig, for "
        "use with virtual reads...\n"
        "MockLog: Done.\n"
    )


def test_compute_average_coverages_badlength():
    bf = pysam.AlignmentFile(BAM)

    # Case 1: BAM has more than FASTA
    with pytest.raises(WeirdError) as ei:
        su.compute_average_coverages(FASTA, {"c1": 22}, bf, False, mock_log)
    assert str(ei.value) == (
        "Contig c1 has length 22 bp, but we saw 23 positions in the alignment "
        "for this contig."
    )

    # Case 2: BAM has less than FASTA
    with pytest.raises(WeirdError) as ei:
        su.compute_average_coverages(FASTA, {"c1": 12345}, bf, False, mock_log)
    assert str(ei.value) == (
        "Contig c1 has length 12,345 bp, but we saw 23 positions in the "
        "alignment for this contig."
    )
