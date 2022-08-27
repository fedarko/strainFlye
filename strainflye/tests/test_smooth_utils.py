import os
import stat
import gzip
import tempfile
import pytest
import pysam
import skbio
import strainflye.smooth_utils as su
from io import StringIO
from strainflye.config import DI_PREF
from strainflye.errors import ParameterError, WeirdError
from strainflye.tests.utils_for_testing import mock_log, mock_log_2


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


def test_get_average_coverages_from_di_single_contig():
    tsv = StringIO(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\n"
        "edge_1\t35.2\t100\t0.5\n"
    )
    assert su.get_average_coverages_from_di({"edge_1": 100}, tsv) == {
        "edge_1": 35.2
    }


def test_get_average_coverages_from_di_multiple_contigs():
    tsv_text = (
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\n"
        "edge_1\t35.2\t100\t0.5\n"
        "edge_2\t12.3\t456\t7.8\n"
        "edge_3\t77.7\t7777\t7.7\n"
    )

    # Case 1: sets of contig names match exactly
    tsv = StringIO(tsv_text)
    assert su.get_average_coverages_from_di(
        {"edge_1": 100, "edge_2": 456, "edge_3": 7777}, tsv
    ) == {"edge_1": 35.2, "edge_2": 12.3, "edge_3": 77.7}

    # Case 2: name2len is a subset of div idx's contigs
    tsv = StringIO(tsv_text)
    assert su.get_average_coverages_from_di(
        {"edge_2": 456, "edge_3": 7777}, tsv
    ) == {"edge_2": 12.3, "edge_3": 77.7}


def test_get_average_coverages_from_di_contig_not_in_di():
    tsv = StringIO(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\n"
        "edge_1\t35.2\t100\t0.5\n"
    )
    with pytest.raises(ParameterError) as ei:
        su.get_average_coverages_from_di({"edge_2": 100}, tsv)
    assert str(ei.value) == (
        "Can't find contig edge_2 in the diversity index file."
    )


def test_get_average_coverages_from_di_differing_contig_lengths():
    tsv = StringIO(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\n"
        "edge_1\t35.2\t100\t0.5\n"
    )
    with pytest.raises(ParameterError) as ei:
        su.get_average_coverages_from_di({"edge_1": 99}, tsv)
    assert str(ei.value) == (
        "Length of contig edge_1 is 99 bp according to the FASTA file, "
        "but the diversity index file says its length is 100 bp."
    )


def test_get_average_coverages_from_di_everything_is_on_fire():
    tsv = StringIO(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\n"
        "edge_1\t35.2\t100\t0.5\n"
    )
    # If both the "missing contig" and "length mismatch" problems happen, then
    # ... well they both can't happen at the same time, at least not for a
    # single contig! I guess one could happen after the other for different
    # contigs, so which error pops up first is dependent on how we traverse
    # through contig_name2len. Not worth testing, probably.
    #
    # But, look, in any case the "missing contig" error is what the user sees
    # if contig_name2len and the diversity index file are completely different,
    # as they are here.
    with pytest.raises(ParameterError) as ei:
        su.get_average_coverages_from_di({"edge_5": 99}, tsv)
    assert str(ei.value) == (
        "Can't find contig edge_5 in the diversity index file."
    )


def test_append_reads():
    # "programming is my passion"
    a_100_times = "A" * 100
    exp_initial_lines = [
        ">r1\n",
        "ACGT\n",
        ">r2\n",
        a_100_times + "\n",
        ">r3\n",
        "CGTAC\n",
    ]
    with tempfile.NamedTemporaryFile() as fh:
        su.append_reads(
            fh.name, {"r1": "ACGT", "r2": a_100_times, "r3": "CGTAC"}
        )
        with gzip.open(fh.name, "rt") as written_fh:
            assert written_fh.readlines() == exp_initial_lines

        # importantly, test that append_reads() *appends* to the end of the
        # file rather than overwriting the stuff that may have already been
        # in it
        su.append_reads(fh.name, {"surprise_its_more_reads": "TACAT"})

        with gzip.open(fh.name, "rt") as written_fh:
            assert written_fh.readlines() == exp_initial_lines + [
                ">surprise_its_more_reads\n",
                "TACAT\n",
            ]


def test_write_virtual_reads_length_disagreement():
    with tempfile.NamedTemporaryFile() as fh:
        with pytest.raises(WeirdError) as ei:
            su.write_virtual_reads(
                "c1",
                skbio.DNA("ACCGT"),
                100,
                [101, 99, 100, 100],
                1,
                0.5,
                fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "len(pos2srcov) == 4 bp, but len(contig_seq) == 5 bp."
        )


def verify_gz_file_empty(fp):
    with gzip.open(fp, "rt") as fh:
        assert fh.readlines() == []


def test_write_virtual_reads_no_low_coverage_positions(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        su.write_virtual_reads(
            "c1",
            skbio.DNA("ACGT"),
            100,
            [101, 99, 100, 100],
            1,
            0.5,
            fh.name,
            mock_log,
        )
        # check expected logging output
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 (average coverage 100.00x, based on the BAM "
            "file) has no low-coverage (\u2264 50.00x) positions (based on "
            "smoothed read coverages). No need to create virtual reads.\n"
        )
        # Make sure that no reads were written out -- the file should be empty,
        # or at least unmodified (by strainFlye) from before
        # write_virtual_reads() was called
        verify_gz_file_empty(fh.name)


def test_write_virtual_reads_one_low_coverage_position(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        # Add some fake data to the file before calling write_virtual_reads()
        # -- this lets us incidentally test that this data is not overwritten,
        # and that the virtual reads occur after this
        with gzip.open(fh.name, "wt") as pre_vr_fh:
            pre_vr_fh.write(">example_fake_smoothed_read\nCATCATCAT\n")

        su.write_virtual_reads(
            "c1",
            skbio.DNA("ACGT"),
            100,
            [101, 30, 100, 100],
            1,
            0.5,
            fh.name,
            mock_log,
        )

        # notably, the average coverage does NOT have to be equal to the
        # average of pos2srcov (and this will probably usually not happen in
        # practice). this is because, as i have tried to make clear in the logs
        # here, the average coverage is based on the BAM, while pos2srcov is
        # just based on coverages by smoothed reads.
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 (average coverage 100.00x, based on the BAM "
            "file) has 1 run(s) of consecutive low-coverage (\u2264 50.00x) "
            "positions (based on smoothed read coverages). Creating virtual "
            "reads...\n"
            "MockLog: Created 70 virtual read(s) total for contig c1.\n"
        )
        with gzip.open(fh.name, "rt") as written_fh:
            curr_readnum = 1
            for linenum, line in enumerate(written_fh):

                # The first two lines are a special case, since they represent
                # the smoothed read that was already there. Shouldn't have been
                # overwritten!
                if linenum == 0:
                    assert line == ">example_fake_smoothed_read\n"
                    continue
                elif linenum == 1:
                    assert line == "CATCATCAT\n"
                    continue

                # Okay, after the first two lines, we can care about the
                # virtual reads
                if linenum % 2 == 0:
                    assert line.startswith(">vr_1_1_")
                    split_line = line.strip().split("_")
                    assert len(split_line) == 4
                    vr_num = split_line[3]
                    assert int(vr_num) == curr_readnum
                else:
                    # we used a vr flank of 1, so we include 1 position before
                    # and after the "run" of this one low-coverage position (C)
                    assert line == "ACG\n"
                    curr_readnum += 1
            # we should've seen exactly 142 lines -- two lines for each of the
            # 70 virtual reads we should've created, plus two for that one
            # smoothed read. And linenum is 0-indexed, since it's managed by
            # enumerate() above.
            assert linenum == 141


def test_write_virtual_reads_multi_run_clamp_and_overlap(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        # two runs: one of (0-indexed) positions [2, 3], another of [6, 7]
        # notably, we use virtual_read_flank = 2 -- so the right run will have
        # to be clamped on the right side. also, there will be some overlap due
        # to these runs being so close together. it is what it is -- will note
        # in the paper
        su.write_virtual_reads(
            "c1",
            skbio.DNA("ACGTCCAC"),
            100,
            [101, 60, 30, 0, 50, 70, 20, 5],
            2,
            0.5,
            fh.name,
            mock_log,
        )
        # average coverage of left  run: (30 + 0) / 2 = 15.00x.
        # average coverage of right run: (20 + 5) / 2 = 12.50x.
        # So, gotta create 85 VRs for left, and
        # round(100 - 12.5) == round(87.5) == 88 VRs for right.
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 (average coverage 100.00x, based on the BAM "
            "file) has 2 run(s) of consecutive low-coverage (\u2264 50.00x) "
            "positions (based on smoothed read coverages). Creating virtual "
            "reads...\n"
            "MockLog: Created 173 virtual read(s) total for contig c1.\n"
        )
        with gzip.open(fh.name, "rt") as written_fh:
            # Go through and verify that the reads look as we expect. The
            # leftmost run's reads should be included first, although please
            # note that this is less a strict guarantee of strainFlye and more
            # an artifact of how convert_to_runs() works.
            #
            # (note that linenum here is not bundled in with enumerate(), we
            # manage it ourselves. sorry this test code is kinda gross.)
            linenum = 0
            curr_readnum = 1
            for line in written_fh:
                # line 170 (0-indexed) is the "cross-over point" -- after this,
                # we get to the reads from the rightmost run
                if linenum == 170:
                    curr_readnum = 1
                if linenum <= 169:
                    if linenum % 2 == 0:
                        assert line.startswith(">vr_2_3_")
                        split_line = line.strip().split("_")
                        assert len(split_line) == 4
                        vr_num = split_line[3]
                        assert int(vr_num) == curr_readnum
                    else:
                        assert line == "ACGTCC\n"
                        curr_readnum += 1
                    linenum += 1
                else:
                    if linenum % 2 == 0:
                        assert line.startswith(">vr_6_7_")
                        split_line = line.strip().split("_")
                        assert len(split_line) == 4
                        vr_num = split_line[3]
                        assert int(vr_num) == curr_readnum
                    else:
                        # should have been clamped on the right side
                        assert line == "CCAC\n"
                        curr_readnum += 1
                    linenum += 1
            assert linenum == 346


def test_write_smoothed_reads_zero_mutations():
    with tempfile.NamedTemporaryFile() as fh:
        with pytest.raises(WeirdError) as ei:
            su.write_smoothed_reads(
                "c1",
                FASTA,
                23,
                {},
                pysam.AlignmentFile(BAM),
                True,
                fh.name,
                mock_log_2,
                mock_log,
            )
        assert str(ei.value) == (
            "Contig c1 has zero mutations. Can't create smoothed reads."
        )


def test_write_smoothed_reads_bad_chunk_size():
    # We don't test for the case where scs is a float, because we don't
    # explicitly check for that in write_smoothed_reads(). I guess we could,
    # but I've made the assumption that things will be of the correct type
    # usually.
    for scs in (-100, -2, -1, 0):
        with tempfile.NamedTemporaryFile() as fh:
            with pytest.raises(WeirdError) as ei:
                su.write_smoothed_reads(
                    "c1",
                    FASTA,
                    23,
                    {10: ("G", "A")},
                    pysam.AlignmentFile(BAM),
                    True,
                    fh.name,
                    mock_log_2,
                    mock_log,
                    sr_chunk_size=scs,
                )
            assert str(ei.value) == (
                f"sr_chunk_size should be a positive integer, but it's {scs}?"
            )


def verify_c1_1mut_smoothedreads(reads_fp):
    with gzip.open(reads_fp, "rt") as written_fh:
        for linenum, line in enumerate(written_fh):
            if linenum % 2 == 0:
                exp_read_num = int((linenum / 2)) + 1
                # since read names go r01, r02, ... r09, r10, r11, ...
                # something something left pad joke goes here
                exp_read_name = "r" + str(exp_read_num).zfill(2)
                # The _1 is because each of these reads is only aligned
                # once to c1. If we'd have supp alignments, then we'd see
                # _2, etc.
                assert line == f">{exp_read_name}_1\n"
            else:
                # Each smoothed read should spell out the exact c1
                # sequence, except for the lone mutation at pos 10. (For
                # the first 5 alignments, there's an A here; for the
                # remaining alignments, there's a G.)
                if linenum <= 9:
                    assert line == "ACTGACACCCAAACCAAACCTAC\n"
                else:
                    assert line == "ACTGACACCCGAACCAAACCTAC\n"
        # Should've seen exactly 24 lines (12 smoothed reads including
        # header + sequence). linenum is zero-indexed.
        assert linenum == 23


def test_write_smoothed_reads_c1_one_mutation_basic(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        # Let's act like c1 only has one mutation, at position 10 (0-indexed).
        # Five of the reads aligned to c1 have an A at pos 10, and seven of the
        # reads aligned to c1 have a G at pos 10. So we should expect to see
        # some very simple smoothed reads reflecting this.
        pos2srcov, contig_seq = su.write_smoothed_reads(
            "c1",
            FASTA,
            23,
            {10: ("G", "A")},
            pysam.AlignmentFile(BAM),
            True,
            fh.name,
            mock_log_2,
            mock_log,
        )
        # First, let's test return values.
        #
        # Each of the 12 reads aligned to c1 can be a smoothed read -- and
        # since these reads all span all of c1, the coverage in c1 should be
        # uniformly set to 12x.
        assert pos2srcov == [12] * 23
        # And we should be able to load c1's sequence without problems.
        assert str(contig_seq) == "ACTGACACCCAAACCAAACCTAC"

        # Second, let's test that the log output looks good.
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 has 1 mutated position(s).\n"
            "MockLog: From the 12 linear alignment(s) to contig c1: created "
            "12 smoothed read(s) and ignored 0 linear alignment(s).\n"
        )

        # Finally, let's verify that the smoothed reads that were written out
        # are correct.
        verify_c1_1mut_smoothedreads(fh.name)


def test_write_smoothed_reads_c1_one_mutation_tiny_buffer(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        pos2srcov, contig_seq = su.write_smoothed_reads(
            "c1",
            FASTA,
            23,
            {10: ("G", "A")},
            pysam.AlignmentFile(BAM),
            True,
            fh.name,
            mock_log_2,
            mock_log,
            # Setting sr_chunk_size to 1 will make us do a lot of write
            # operations, but we should get the same results
            sr_chunk_size=1,
        )
        assert pos2srcov == [12] * 23
        assert str(contig_seq) == "ACTGACACCCAAACCAAACCTAC"
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 has 1 mutated position(s).\n"
            "MockLog: From the 12 linear alignment(s) to contig c1: created "
            "12 smoothed read(s) and ignored 0 linear alignment(s).\n"
        )
        verify_c1_1mut_smoothedreads(fh.name)


def test_write_smoothed_reads_c1_one_mutation_no_pos2srcov(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        pos2srcov, contig_seq = su.write_smoothed_reads(
            "c1",
            FASTA,
            23,
            {10: ("G", "A")},
            pysam.AlignmentFile(BAM),
            False,
            fh.name,
            mock_log_2,
            mock_log,
        )
        assert pos2srcov is None
        assert str(contig_seq) == "ACTGACACCCAAACCAAACCTAC"
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 has 1 mutated position(s).\n"
            "MockLog: From the 12 linear alignment(s) to contig c1: created "
            "12 smoothed read(s) and ignored 0 linear alignment(s).\n"
        )
        verify_c1_1mut_smoothedreads(fh.name)


def test_write_smoothed_reads_bad_length():
    with tempfile.NamedTemporaryFile() as fh:
        with pytest.raises(WeirdError) as ei:
            pos2srcov, contig_seq = su.write_smoothed_reads(
                "c1",
                FASTA,
                24,
                {10: ("G", "A")},
                pysam.AlignmentFile(BAM),
                False,
                fh.name,
                mock_log_2,
                mock_log,
            )
        assert str(ei.value) == (
            "len(contig_seq) == 23 bp, but contig_len == 24 bp."
        )


def test_write_smoothed_reads_c1_one_mutation_all_alns_ignored(capsys):
    with tempfile.NamedTemporaryFile() as fh:
        pos2srcov, contig_seq = su.write_smoothed_reads(
            "c1",
            FASTA,
            23,
            # All of c1's alignments to 10 have an A or a G. If we act like
            # this mutation doesn't include either of these nts, then we'll
            # ignore all of these alignments.
            {10: ("C", "T")},
            pysam.AlignmentFile(BAM),
            False,
            fh.name,
            mock_log_2,
            mock_log,
        )
        assert pos2srcov is None
        assert contig_seq is None
        # MockLog2 (aka the non-verbose-only logger) should be used for this
        # warning, because this is super weird and should hopefully be uncommon
        assert capsys.readouterr().out == (
            "MockLog: Contig c1 has 1 mutated position(s).\n"
            "MockLog2: Ignored all 12 linear alignments to contig c1: couldn't "
            "generate any smoothed reads. Ignoring this contig.\n"
        )
        # Nothing shoulda gotten written out
        verify_gz_file_empty(fh.name)
