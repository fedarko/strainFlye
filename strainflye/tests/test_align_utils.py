import os
import subprocess
import pysam
import pytest
import strainflye.align_utils as au
from strainflye.tests.utils_for_testing import mock_log
from strainflye.errors import SequencingDataError, GraphParsingError

TI_DIR = os.path.join("strainflye", "tests", "inputs")
S1 = os.path.join(TI_DIR, "sample1.gfa")


def test_get_coords_good():
    samfile = pysam.AlignmentFile(os.path.join(TI_DIR, "camp-sf.sam"), "r")
    for aln in samfile.fetch():
        assert au.get_coords(aln) == (1162486, 1259415)


def test_check_contigs_in_graph_good():
    au.check_contigs_in_graph(
        {"1": 8, "2": 10, "3": 21, "4": 7, "5": 8, "6": 4}, S1
    )


def test_check_contigs_in_graph_missing():
    bad_inputs = [
        # Contig "2" is in the graph but not the fasta
        {"1": 8, "3": 21, "4": 7, "5": 8, "6": 4},
        # Contig "99" is in the fasta but not the graph
        {"1": 8, "2": 10, "3": 21, "4": 7, "5": 8, "6": 4, "99": 100},
        {},
        {"1": 8, "99": 100},
    ]

    for n2l in bad_inputs:
        with pytest.raises(SequencingDataError) as ei:
            au.check_contigs_in_graph(n2l, S1)
        assert str(ei.value) == (
            "Segment names in the GFA file don't match contig names in the "
            "FASTA file."
        )


def test_check_contigs_in_graph_empty_gfa():
    # Fortunately, load_graph() already throws an error if the GFA file is
    # empty -- so this will result in a different error and message
    gf = os.path.join(TI_DIR, "sample1-empty.gfa")
    with pytest.raises(GraphParsingError) as ei:
        au.check_contigs_in_graph({"1": 8}, gf)
    assert str(ei.value) == f"Less than 1 segment(s) are given in {gf}."


def test_check_contigs_in_graph_diff_lengths():
    with pytest.raises(SequencingDataError) as ei:
        au.check_contigs_in_graph(
            {"1": 9, "2": 10, "3": 21, "4": 7, "5": 8, "6": 4}, S1
        )
    assert str(ei.value) == (
        "Contig 1 has length 9 in the FASTA file, but length 8 in the GFA "
        "file."
    )

    with pytest.raises(SequencingDataError) as ei:
        au.check_contigs_in_graph(
            {"1": 8, "2": 10, "3": 21, "4": 7, "5": 8, "6": 4000}, S1
        )
    assert str(ei.value) == (
        "Contig 6 has length 4,000 in the FASTA file, but length 4 in the GFA "
        "file."
    )


def test_filter_osa_reads_basic(capsys, tmp_path):
    sam_fp = os.path.join(tmp_path, "aln.sam")
    with open(sam_fp, "w") as fh:
        # r02 should be filtered out, since it overlaps with itself at position
        # 5 on the reference. everything else stays. (r12 is the sort of thing
        # the partially-mapped read filter should remove, and that'd happen
        # after the OSA filter step.)
        #
        # Diagram on c1, for reference:
        #
        #                  11111111112222
        #         12345678901234567890123
        # contig: ACTGACACCCAAACCAAACCTAC
        #    r01: ACTGACACCCAAACCAAACCTAC
        #    r02: ACTGA
        #    r02:     ACACCC
        #    r04: ACT
        #    r04:                   CCTAC
        #    r12: TAAAAAGGGGGG
        #
        # This illustrates that r02 is the only read with an OSA.
        fh.write(
            "@HD	VN:1.6	SO:coordinate\n"
            "@SQ	SN:c1	LN:23\n"
            "@SQ	SN:c2	LN:12\n"
            "@SQ	SN:c4	LN:100\n"
            "r01	0	c1	1	30	23M	*	0	0	ACTGACACCCAAACCAAACCTAC	*\n"
            "r02	0	c1	1	30	5M5S	*	0	0	ACTGACACCC	*\n"
            "r04	2048	c1	1	30	3M5S	*	0	0	ACTCCTAC	*\n"
            "r12	0	c1	1	30	12M	*	0	0	TAAAAAGGGGGG	*\n"
            "r02	2048	c1	5	30	4S6M	*	0	0	ACTGACACCC	*\n"
            "r04	0	c1	19	30	3S5M	*	0	0	ACTCCTAC	*\n"
            "r12	2048	c2	1	30	12M	*	0	0	TAAAAAGGGGGG	*\n"
            "r13	0	c2	1	30	12M	*	0	0	TAAAAAGGGGGG	*\n"
            "r14	0	c2	1	30	12M	*	0	0	TAAAAAGGGGGG	*\n"
        )
    in_bam_fp = os.path.join(tmp_path, "aln.bam")
    subprocess.run(
        ["samtools", "view", "-b", sam_fp, "-o", in_bam_fp], check=True
    )
    # i mean i guess let's incidentally test this also
    au.index_bam(in_bam_fp, "test BAM", mock_log)
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Indexing the test BAM...\n"
        "MockLog: Done indexing the test BAM.\n"
    )

    out_bam_fp = os.path.join(tmp_path, "osa-filtered.bam")

    # okay now run the OSA filter
    au.filter_osa_reads(in_bam_fp, out_bam_fp, mock_log, True)

    # verify the logged output looks good (we do this first, before inspecting
    # the actual output BAM file, because the indexing functions also log
    # output; it's simplest to just do stuff in this order).
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Filtering reads with overlapping supplementary "
        "alignments (OSAs)...\n"
        "MockLog: "
        "OSA filter pass 1/2: on contig c1 (1 / 3 contigs = 33.33%).\n"
        "MockLog: There are 6 linear alignment(s) (from 4 unique read(s)) to "
        "contig c1.\n"
        "MockLog: 1 / 4 (25.00%) of these unique read(s) have OSAs.\n"
        "MockLog: "
        "OSA filter pass 1/2: on contig c2 (2 / 3 contigs = 66.67%).\n"
        "MockLog: There are 3 linear alignment(s) (from 3 unique read(s)) to "
        "contig c2.\n"
        "MockLog: 0 / 3 (0.00%) of these unique read(s) have OSAs.\n"
        "MockLog: "
        "OSA filter pass 1/2: on contig c4 (3 / 3 contigs = 100.00%).\n"
        "MockLog: Nothing is aligned to contig c4! Ignoring this contig.\n"
        "MockLog: Done with pass 1 of the OSA filter; moving on to pass 2...\n"
        "MockLog: "
        "OSA filter pass 2/2: on contig c1 (1 / 3 contigs = 33.33%).\n"
        "MockLog: 4 / 6 (66.67%) linear aln(s) retained in contig c1.\n"
        "MockLog: "
        "OSA filter pass 2/2: on contig c2 (2 / 3 contigs = 66.67%).\n"
        "MockLog: 3 / 3 (100.00%) linear aln(s) retained in contig c2.\n"
        "MockLog: Done filtering reads with overlapping supplementary "
        "alignments.\n"
    )

    # Inspect the OSA-filtered BAM -- verify it's actually correct!

    # (we have to index it to load it in pysam, tho)
    au.index_bam(out_bam_fp, "OSA-filtered test BAM", mock_log)
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Indexing the OSA-filtered test BAM...\n"
        "MockLog: Done indexing the OSA-filtered test BAM.\n"
    )

    bf = pysam.AlignmentFile(out_bam_fp, "rb")

    # check c1 first
    # (I'm pretty sure this order of iteration through the linear alignments is
    # guaranteed due to the order staying the same from SAM --> BAM, but don't
    # quote me on that -- sorry if this test starts breaking 20 years from
    # now...)
    exp_read_names = ["r01", "r04", "r12", "r04"]
    exp_qaln_seqs = ["ACTGACACCCAAACCAAACCTAC", "ACT", "TAAAAAGGGGGG", "CCTAC"]
    obs_read_names = []
    obs_qaln_seqs = []

    for linearaln in bf.fetch("c1"):
        obs_read_names.append(linearaln.query_name)
        # sanity check that pysam understands CIGAR strings!
        obs_qaln_seqs.append(linearaln.query_alignment_sequence)

    assert obs_read_names == exp_read_names
    assert obs_qaln_seqs == exp_qaln_seqs

    # check c2
    exp_read_names = ["r12", "r13", "r14"]
    exp_qaln_seqs = ["TAAAAAGGGGGG"] * 3
    obs_read_names = []
    obs_qaln_seqs = []
    for linearaln in bf.fetch("c2"):
        obs_read_names.append(linearaln.query_name)
        obs_qaln_seqs.append(linearaln.query_alignment_sequence)

    assert obs_read_names == exp_read_names
    assert obs_qaln_seqs == exp_qaln_seqs

    # check c4
    obs_read_names = []
    for linearaln in bf.fetch("c4"):
        obs_read_names.append(linearaln.query_name)
    assert obs_read_names == []
