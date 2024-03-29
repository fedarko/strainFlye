import os
import copy
import pickle
import pysam
import networkx as nx
import pytest
import strainflye.link_utils as lu
from collections import defaultdict
from pytest import approx
from strainflye.config import (
    POS_FILE_LBL,
    POSPAIR_FILE_LBL,
    LG_DOT_HEADER_ATTRS,
)
from strainflye.errors import ParameterError, WeirdError
from strainflye.tests.utils_for_testing import mock_log


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BCF = os.path.join(IN_DIR, "call-r-min3-di12345", "naive-calls.bcf")
BAM = os.path.join(IN_DIR, "alignment.bam")
DEGEN_BAM = os.path.join(IN_DIR, "degen.bam")
SECONDARY_BAM = os.path.join(IN_DIR, "c4-and-secondary.bam")
NO_C1_READS_BAM = os.path.join(IN_DIR, "no-c1-reads.bam")


def test_gen_ddi():
    assert lu.gen_ddi() == defaultdict(int)


def test_get_readname2pos2nt_good():
    bf = pysam.AlignmentFile(BAM, "rb")

    # Note that A/C/G/T are encoded in the inner dicts as 0/1/2/3
    final_few_read_haplotypes = {0: 0, 3: 2, 22: 1}

    assert lu.get_readname2pos2nt(bf, "c1", [0, 3, 22]) == {
        "r01": {0: 0, 3: 2, 22: 1},
        "r02": {0: 0, 3: 3, 22: 1},
        "r03": {0: 0, 3: 1, 22: 1},
        "r04": {0: 0, 3: 0, 22: 1},
        "r05": {0: 0, 3: 1, 22: 1},
        "r06": {0: 0, 3: 3, 22: 1},
        "r07": {0: 0, 3: 3, 22: 1},
        "r08": final_few_read_haplotypes,
        "r09": final_few_read_haplotypes,
        "r10": final_few_read_haplotypes,
        "r11": final_few_read_haplotypes,
        "r12": final_few_read_haplotypes,
    }


def test_get_readname2pos2nt_onepos_with_nonspanning_read():
    bf = pysam.AlignmentFile(BAM, "rb")
    r2p2t = lu.get_readname2pos2nt(bf, "c2", [0])

    assert r2p2t == {
        "r13": {0: 3},
        "r14": {0: 3},
        "r15": {0: 3},
        "r16": {0: 3},
        "r17": {0: 3},
        "r18": {0: 3},
        "r19": {0: 3},
        "r20": {0: 3},
        "r21": {0: 3},
        "r22": {0: 3},
    }
    # read r23 (which skips over the first position in c2) shouldn't be
    # included in r2p2t; however, r2p2t is a defaultdict, so we should be able
    # to try to access r23 without any errors. (this behavior makes sense
    # because it implies that this read simply just doesn't span any positions
    # of interest -- by avoiding storing these currently-"unhelpful" reads in
    # r2p2t, we save some space.)
    assert r2p2t["r23"] == {}


def test_get_readname2pos2nt_degen_doesnt_span():
    bf = pysam.AlignmentFile(DEGEN_BAM, "rb")
    # All reads have a "C" aligned to zero-indexed position 8, aside from reads
    # r06, r07, and r08 (which in "degen.bam" have a "N" aligned to this
    # position). So, these three reads should not be counted as spanning this
    # position, even if they do span the position right next to it (9).
    assert lu.get_readname2pos2nt(bf, "c1", [8, 9]) == {
        "r01": {8: 1, 9: 1},
        "r02": {8: 1, 9: 1},
        "r03": {8: 1, 9: 1},
        "r04": {8: 1, 9: 1},
        "r05": {8: 1, 9: 1},
        "r06": {9: 1},
        "r07": {9: 1},
        "r08": {9: 1},
        "r09": {8: 1, 9: 1},
        "r10": {8: 1, 9: 1},
        "r11": {8: 1, 9: 1},
        "r12": {8: 1, 9: 1},
    }


def test_get_readname2pos2nt_err_on_read_self_overlap():
    bf = pysam.AlignmentFile(SECONDARY_BAM, "rb")
    with pytest.raises(ParameterError) as ei:
        lu.get_readname2pos2nt(bf, "c1", [8, 9])
    assert str(ei.value) == (
        "Read r11 is aligned to the 0-indexed position 8 in contig c1 "
        "multiple times. Make sure that you have removed secondary "
        "alignments and overlapping supplementary alignments from your "
        'alignment file; you can use "strainFlye align" to compute '
        "such an alignment."
    )


def test_get_pos_nt_info_just_one_pos():
    # (We already have a simple "good" test of this function as a doctest in
    # link_utils, just to make the docstring there easier to read by itself.)
    pos2nt2ct, pospair2ntpair2ct = lu.get_pos_nt_info(
        {"r1": {0: 3}, "r2": {0: 1}, "r3": {0: 1}}
    )
    # Note that positions get converted to 1-indexing
    assert pos2nt2ct == {
        1: {1: 2, 3: 1},
    }
    assert pospair2ntpair2ct == {}


def check_c3_nt_info(tmp_path):
    fp_c3p = tmp_path / f"c3_{POS_FILE_LBL}.pickle"
    fp_c3pp = tmp_path / f"c3_{POSPAIR_FILE_LBL}.pickle"

    with open(fp_c3p, "rb") as f:
        c3p = pickle.load(f)
        assert c3p == {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}}

    with open(fp_c3pp, "rb") as f:
        c3pp = pickle.load(f)
        assert c3pp == {
            (7, 8): {
                # (A, T)
                (0, 3): 6,
                # (A, C)
                (0, 1): 1,
                # (T, C)
                (3, 1): 2,
                # (T, T)
                (3, 3): 4,
            }
        }


def test_run_nt_good(capsys, tmp_path):
    lu.run_nt(FASTA, BAM, BCF, tmp_path, True, mock_log)

    fp_c1p = tmp_path / f"c1_{POS_FILE_LBL}.pickle"
    fp_c1pp = tmp_path / f"c1_{POSPAIR_FILE_LBL}.pickle"

    with open(fp_c1p, "rb") as f:
        c1p = pickle.load(f)
        assert c1p == {
            4: {0: 1, 1: 2, 2: 6, 3: 3},
            11: {0: 5, 2: 7},
            13: {0: 5, 2: 7},
        }

    with open(fp_c1pp, "rb") as f:
        c1pp = pickle.load(f)
        # since c1's mutations at positions 11 and 13 are the same across all
        # reads, the pair data for (4, 11) and (4, 13) is identical
        first_two_pair_data = {
            # (A, A)
            (0, 0): 1,
            # (C, A)
            (1, 0): 2,
            # (G, A)
            (2, 0): 1,
            # (G, G)
            (2, 2): 5,
            # (T, A)
            (3, 0): 1,
            # (T, G)
            (3, 2): 2,
        }
        assert c1pp == {
            (4, 11): first_two_pair_data,
            (4, 13): first_two_pair_data,
            # this one is simple, because these mutations are the same across
            # all reads (ofc it's worth noting this is a simple tiny example
            # used for testing and real data will be more complex)
            (11, 13): {
                (0, 0): 5,
                (2, 2): 7,
            },
        }

    check_c3_nt_info(tmp_path)

    exp_out = (
        (
            "PREFIX\nMockLog: Loading and checking FASTA, BAM, and BCF "
            "files...\n"
            "MockLog: The FASTA file describes 3 contig(s).\n"
            "MockLog: All FASTA contig(s) are included in the BAM file (this "
            "BAM file has 3 reference(s)).\n"
            "MockLog: All FASTA contig(s) are included in the BCF file (the "
            "header of this BCF file describes 3 contig(s)).\n"
            "MockLog: The lengths of all contig(s) in the FASTA file match "
            "the corresponding lengths in the BAM and BCF files.\n"
            "MockLog: So far, these files seem good.\n"
        )
        + (
            "PREFIX\nMockLog: Going through contigs and computing nucleotide "
            "(co-)occurrence information...\n"
        )
        + (
            "MockLog: On contig c1 (23 bp) (1 / 3 contigs = 33.33%).\n"
            "MockLog: Contig c1 has 3 mutated position(s). Going through "
            "alignments to this contig...\n"
            "MockLog: Found 12 read(s) spanning \u2265 1 mutated position in "
            "contig c1. Computing (co-)occurrence info...\n"
            "MockLog: Wrote out (co-)occurrence info for contig c1.\n"
        )
        + (
            "MockLog: On contig c2 (12 bp) (2 / 3 contigs = 66.67%).\n"
            "MockLog: Contig c2 has no mutations; ignoring this contig.\n"
        )
        + (
            "MockLog: On contig c3 (16 bp) (3 / 3 contigs = 100.00%).\n"
            "MockLog: Contig c3 has 2 mutated position(s). Going through "
            "alignments to this contig...\n"
            "MockLog: Found 13 read(s) spanning \u2265 1 mutated position in "
            "contig c3. Computing (co-)occurrence info...\n"
            "MockLog: Wrote out (co-)occurrence info for contig c3.\n"
        )
        + ("MockLog: Done.\n")
    )

    assert capsys.readouterr().out == exp_out


def test_run_nt_no_reads_spanning_muts(capsys, tmp_path):
    # use the same contigs and mutation calls as test_run_nt_good(), but remove
    # all linear alignments to contig c1. This causes all mutations in c1 to be
    # uncovered, so we shouldn't write out any info for it.
    lu.run_nt(FASTA, NO_C1_READS_BAM, BCF, tmp_path, False, mock_log)

    # We should *only* see c3 info in the output directory -- nothing for c1.
    # (And nothing for c2, for that matter, but that was already the case...)
    check_c3_nt_info(tmp_path)
    assert len(os.listdir(tmp_path)) == 2

    # we turned off "verbose", so we should see less output -- but we should
    # still see the warning about this weird case.
    exp_out = (
        (
            "PREFIX\nMockLog: Loading and checking FASTA, BAM, and BCF "
            "files...\n"
            "MockLog: The FASTA file describes 3 contig(s).\n"
            "MockLog: All FASTA contig(s) are included in the BAM file (this "
            "BAM file has 3 reference(s)).\n"
            "MockLog: All FASTA contig(s) are included in the BCF file (the "
            "header of this BCF file describes 3 contig(s)).\n"
            "MockLog: The lengths of all contig(s) in the FASTA file match "
            "the corresponding lengths in the BAM and BCF files.\n"
            "MockLog: So far, these files seem good.\n"
        )
        + (
            "PREFIX\nMockLog: Going through contigs and computing nucleotide "
            "(co-)occurrence information...\n"
        )
        + (
            "MockLog: Warning: Contig c1 has no reads spanning any of its "
            "mutated positions; ignoring this contig.\n"
        )
        + ("MockLog: Done.\n")
    )

    assert capsys.readouterr().out == exp_out


def test_make_linkgraph_good(capsys):
    g = lu.make_linkgraph(
        {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}},
        {
            (7, 8): {
                # (A, T)
                (0, 3): 6,
                # (A, C)
                (0, 1): 1,
                # (T, C)
                (3, 1): 2,
                # (T, T)
                (3, 3): 4,
            }
        },
        1,
        1,
        0,
        "c3",
        mock_log,
    )

    check_c3_linkgraph(g)

    assert capsys.readouterr().out == (
        "MockLog: Creating a link graph for contig c3...\n"
        "MockLog: The link graph for contig c3 has 4 node(s).\n"
        "MockLog: The link graph for contig c3 has 4 edge(s).\n"
    )


def test_make_linkgraph_uncovered_pos():
    with pytest.raises(WeirdError) as ei:
        lu.make_linkgraph(
            {7: {0: 0, 3: 0}, 8: {1: 3, 3: 10}},
            {
                (7, 8): {
                    # (A, T)
                    (0, 3): 6,
                }
            },
            1,
            1,
            0,
            "c3",
            mock_log,
        )
    assert str(ei.value) == "Coverage at pos 7 in contig c3 is 0?"


def test_make_linkgraph_lowcov_nt():
    # Let's set min_nt_ct to 4 -- this will ensure that the node (8, 1) will be
    # excluded from the graph
    g = lu.make_linkgraph(
        {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}},
        {
            (7, 8): {
                # (A, T)
                (0, 3): 6,
                # (A, C)
                (0, 1): 1,
                # (T, C)
                (3, 1): 2,
                # (T, T)
                (3, 3): 4,
            }
        },
        4,
        1,
        0,
        "c3",
        mock_log,
    )
    assert len(g.nodes) == 3
    assert set(g.nodes) == set([(7, 0), (7, 3), (8, 3)])
    assert g.nodes[(7, 0)] == {"ct": 7, "freq": 7 / 13}
    assert g.nodes[(7, 3)] == {"ct": 6, "freq": 6 / 13}
    # Importantly, the frequency of (8, 3) remains unchanged -- we should not
    # forget about the existence of (8, 1) entirely
    assert g.nodes[(8, 3)] == {"ct": 10, "freq": 10 / 13}

    assert len(g.edges) == 2
    assert g.edges[(7, 0), (8, 3)] == {"link": 6 / 10}
    assert g.edges[(7, 3), (8, 3)] == {"link": 4 / 10}


def test_make_linkgraph_non_spanned_pospair():
    with pytest.raises(WeirdError) as ei:
        # weird example where no reads span both 7 and 8. this is technically
        # possible if you have, like, the shittiest sequencing technology known
        # to humanity, i guess. can i patent that?
        lu.make_linkgraph(
            {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}},
            {
                (7, 8): {
                    (0, 3): 0,
                    (0, 1): 0,
                    (3, 1): 0,
                    (3, 3): 0,
                }
            },
            1,
            1,
            0,
            "c3",
            mock_log,
        )
    assert str(ei.value) == (
        "Number of spanning reads for position pair (7, 8) in contig c3 is 0?"
    )


def test_make_linkgraph_low_span_pospair():
    # Let's set min_span to 14, which will preclude any edges between 7 and 8
    # (in total, 13 reads span these positions)
    g = lu.make_linkgraph(
        {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}},
        {
            (7, 8): {
                # (A, T)
                (0, 3): 6,
                # (A, C)
                (0, 1): 1,
                # (T, C)
                (3, 1): 2,
                # (T, T)
                (3, 3): 4,
            }
        },
        1,
        14,
        0,
        "c3",
        mock_log,
    )
    assert len(g.nodes) == 4
    assert set(g.nodes) == set([(7, 0), (7, 3), (8, 1), (8, 3)])
    assert g.nodes[(7, 0)] == {"ct": 7, "freq": 7 / 13}
    assert g.nodes[(7, 3)] == {"ct": 6, "freq": 6 / 13}
    assert g.nodes[(8, 1)] == {"ct": 3, "freq": 3 / 13}
    assert g.nodes[(8, 3)] == {"ct": 10, "freq": 10 / 13}

    # the main thing we check for
    assert len(g.edges) == 0


def test_make_linkgraph_low_link_ntpairs():
    # Let's set low_link to 1 / 3. This will result in two edges:
    #
    # (7, 0) -- (8, 1) (link = 1 / 7)
    #
    # (7, 3) -- (8, 1) (link = 2 / 6 = 1 / 3, and the low_link check uses >,
    #                   although note that floating point stuff makes this
    #                   check imperfect so whatevs)
    #
    # ... not being included in the graph.
    g = lu.make_linkgraph(
        {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}},
        {
            (7, 8): {
                # (A, T)
                (0, 3): 6,
                # (A, C)
                (0, 1): 1,
                # (T, C)
                (3, 1): 2,
                # (T, T)
                (3, 3): 4,
            }
        },
        1,
        1,
        1 / 3,
        "c3",
        mock_log,
    )
    assert len(g.nodes) == 4
    assert set(g.nodes) == set([(7, 0), (7, 3), (8, 1), (8, 3)])
    assert g.nodes[(7, 0)] == {"ct": 7, "freq": 7 / 13}
    assert g.nodes[(7, 3)] == {"ct": 6, "freq": 6 / 13}
    assert g.nodes[(8, 1)] == {"ct": 3, "freq": 3 / 13}
    assert g.nodes[(8, 3)] == {"ct": 10, "freq": 10 / 13}

    assert len(g.edges) == 2
    assert g.edges[(7, 0), (8, 3)] == {"link": 6 / 10}
    assert g.edges[(7, 3), (8, 3)] == {"link": 4 / 10}


def test_make_linkgraph_zeroct_ntpair():
    with pytest.raises(WeirdError) as ei:
        # Set the count of (7, 3) -- (8, 3) to zero. In practice, we shouldn't
        # include non-co-occurring pairs of nucleotides in pospair2ntpair2ct.
        lu.make_linkgraph(
            {7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}},
            {
                (7, 8): {
                    (0, 3): 6,
                    (0, 1): 1,
                    (3, 1): 2,
                    (3, 3): 0,
                }
            },
            1,
            1,
            0,
            "c3",
            mock_log,
        )
    assert str(ei.value) == (
        "Numerator of link() for nodes (7, 3) and (8, 3) in contig c3 is 0?"
    )


def test_make_linkgraph_zero_den():
    with pytest.raises(WeirdError) as ei:
        # Set both (7, 3) and (8, 3) to have zero counts. These positions will
        # still have positive coverage, but when we try to create an edge
        # between these positions we'll see that the max count across either of
        # them is zero.
        #
        # Notably, to do this, we have to commit a sin -- we have to set
        # min_nt_ct to zero, which is straight-up impossible using strainFlye's
        # CLI at the moment. The reason for this is that otherwise (if we left
        # min_nt_ct at one) then these nodes with zero counts would get dropped
        # from the graph automatically.
        #
        # So, this test is impossible to replicate in practice (at least not
        # without doing some weird stuff to bypass Click's IntRange checks),
        # but whatevs. This test is just "this should never happen" on
        # steroids, I guess.
        lu.make_linkgraph(
            {7: {0: 7, 3: 0}, 8: {1: 3, 3: 0}},
            {
                (7, 8): {
                    (0, 3): 6,
                    (0, 1): 1,
                    (3, 1): 2,
                    (3, 3): 4,
                }
            },
            0,
            1,
            0,
            "c3",
            mock_log,
        )
    assert str(ei.value) == (
        "Denominator of link() for nodes (7, 3) and (8, 3) in contig c3 is 0?"
    )


def test_write_linkgraph_to_dot_good(tmp_path):
    g = nx.Graph()

    g.add_node((100, 1), ct=5, freq=0.25)
    g.add_node((100, 2), ct=15, freq=0.75)
    g.add_node((200, 0), ct=20, freq=1)

    g.add_edge((100, 1), (200, 0), link=0.25)
    g.add_edge((100, 2), (200, 0), link=0.75)

    lu.write_linkgraph_to_dot(g, tmp_path, "borgar", 0.01, 5)

    with open(tmp_path / "borgar_linkgraph.gv", "r") as f:
        dot_txt = f.read()

    # As of writing, this works -- however, if this test starts failing, this
    # may be due to NetworkX changing the order of its iteration through nodes
    # and/or edges. We could fix this in the actual code by sorting nodes/edges
    # before iteration -- or, we could just fix this in the test by testing
    # that the sets of strings representing each line exist in both files.
    # Something like that. (But I'm not gonna do this until it's necessary.)
    assert dot_txt == (
        "graph {\n"
        f"{LG_DOT_HEADER_ATTRS}"
        '  "(100, 1)" [label="100 (C)\\n5x (25.00%)"];\n'
        '  "(100, 2)" [label="100 (G)\\n15x (75.00%)"];\n'
        '  "(200, 0)" [label="200 (A)\\n20x (100.00%)"];\n'
        '  "(100, 1)" -- "(200, 0)" [penwidth=1.25];\n'
        '  "(100, 2)" -- "(200, 0)" [penwidth=3.75];\n'
        "}"
    )


def test_write_linkgraph_to_dot_one_node_big_nums(tmp_path):
    # checks that 1-node link graphs get written out ok,
    # and that large numbers result in commas included in the node label
    g = nx.Graph()
    g.add_node((12345, 3), ct=56789, freq=0.56789)
    lu.write_linkgraph_to_dot(g, tmp_path, "lonely", 0.01, 5)

    with open(tmp_path / "lonely_linkgraph.gv", "r") as f:
        dot_txt = f.read()
    assert dot_txt == (
        "graph {\n"
        f"{LG_DOT_HEADER_ATTRS}"
        '  "(12345, 3)" [label="12,345 (T)\\n56,789x (56.79%)"];\n'
        "}"
    )


def test_write_linkgraph_to_dot_penwidth_clamp(tmp_path):
    g = nx.Graph()
    g.add_node((12345, 3), ct=56789, freq=0.56789)
    g.add_node((12346, 0), ct=2, freq=1)
    g.add_edge((12345, 3), (12346, 0), link=(2 / 56789))
    lu.write_linkgraph_to_dot(g, tmp_path, "clampy", 0.01, 5)

    with open(tmp_path / "clampy_linkgraph.gv", "r") as f:
        dot_txt = f.read()
    # Verify that the penwidth of the edge is clamped
    assert dot_txt == (
        "graph {\n"
        f"{LG_DOT_HEADER_ATTRS}"
        '  "(12345, 3)" [label="12,345 (T)\\n56,789x (56.79%)"];\n'
        '  "(12346, 0)" [label="12,346 (A)\\n2x (100.00%)"];\n'
        '  "(12345, 3)" -- "(12346, 0)" [penwidth=0.01];\n'
        "}"
    )


def check_c3_linkgraph(g):
    # assumes that we included all possible nodes and edges
    assert len(g.nodes) == 4
    assert set(g.nodes) == set([(7, 0), (7, 3), (8, 1), (8, 3)])
    assert g.nodes[(7, 0)] == {"ct": 7, "freq": 7 / 13}
    assert g.nodes[(7, 3)] == {"ct": 6, "freq": 6 / 13}
    assert g.nodes[(8, 1)] == {"ct": 3, "freq": 3 / 13}
    assert g.nodes[(8, 3)] == {"ct": 10, "freq": 10 / 13}

    assert len(g.edges) == 4
    assert g.edges[(7, 0), (8, 1)] == {"link": 1 / 7}
    assert g.edges[(7, 0), (8, 3)] == {"link": 6 / 10}
    assert g.edges[(7, 3), (8, 1)] == {"link": 2 / 6}
    assert g.edges[(7, 3), (8, 3)] == {"link": 4 / 10}


def test_run_graph_good_verbose_nx(capsys, tmp_path):
    # Write out "inputs", mocking the output of run_nt()
    ndir = tmp_path / "ndir"
    os.makedirs(ndir)

    with open(ndir / f"c3_{POS_FILE_LBL}.pickle", "wb") as f:
        pickle.dump({7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}}, f)

    with open(ndir / f"c3_{POSPAIR_FILE_LBL}.pickle", "wb") as f:
        pickle.dump({(7, 8): {(0, 3): 6, (0, 1): 1, (3, 1): 2, (3, 3): 4}}, f)

    gdir = tmp_path / "gdir"
    lu.run_graph(ndir, 1, 1, 0, "nx", 0.01, 5, gdir, True, mock_log)

    with open(gdir / "c3_linkgraph.pickle", "rb") as f:
        g = pickle.load(f)

    check_c3_linkgraph(g)

    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Going through (co-)occurrence information and "
        "creating link graphs...\n"
        "MockLog: Creating a link graph for contig c3...\n"
        "MockLog: The link graph for contig c3 has 4 node(s).\n"
        "MockLog: The link graph for contig c3 has 4 edge(s).\n"
        'MockLog: Wrote out a link graph (format: "nx") for contig c3.\n'
        "MockLog: Done.\n"
    )


def test_run_graph_empty_verbose_dot(capsys, tmp_path):
    # Write out "inputs", mocking the output of run_nt()
    ndir = tmp_path / "ndir"
    os.makedirs(ndir)

    with open(ndir / f"c3_{POS_FILE_LBL}.pickle", "wb") as f:
        pickle.dump({7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}}, f)

    with open(ndir / f"c3_{POSPAIR_FILE_LBL}.pickle", "wb") as f:
        pickle.dump({(7, 8): {(0, 3): 6, (0, 1): 1, (3, 1): 2, (3, 3): 4}}, f)

    gdir = tmp_path / "gdir"

    # set min_nt_ct to 1,000 -- this prevents any of the nodes from being in
    # the graph, and triggers a unique branch of the code (it doesn't bother
    # trying to create edges and just spits out an empty graph)
    lu.run_graph(ndir, 1000, 1, 0, "dot", 0.01, 5, gdir, True, mock_log)

    with open(gdir / "c3_linkgraph.gv", "r") as f:
        dot = f.read()

    # make sure the graph is empty
    assert dot == ("graph {\n" f"{LG_DOT_HEADER_ATTRS}" "}")

    # special logging message
    # (not a warning b/c this might happen a lot if the user wants to focus on
    # super high cov MAGs for example maybe idk)
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Going through (co-)occurrence information and "
        "creating link graphs...\n"
        "MockLog: Creating a link graph for contig c3...\n"
        "MockLog: The link graph for contig c3 has no nodes, so it's "
        "empty. Try lowering --min-nt-count (it's set to 1,000).\n"
        'MockLog: Wrote out a link graph (format: "dot") for contig c3.\n'
        "MockLog: Done.\n"
    )


def test_run_graph_no_nt_info_at_all(tmp_path):
    ndir = tmp_path / "ndir"
    os.makedirs(ndir)
    # (no need to make gdir exist, because run_graph() will -- or at least, it
    # would if any co-occurrence information existed.)
    gdir = tmp_path / "gdir"

    # First, check that this raises an error.
    with pytest.raises(FileNotFoundError) as ei:
        lu.run_graph(ndir, 1, 1, 0, "nx", 0.01, 5, gdir, True, mock_log)
    assert str(ei.value) == (
        f"Didn't find any (co-)occurrence information in {ndir}."
    )

    # Second, check that this didn't create gdir yet.
    assert not os.path.exists(gdir)


def test_run_graph_only_pos_info_file(tmp_path):
    ndir = tmp_path / "ndir"
    os.makedirs(ndir)

    pf = ndir / f"c3_{POS_FILE_LBL}.pickle"
    ppf = ndir / f"c3_{POSPAIR_FILE_LBL}.pickle"
    with open(pf, "wb") as f:
        pickle.dump({7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}}, f)

    gdir = tmp_path / "gdir"

    with pytest.raises(FileNotFoundError) as ei:
        lu.run_graph(ndir, 1, 1, 0, "nx", 0.01, 5, gdir, True, mock_log)
    assert str(ei.value) == f"Found file {pf}, but not file {ppf}."

    assert not os.path.exists(gdir)


def test_run_graph_only_pospair_info_file(tmp_path):
    # this doesn't trigger a unique error message, since we use the presence of
    # the "position" info file as a sign to look for a "position pair" info
    # file. however, i expect this case should be rare in practice (should only
    # happen if you start deleting stuff from the output of "link nt"...)
    ndir = tmp_path / "ndir"
    os.makedirs(ndir)

    ppf = ndir / f"c3_{POSPAIR_FILE_LBL}.pickle"
    with open(ppf, "wb") as f:
        pickle.dump({(7, 8): {(0, 3): 6, (0, 1): 1, (3, 1): 2, (3, 3): 4}}, f)

    gdir = tmp_path / "gdir"

    with pytest.raises(FileNotFoundError) as ei:
        lu.run_graph(ndir, 1, 1, 0, "nx", 0.01, 5, gdir, True, mock_log)
    assert str(ei.value) == (
        f"Didn't find any (co-)occurrence information in {ndir}."
    )

    assert not os.path.exists(gdir)


def test_run_graph_bad_output_format(tmp_path):
    ndir = tmp_path / "ndir"
    os.makedirs(ndir)

    with open(ndir / f"c3_{POS_FILE_LBL}.pickle", "wb") as f:
        pickle.dump({7: {0: 7, 3: 6}, 8: {1: 3, 3: 10}}, f)

    with open(ndir / f"c3_{POSPAIR_FILE_LBL}.pickle", "wb") as f:
        pickle.dump({(7, 8): {(0, 3): 6, (0, 1): 1, (3, 1): 2, (3, 3): 4}}, f)

    gdir = tmp_path / "gdir"

    with pytest.raises(WeirdError) as ei:
        lu.run_graph(ndir, 1, 1, 0, "Sus", 0.01, 5, gdir, True, mock_log)
    assert str(ei.value) == 'Unrecognized output format: "Sus"'


def check_edge_existence(exp_edges, obs_edge_lines):
    # utility function: check for the presence of "expected edge" 3-tuples of
    # (left node, right node, penwidth) in a set of edge lines in a DOT file.

    # we'll remove stuff from this set, so make a copy first to avoid breaking
    # things unexpectedly for the caller
    obs_edge_lines_copy = copy.deepcopy(obs_edge_lines)

    # this is a terrible inefficient way of doing this check but there should
    # only be a few elements in either set so efficiency doesn't really matter
    for e in exp_edges:
        found_matching_edge_line = False
        for el in obs_edge_lines_copy:
            if str(e[0]) in el and str(e[1]) in el:
                assert el.startswith(f'  "{e[0]}" -- "{e[1]}" [penwidth=')
                # the [:-3] slices off the ];\n (the \n is one character)
                # we're making the assumption that there are no edge attrs
                # besides penwidth, ofc
                penwidth = float(el.split("=")[1][:-3])
                assert penwidth == approx(e[2])
                obs_edge_lines_copy.remove(el)
                found_matching_edge_line = True
                break
        if not found_matching_edge_line:
            raise AssertionError(
                f"Didn't find a matching edge line for edge {e}"
            )
    assert len(obs_edge_lines_copy) == 0


def test_link_graph_integration(tmp_path):
    ndir = tmp_path / "ndir"
    lu.run_nt(FASTA, BAM, BCF, ndir, False, mock_log)

    gdir = tmp_path / "gdir"
    lu.run_graph(ndir, 1, 1, 0, "dot", 0.01, 5, gdir, False, mock_log)

    # there shouldn't be a c2 link graph b/c it has no mutations
    assert sorted(os.listdir(gdir)) == ["c1_linkgraph.gv", "c3_linkgraph.gv"]

    # First, check the c3 graph -- same data as above

    # ignore the order in which we see node / edge lines, to make testing
    # more consistent. i guess we could make this stricter by checking that the
    # observed file at least starts with the graph { line and ends with the }
    # line, but we've already tested the DOT-exporting function above so that
    with open(gdir / "c3_linkgraph.gv", "r") as f:
        obs_dot_lines = set(f.readlines())

    # isn't a priority.
    exp_nonedge_lines = set(
        [
            "graph {\n",
            '  "(7, 0)" [label="7 (A)\\n7x (53.85%)"];\n',
            '  "(7, 3)" [label="7 (T)\\n6x (46.15%)"];\n',
            '  "(8, 1)" [label="8 (C)\\n3x (23.08%)"];\n',
            '  "(8, 3)" [label="8 (T)\\n10x (76.92%)"];\n',
            "}",
        ]
    ) | set([x + "\n" for x in LG_DOT_HEADER_ATTRS.splitlines()])

    assert exp_nonedge_lines.issubset(obs_dot_lines)

    # i guess there's a possibility that newline crap makes this set have 5
    # elements (including an extra line at the end of the file). if so it
    # should be simple to adjust this test to accommodate that (but i don't
    # wannnnna)
    obs_edge_lines = obs_dot_lines - exp_nonedge_lines
    assert len(obs_edge_lines) == 4

    # checking the presence of edges is a bit trickier because penwidths are
    # not formatted in a consistent way (e.g. always 2 digits, like the
    # percentages above). so we do some custom stuff
    exp_edges = [
        ((7, 0), (8, 1), (1 / 7) * 5),
        ((7, 0), (8, 3), 3),
        ((7, 3), (8, 1), (2 / 6) * 5),
        ((7, 3), (8, 3), 2),
    ]

    check_edge_existence(exp_edges, obs_edge_lines)

    # OK now check the c1 graph

    with open(gdir / "c1_linkgraph.gv", "r") as f:
        obs_dot_lines = set(f.readlines())

    exp_nonedge_lines = set(
        [
            "graph {\n",
            '  "(4, 0)" [label="4 (A)\\n1x (8.33%)"];\n',
            '  "(4, 1)" [label="4 (C)\\n2x (16.67%)"];\n',
            '  "(4, 2)" [label="4 (G)\\n6x (50.00%)"];\n',
            '  "(4, 3)" [label="4 (T)\\n3x (25.00%)"];\n',
            '  "(11, 0)" [label="11 (A)\\n5x (41.67%)"];\n',
            '  "(11, 2)" [label="11 (G)\\n7x (58.33%)"];\n',
            '  "(13, 0)" [label="13 (A)\\n5x (41.67%)"];\n',
            '  "(13, 2)" [label="13 (G)\\n7x (58.33%)"];\n',
            "}",
        ]
    ) | set([x + "\n" for x in LG_DOT_HEADER_ATTRS.splitlines()])

    assert exp_nonedge_lines.issubset(obs_dot_lines)

    obs_edge_lines = obs_dot_lines - exp_nonedge_lines
    assert len(obs_edge_lines) == 14

    # see test_run_nt() above for some context on this insanity
    exp_edges = [
        ((4, 0), (11, 0), (1 / 5) * 5),
        ((4, 1), (11, 0), (2 / 5) * 5),
        ((4, 2), (11, 0), (1 / 6) * 5),
        ((4, 2), (11, 2), (5 / 7) * 5),
        ((4, 3), (11, 0), (1 / 5) * 5),
        ((4, 3), (11, 2), (2 / 7) * 5),
        # the edges between the position 4 and 11 nodes exactly match the edges
        # between the position 4 and 13 nodes
        ((4, 0), (13, 0), (1 / 5) * 5),
        ((4, 1), (13, 0), (2 / 5) * 5),
        ((4, 2), (13, 0), (1 / 6) * 5),
        ((4, 2), (13, 2), (5 / 7) * 5),
        ((4, 3), (13, 0), (1 / 5) * 5),
        ((4, 3), (13, 2), (2 / 7) * 5),
        # positions 11 and 13 always match
        ((11, 0), (13, 0), 5),
        ((11, 2), (13, 2), 5),
    ]

    check_edge_existence(exp_edges, obs_edge_lines)
