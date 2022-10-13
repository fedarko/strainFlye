import os
import pickle
import pysam
import pytest
import strainflye.link_utils as lu
from collections import defaultdict
from strainflye.config import POS_FILE_LBL, POSPAIR_FILE_LBL
from strainflye.errors import ParameterError
from strainflye.tests.utils_for_testing import mock_log


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BCF = os.path.join(IN_DIR, "call-r-min3-di12345", "naive-calls.bcf")
BAM = os.path.join(IN_DIR, "alignment.bam")
DEGEN_BAM = os.path.join(IN_DIR, "degen.bam")
SECONDARY_BAM = os.path.join(IN_DIR, "c4-and-secondary.bam")


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


def test_run_nt(capsys, tmp_path):
    output_dir = os.path.join(tmp_path, "run-nt-out")
    lu.run_nt(FASTA, BAM, BCF, output_dir, True, mock_log)

    fp_c1p = os.path.join(output_dir, f"c1_{POS_FILE_LBL}.pickle")
    fp_c1pp = os.path.join(output_dir, f"c1_{POSPAIR_FILE_LBL}.pickle")
    fp_c3p = os.path.join(output_dir, f"c3_{POS_FILE_LBL}.pickle")
    fp_c3pp = os.path.join(output_dir, f"c3_{POSPAIR_FILE_LBL}.pickle")

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
            # (T, A)
            (3, 0): 1,
            # (T, G)
            (3, 2): 2,
            # (G, G)
            (2, 2): 5,
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