import os
import pysam
import pytest
import strainflye.link_utils as lu
from strainflye.errors import ParameterError


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
BAM = os.path.join(IN_DIR, "alignment.bam")
DEGEN_BAM = os.path.join(IN_DIR, "degen.bam")
SECONDARY_BAM = os.path.join(IN_DIR, "c4-and-secondary.bam")


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
