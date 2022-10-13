import os
import pysam
import strainflye.link_utils as lu


# import stat
# import gzip
# import tempfile
# import contextlib
# import pytest
# import skbio
# from io import StringIO
# from strainflye.config import DI_PREF, DEFAULT_LJA_PARAMS
# from strainflye.errors import ParameterError, WeirdError
# from strainflye.tests.utils_for_testing import mock_log, mock_log_2


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
# FASTA = os.path.join(IN_DIR, "contigs.fasta")
BAM = os.path.join(IN_DIR, "alignment.bam")
# BCF = os.path.join(IN_DIR, "call-r-min3-di12345", "naive-calls.bcf")


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
