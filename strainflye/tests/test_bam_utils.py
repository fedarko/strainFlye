import os
import pysam
import strainflye.bam_utils as bu

TI_DIR = os.path.join("strainflye", "tests", "inputs")


def test_get_coords_good():
    # TODO FIX
    samfile = pysam.AlignmentFile(os.path.join(TI_DIR, "camp-sf.sam"), "r")
    for aln in samfile.fetch():
        assert bu.get_coords(aln) == (1162486, 1259415)
