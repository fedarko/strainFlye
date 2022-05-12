import pysam
from strainflye.align_utils import get_coords


def test_get_coords():
    samfile = pysam.AlignmentFile("strainflye/tests/inputs/camp-sf.sam", "r")
    for aln in samfile.fetch():
        assert get_coords(aln) == (1162486, 1259415)
