import os
import pysam
import pytest
import strainflye.align_utils as au
from strainflye.errors import SequencingDataError, GraphParsingError

TI_DIR = os.path.join("strainflye", "tests", "inputs")
S1 = os.path.join(TI_DIR, "sample1.gfa")

def test_get_coords_good():
    samfile = pysam.AlignmentFile(os.path.join(TI_DIR, "camp-sf.sam"), "r")
    for aln in samfile.fetch():
        assert au.get_coords(aln) == (1162486, 1259415)


def test_check_contigs_in_graph_good():
    au.check_contigs_in_graph(
        {"1": 8, "2": 10, "3": 21, "4": 7, "5": 8, "6": 4},
        S1
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
            au.check_contigs_in_graph( n2l, S1)
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
