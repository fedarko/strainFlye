import pytest
import networkx as nx
import strainflye.graph_utils as gu
from io import StringIO


def check_sample1_graph(g):
    # this should be an undirected graph!
    assert type(g) == nx.Graph
    assert sorted(g.nodes()) == ["1", "2", "3", "4", "5", "6"]
    # This is a pretty brittle way of doing this but it works
    assert sorted(g.edges()) == [
        ("1", "2"),
        ("2", "3"),
        ("3", "4"),
        ("4", "5"),
    ]
    assert g.nodes["1"]["length"] == 8
    assert g.nodes["2"]["length"] == 10
    assert g.nodes["3"]["length"] == 21
    assert g.nodes["4"]["length"] == 7
    assert g.nodes["5"]["length"] == 8
    assert g.nodes["6"]["length"] == 4


def test_load_gfa():
    g = gu.load_gfa("strainflye/tests/inputs/sample1.gfa")
    check_sample1_graph(g)


def test_load_gfa_noseq():
    g = gu.load_gfa("strainflye/tests/inputs/sample1-noseq.gfa")
    check_sample1_graph(g)


def test_load_gfa_duplen():
    # based on previous test code I wrote at
    # https://github.com/marbl/MetagenomeScope/blob/master/metagenomescope/tests/assembly_graph_parser/utils.py
    with pytest.raises(gu.GraphParsingError) as errorinfo:
        gu.load_gfa("strainflye/tests/inputs/sample1-duplen.gfa")
    assert "Duplicate length for segment 3" == str(errorinfo.value)


def test_load_gfa_nolen():
    with pytest.raises(gu.GraphParsingError) as errorinfo:
        gu.load_gfa("strainflye/tests/inputs/sample1-nolen.gfa")
    assert "No length given for segment 4" == str(errorinfo.value)


def test_gfa_to_fasta():
    """A 'normal' test case for this."""
    sio = StringIO()
    gu.gfa_to_fasta("strainflye/tests/inputs/sample1.gfa", sio)

    assert sio.getvalue() == (
        ">1\nCGATGCAA\n"
        ">2\nTGCAAAGTAC\n"
        ">3\nTGCAACGTATAGACTTGTCAC\n"
        ">4\nTATATGC\n"
        ">5\nCGATGATA\n"
        ">6\nATGA\n"
    )


def test_gfa_to_fasta_noseq():
    sio = StringIO()
    with pytest.raises(gu.GraphParsingError) as errorinfo:
        gu.gfa_to_fasta("strainflye/tests/inputs/sample1-noseq.gfa", sio)

    # All of the segments in this particular test GFA have no sequences given,
    # so the first one triggering this error should be the one mentioned in the
    # message. ... Arguably, this error message should be comprehensive, but
    # that's a lot of work and probably not worth the hassle.
    assert "No sequence given for segment 1" == str(errorinfo.value)

    # Nothing should have been written to the StringIO representing the output
    # FASTA file yet
    assert sio.getvalue() == ""


def test_gfa_to_fasta_nosegs():
    sio = StringIO()
    with pytest.raises(gu.GraphParsingError) as errorinfo:
        gu.gfa_to_fasta("strainflye/tests/inputs/sample1-empty.gfa", sio)

    assert "Didn't see any segments in the GFA file?" == str(errorinfo.value)
    assert sio.getvalue() == ""
