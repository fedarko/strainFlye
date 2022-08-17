import tempfile
import pytest
import networkx as nx
import strainflye.graph_utils as gu
from io import StringIO
from strainflye.errors import GraphParsingError


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
    with pytest.raises(GraphParsingError) as errorinfo:
        gu.load_gfa("strainflye/tests/inputs/sample1-duplen.gfa")
    assert "Duplicate length for segment 3" == str(errorinfo.value)


def test_load_gfa_nolen():
    with pytest.raises(GraphParsingError) as errorinfo:
        gu.load_gfa("strainflye/tests/inputs/sample1-nolen.gfa")
    assert "No length given for segment 4" == str(errorinfo.value)


def test_load_gfa_nonunique_seq():
    # Based loosely on
    # https://github.com/marbl/MetagenomeScope/blob/master/metagenomescope/tests/assembly_graph_parser/utils.py
    fh = tempfile.NamedTemporaryFile(suffix=".gfa")
    # There is probably a way to write directly to the filehandler, but idk
    # how and this is easiest
    with open(fh.name, "w") as f:
        f.write(
            "H\tVN:Z:1.0\n"
            "S\t1\t*\tLN:i:8\n"
            "S\t2\t*\tLN:i:10\n"
            "S\t1\t*\tLN:i:9\n"
        )

    with pytest.raises(GraphParsingError) as errorinfo:
        gu.load_gfa(fh.name)
    assert "Segment 1 is defined multiple times" == str(errorinfo.value)


def test_load_gfa_empty():
    gf = "strainflye/tests/inputs/sample1-empty.gfa"
    with pytest.raises(GraphParsingError) as errorinfo:
        gu.load_gfa(gf)
    assert f"Less than 1 segment(s) are given in {gf}." == str(errorinfo.value)


def test_load_gfa_min_num_nodes():
    gf = "strainflye/tests/inputs/sample1.gfa"
    for mnn in [7, 8, 9, 10, 14, 20, 1000, 10000]:
        with pytest.raises(GraphParsingError) as ei:
            gu.load_gfa(gf, min_num_nodes=mnn)
        assert f"Less than {mnn:,} segment(s) are given in {gf}." == str(
            ei.value
        )
    for mnn in range(1, 7):
        g = gu.load_gfa(gf, min_num_nodes=mnn)
        check_sample1_graph(g)


def check_sample1_fasta(fasta_text, num_seqs):
    assert fasta_text == (
        ">1\nCGATGCAA\n"
        ">2\nTGCAAAGTAC\n"
        ">3\nTGCAACGTATAGACTTGTCAC\n"
        ">4\nTATATGC\n"
        ">5\nCGATGATA\n"
        ">6\nATGA\n"
    )

    assert num_seqs == 6


def test_gfa_to_fasta():
    """A 'normal' test case for this."""
    sio = StringIO()
    num_seqs = gu.gfa_to_fasta("strainflye/tests/inputs/sample1.gfa", sio)
    check_sample1_fasta(sio.getvalue(), num_seqs)


def test_gfa_to_fasta_smallchunksize():
    """Verifies that using a small chunk size works ok."""
    sio = StringIO()
    num_seqs = gu.gfa_to_fasta(
        "strainflye/tests/inputs/sample1.gfa", sio, chunk_size=2
    )
    check_sample1_fasta(sio.getvalue(), num_seqs)


def test_gfa_to_fasta_noseq():
    sio = StringIO()
    with pytest.raises(GraphParsingError) as errorinfo:
        gu.gfa_to_fasta("strainflye/tests/inputs/sample1-noseq.gfa", sio)

    # All of the segments in this particular test GFA have no sequences given,
    # so the first one triggering this error should be the one mentioned in the
    # message. ... Arguably, this error message should be comprehensive, but
    # that's a lot of work and probably not worth the hassle.
    assert "No sequence given for segment 1" == str(errorinfo.value)

    # Nothing should have been written to the StringIO representing the output
    # FASTA file yet -- however, since we write to the FASTA file periodically,
    # this behavior is only guaranteed here because segment 1 is the first one
    # given in the GFA file
    assert sio.getvalue() == ""


def test_gfa_to_fasta_nosegs():
    sio = StringIO()
    with pytest.raises(GraphParsingError) as errorinfo:
        gu.gfa_to_fasta("strainflye/tests/inputs/sample1-empty.gfa", sio)

    assert "Didn't see any segments in the GFA file?" == str(errorinfo.value)
    assert sio.getvalue() == ""
