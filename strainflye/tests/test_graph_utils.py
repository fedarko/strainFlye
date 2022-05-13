import pytest
import networkx as nx
import strainflye.graph_utils as gu


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
    assert "Duplicate length for segment 3" in str(errorinfo.value)


def test_load_gfa_nolen():
    with pytest.raises(gu.GraphParsingError) as errorinfo:
        gu.load_gfa("strainflye/tests/inputs/sample1-nolen.gfa")
    assert "No length given for segment 4" in str(errorinfo.value)
