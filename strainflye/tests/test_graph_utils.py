import networkx as nx
from strainflye.graph_utils import load_gfa


def test_load_gfa():
    g = load_gfa("strainflye/tests/inputs/sample1.gfa")
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
