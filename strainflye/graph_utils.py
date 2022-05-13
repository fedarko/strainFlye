import networkx as nx


def load_gfa(gfa_fp):
    """Quick and dirty function that loads GFA files in NetworkX.

    GfaPy was being really slow on the SheepGut graph, which is why I wrote
    this.

    Parameters
    ----------
    gfa_fp: str
        Filepath to a GFA 1 file. Eventually I should ditch this function and
        just use Gfapy to load this graph, but
        https://github.com/ggonnella/gfapy/issues/25 is blocking that.

    Returns
    -------
    graph: nx.Graph
        Undirected graph representing the topology of the graph contained in
        the input GFA file. Note that we ignore segment sequences (only storing
        information about the nodes' lengths), and that we ignore
        directionality entirely.

    Raises
    ------
    ValueError
        If any segment has either no length given for it (sequence is * and the
        LN tag is missing) or multiple lengths given for it (sequence is
        defined, and there is an LN tag given).

        The latter of these conditions is a biiit overly strict, since -- for
        some arbitrary segment S -- if a sequence is defined for S, and LN is
        given for s, and if these lengths match up, then this is still
        understandable. But I'm pessimistic, and if we accept this then
        eventually we're gonna start seeing cases where the two lengths
        disagree, and that way lies madness.
    """

    # We ignore directionality for right now.
    graph = nx.Graph()

    with open(gfa_fp, "r") as gfafile:
        for line in gfafile:

            if line[0] == "S":
                parts = line.strip().split("\t")
                node_name = parts[1]
                if parts[2] == "*":
                    node_len = None
                else:
                    node_len = len(parts[2])

                # Parse GFA tags, only caring about length
                extra_data = parts[3:]
                for tag in extra_data:
                    # This is overcautious, probably, but whatevs
                    if tag.startswith("LN:i:"):
                        if node_len is None:
                            node_len = int(tag[5:])
                        else:
                            raise ValueError(
                                "Duplicate length for node {node_name}"
                            )

                if node_len is None:
                    raise ValueError("No length given for node {node_name}.")

                graph.add_node(node_name, length=node_len)

            # Parse links between sequences.
            # Each link line looks something like:
            #
            # L    edge_74133    -    edge_71431    -    0M
            #
            # (where the whitespace gaps are each a single tab character)
            elif line[0] == "L":
                parts = line.strip().split("\t")
                src = parts[1]
                snk = parts[3]
                # As mentioned, we ignore directionality for this specific
                # application. So, e.g.:
                # A+ -> B+
                # A+ -> B-
                # A- -> B+
                # A- -> B-
                # ... would all get treated as the same edge. "Duplicate" edges
                # are implicitly ignored by networkx.
                graph.add_edge(src, snk)
    return graph
