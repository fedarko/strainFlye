import networkx as nx


def load_gfa(filepath):
    """Quick and dirty function that loads GFA files in NetworkX.

    GfaPy was being really slow on the SheepGut graph, which is why I wrote
    this.
    """

    # We ignore directionality for right now.
    graph = nx.Graph()

    with open(filepath, "r") as gfafile:
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
