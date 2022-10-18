# Utilities for dealing with assembly graph files (at the moment, just GFA).

import networkx as nx
from .errors import GraphParsingError, WeirdError


def load_gfa(gfa_fp, min_num_nodes=1):
    """Quick and dirty function that loads GFA files in NetworkX.

    GfaPy was being really slow on the SheepGut graph, which is why I wrote
    this.

    Parameters
    ----------
    gfa_fp: str
        Filepath to a GFA 1 file. Eventually I should ditch this function and
        just use Gfapy to load this graph, but
        https://github.com/ggonnella/gfapy/issues/25 is blocking that.

    min_num_nodes: int
        Minimum number of segments that must be described in the GFA file.
        We'll raise an error if this condition is not met.

    Returns
    -------
    graph: nx.Graph
        Undirected graph representing the topology of the graph contained in
        the input GFA file. Note that we ignore segment sequences (only storing
        information about the nodes' lengths), and that we ignore
        directionality entirely.

    Raises
    ------
    GraphParsingError
        If any segment has either no length given for it (sequence is * and the
        LN tag is missing) or multiple lengths given for it (sequence is
        defined, and there is an LN tag given).

        The latter of these conditions is a biiit overly strict, since -- for
        some arbitrary segment S -- if a sequence is defined for S, and LN is
        given for s, and if these lengths match up, then this is still
        understandable. But I'm pessimistic, and if we accept this then
        eventually we're gonna start seeing cases where the two lengths
        disagree, and that way lies madness.

        Also raised if the graph has less than min_num_nodes nodes (segments).
    """

    # We ignore directionality for right now.
    graph = nx.Graph()

    num_nodes = 0
    with open(gfa_fp, "r") as gfafile:
        for line in gfafile:

            if line[0] == "S":
                parts = line.strip().split("\t")
                node_name = parts[1]
                if parts[2] == "*":
                    node_len = None
                else:
                    node_len = len(parts[2])

                # Parse GFA tags, only caring (for now, at least) about length
                extra_data = parts[3:]
                for tag in extra_data:
                    if tag.startswith("LN:i:"):
                        if node_len is None:
                            node_len = int(tag[5:])
                        else:
                            # Initially, I always threw an error in this case,
                            # but some real assembly graphs (e.g. hifiasm-meta
                            # output, as of writing) have both LN tags and
                            # sequences given. So now I only throw an error if
                            # the multiple sources of information on length
                            # disagree.
                            tag_len = int(tag[5:])
                            if tag_len != node_len:
                                raise GraphParsingError(
                                    "Contradictory lengths found for segment "
                                    f"{node_name}: sequence length is "
                                    f"{node_len:,} bp, but the LN tag says "
                                    f"the length is {tag_len:,} bp."
                                )

                if node_len is None:
                    raise GraphParsingError(
                        f"No length given for segment {node_name}"
                    )

                if node_name in graph.nodes:
                    raise GraphParsingError(
                        f"Segment {node_name} is defined multiple times"
                    )

                graph.add_node(node_name, length=node_len)
                num_nodes += 1

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

    if num_nodes < min_num_nodes:
        raise GraphParsingError(
            f"Less than {min_num_nodes:,} segment(s) are given in {gfa_fp}."
        )
    return graph


def gfa_to_fasta(gfa_fp, fasta_file, chunk_size=500):
    """Converts a GFA file to a FASTA file.

    Parameters
    ----------
    gfa_fp: str
        Filepath to a GFA 1 file.

    fasta_file: typing.TextIO
        Opened IO stream to which an output FASTA file (containing the
        sequences of the segments in the input GFA file) will be written.

        To clarify, this should be something like an opened file handle
        or an opened io.StringIO object. (The main use-case for the io.StringIO
        thing is to simplify testing.)

        See https://stackoverflow.com/a/38569536 for details on this
        sort of abstraction.

    chunk_size: int
        This should be a positive integer. After seeing this many segments,
        we'll write out their sequences to fasta_file. This way, we avoid
        storing every sequence in the GFA file in memory at once.

        Note that there's a tradeoff here, since writing things out to the
        FASTA file too frequently will slow things down. The default should be
        reasonable-ish, but feel free to adjust it if desired.

    Returns
    -------
    num_seqs: int
        The number of segments (sequences) contained in the input GFA (output
        FASTA) file.

    Raises
    ------
    GraphParsingError
        If any segment has no sequence given (i.e. the provided sequence is *),
        or if there are no segments contained in the GFA file.

    WeirdError
        If chunk_size is less than 1.
    """
    if chunk_size < 1:
        raise WeirdError(
            f"chunk_size should be a positive integer, but it's {chunk_size}?"
        )
    num_seqs = 0
    fout = ""
    with open(gfa_fp, "r") as gfa_file:
        for line in gfa_file:
            if line.startswith("S\t"):
                split = line.strip().split("\t")
                seq = split[2]
                if seq == "*":
                    raise GraphParsingError(
                        f"No sequence given for segment {split[1]}"
                    )
                fout += f">{split[1]}\n{split[2]}\n"
                num_seqs += 1
                if num_seqs % chunk_size == 0:
                    fasta_file.write(fout)
                    fout = ""
    if fout != "":
        fasta_file.write(fout)

    if num_seqs == 0:
        raise GraphParsingError("Didn't see any segments in the GFA file?")

    return num_seqs
