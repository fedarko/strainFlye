# Utilities for strainFlye's alignment step.

import pysam
import subprocess
from itertools import combinations
from collections import defaultdict


def index_bam(in_bam, bam_descriptor, fancylog):
    """Indexes a BAM file using samtools.

    This creates a .bam.bai file in the same location as the BAM file.

    Parameters
    ----------
    in_bam: str
        Location of the BAM file to be indexed.

    bam_descriptor: str
        Description of the in_bam file, to be used in the logging message.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    fancylog(f"Indexing the {bam_descriptor}...")
    subprocess.run(["samtools", "index", in_bam])
    fancylog(f"Done indexing the {bam_descriptor}.")


def get_quartet(alnseg):
    """Returns a tuple of information distinguishing a linear alignment.

    Parameters
    ----------
    alnseg: pysam.AlignedSegment

    Returns
    -------
    (s, e, mq, st): (int, int, int, str)
        Segment start, end, mapping quality, and to_string() output.

        The start and end are both inclusive, to simplify comparison of
        alignment ranges for detecting overlaps.

        The reason we include the fourth element (to_string()) is to make
        it easier to distinguish linear alignments from the same read. It
        is very unlikely (but still possible I guess) that multiple
        alignments from a read will have identical QUAL values AND
        identical tags, both of which are included in to_string() as of
        writing.

    Raises
    ------
    ValueError
        If the segment's start is greater than its end (both in inclusive
        coordinates). (If this ends up being a problem in practice, maybe
        because there of reverse-mapped reads or something (???), then this
        could probs be modified to just reverse the start and end in this
        case.)
    """
    # We add 1 to the end since this is a half-open interval -- we want
    # the coordinates we use for computing overlap to be completely
    # inclusive intervals
    s = linearaln.reference_start
    e = linearaln.reference_end + 1
    if s > e:
        raise ValueError(
            f"Malformed linear alignment coordinates: start {s}, end {e}"
        )
    mq = linearaln.mapping_quality
    st = linearaln.to_string()
    return (s, e, mq, st)


def filter_osa_reads(in_bam, out_bam, fancylog):
    """Filters a BAM file to remove all reads with OSAs.

    A read has an OSA (overlapping supplementary alignment) if any two
    supplementary alignments of this read are aligned to coordinates on a
    reference sequence that overlap.

    This assumes that this BAM file does not contain secondary alignments.
    Supplementary alignments should have been left in, but secondary alignments
    will mess this up!

    Parameters
    ----------
    in_bam: str
        Location of the BAM file on which we will perform filtering. This BAM
        file should already have been sorted and indexed.

    out_bam: str
        Location to which an output BAM file (containing no reads with OSAs)
        will be written.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    bf = pysam.AlignmentFile(in_bam, "rb")

    # Keeps track of the names of reads with OSAs.
    # We'll filter these reads out completely, so that they are not represented
    # in any of the alignments remaining in the BAM file.
    #
    # I was initially going to maintain a sorted array for this, so that we
    # could use binary search to quickly check if reads seen in the final
    # "pass" over the BAM file had OSAs, but it turns out that in Python using
    # sets is probably a better idea (or at the very least a simpler one):
    # https://stackoverflow.com/a/212971
    reads_with_osa = set()

    for si, seq in enumerate(bf.references, 1):
        pct = 100 * (si / bf.nreferences)
        fancylog(
            f"OSA filter pass 1/2: on seq {seq} ({si:,} / {bf.nreferences:,}) "
            f"({pct:.2f}%)."
        )

        # Identify all linear alignments of each read to this sequence
        readname2CoordsAndMQ = defaultdict(list)
        num_linear_alns = 0
        for linearaln in bf.fetch(seq):
            rn = linearaln.query_name
            alndetails = get_quartet(linearaln)
            if alndetails in readname2CoordsAndMQ[rn]:
                raise ValueError(
                    f"Indistinguishable linear alignments to seq {seq} with read "
                    f"name {rn}: multiple reads share (start, end, mapq, "
                    f"to_string()) of {alndetails}"
                )
            readname2CoordsAndMQ[rn].append(alndetails)
            num_linear_alns += 1

        # The number of unique reads is just the number of keys in this dict
        num_reads = len(readname2CoordsAndMQ)

        # print(f"\t{num_reads:,} unique read(s), {num_linear_alns:,} linear aln(s)...")

        # Identify overlapping alignments from the same read
        n_reads_w_osa_in_seq = 0
        for rn in readname2CoordsAndMQ:
            alns = readname2CoordsAndMQ[rn]
            if len(alns) > 1:
                # Okay, so this particular read has multiple supplementary
                # alignments to this sequence. Check if they overlap.

                for (a1, a2) in combinations(alns, 2):
                    # Efficiently test for overlap between two ranges:
                    # https://stackoverflow.com/a/3269471
                    if a1[0] <= a2[1] and a2[0] <= a1[1]:
                        # Okay, these two alignments of this read overlap.
                        reads_with_osa.add(rn)
                        n_reads_w_osa_in_seq += 1
                        break

        fancylog(
            f"{n_reads_w_osa_in_seq:,} read(s) with OSA(s) in seq {seq}..."
        )

    # Now we've made note of all reads with OSAs across *all* sequences in the
    # alignments. We can make another pass through and output all reads without
    # OSAs into a new BAM file.

    # Output BAM file (filtered to remove reads with overlapping supplementary
    # alignments, aka OSAs)
    of = pysam.AlignmentFile(out_bam, "wb", template=bf)

    # TODO: maybe generalize this iteration code into a generator or something
    # to limit code reuse
    for si, seq in enumerate(bf.references, 1):
        pct = 100 * (si / bf.nreferences)
        fancylog(
            f"OSA filter pass 2/2: on seq {seq} ({si:,} / {bf.nreferences:,}) "
            f"({pct:.2f}%)."
        )

        num_alns_retained = 0
        num_alns_filtered = 0
        for linearaln in bf.fetch(seq):
            rn = linearaln.query_name
            # If this read has OSAs anywhere in the alignment, don't include it in
            # the output BAM file. Otherwise, *do* include it!
            if rn in reads_with_osa:
                num_alns_filtered += 1
            else:
                of.write(linearaln)
                num_alns_retained += 1

        num_alns_total = num_alns_retained + num_alns_filtered
        if num_alns_total > 0:
            apct = 100 * (num_alns_retained / num_alns_total)
        else:
            apct = float("inf")
        fancylog(
            f"{num_alns_retained:,} / {num_alns_total:,} ({apct:.2f}%) "
            "linear aln(s) retained in seq {seq}."
        )

    bf.close()
    of.close()

    fancylog("Done filtering out overlapping supplementary alignments.")


###############################################################################


def filter_pm_reads():
    pass


#! /usr/bin/env python3
# Filters reads that are less than some percentage (MIN_PERCENT_ALIGNED)
# aligned to an edge in the graph, or to other edges adjacent to this edge in
# the repeat graph.
#
# Realized after the fact that this bears some resemblance to samclip
# (https://github.com/tseemann/samclip), although this differs a bit in the
# sort of alignments this allows to pass the filter and the sort of information
# it takes into account. These are probs ultimately minor distinctions, tho.

import re
import pysam
from collections import defaultdict
from utils import load_gfa

# This is a percentage (value in [0, 1]). Reads where less than this
# percentage of the sequence is aligned to a given edge are filtered out of
# the BAM file.
MIN_PERCENT_ALIGNED = 0.9

print(
    "Filtering out reads mapped to other edges and partially-mapped reads..."
)

# flush stdout to force this to show up even though it's not a newline
print("Loading graph...", end=" ", flush=True)
with open("../config/input-graph", "r") as igfile:
    graph_filename = next(igfile).strip()

graph = load_gfa(graph_filename)
print("Done.")

# Input BAM file (contains alignments to all edges in the graph)
bf = pysam.AlignmentFile(
    "output/overlap-supp-aln-filtered-and-sorted-aln.bam", "rb"
)
# Output BAM file (will just contain alignments that pass the filter)
of = pysam.AlignmentFile("output/pmread-filtered-aln.bam", "wb", template=bf)

# Per the SAM v1 specification, top of page 8:
# "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ."
#
# These are the five CIGAR operations that consume character(s) from the query
# sequence (i.e. a read).
#
# We ignore S (soft clipping), since this indicates that a given position in a
# read is not matched to an edge; and we ignore I (insertion), since this also
# indicates that a position in a read is not really "matched" anywhere on the
# edge (although I guess you could argue that including insertions here may be
# useful; it probably isn't a big deal either way).
#
# By just looking at M, X, and = occurrences, we can get a count for each
# alignment of the number of bases "(mis)matched" to an edge.
# Since we have already filtered secondary alignments and overlapping
# supplementary alignments, we can sum these match counts across all
# alignments of a read and divide by the read length to get the approx
# percentage (see above for slight caveats) of bases in the read aligned to
# an edge or group of edges.
matches = re.compile("(\d+)[MX=]")


def check_and_update_alignment(aln, readname2len, readname2matchct, edge_name):
    """Updates readname2len and readname2matchct based on an AlignedSegment.

    Note that readname2len and readname2matchct are both defined relative to
    a single edge in the graph (these aren't universal structures).
    """
    # Ensure that read length is consistent across all alignments involving
    # this read; also, make a record of previously unseen reads' lengths.
    # (Apparently query length is dependent on the actual alignment of this
    # read, so we use infer_read_length() instead because we care about the
    # actual length of the read)
    if aln.query_name in readname2len:
        if readname2len[aln.query_name] != aln.infer_read_length():
            raise ValueError(
                "Inconsistent read lengths across alignments: {}: {}, {}".format(
                    aln.query_name,
                    readname2len[aln.query_name],
                    aln.infer_read_length(),
                )
            )
    else:
        readname2len[aln.query_name] = aln.infer_read_length()

    # Each AlignedSegment returned by fetch(edge) should pertain to that
    # specific edge sequence -- this lets know that the CIGAR string of
    # this segment applies to the edge
    if aln.reference_name != edge_name:
        raise ValueError(
            "Alignment reference name, {}, isn't {} as expected".format(
                aln.reference_name, edge_name
            )
        )

    # The meat of this -- parse the CIGAR string of this alignment and
    # count all (mis)match operations, updating a defaultdict.
    allmatches = matches.findall(aln.cigarstring)
    if allmatches:
        matchct = sum([int(c) for c in allmatches])
        readname2matchct[aln.query_name] += matchct
    else:
        # Raise an error if an alignment of this read does not involve
        # any (mis)match operations at all. This *could* happen in practice,
        # I guess, but if it does something is likely wrong. If this check
        # needs to be removed in the future, then this block could just be
        # commented out or replaced with a "pass" statement or something.
        raise ValueError(
            "No match chars (M/X/=) found in read {} CIGAR: {}".format(
                aln.query_name, aln.cigarstring
            )
        )


# Figure out all reads that are aligned to each edge to focus on
for edge_to_focus_on in graph.nodes():
    # The loaded graph doesn't use the "edge_" prefix for nodes, but the
    # alignment does -- add back in and use it wherever we have to deal with
    # pysam here
    focused_seq = f"edge_{edge_to_focus_on}"
    print(f"Looking at edge {edge_to_focus_on}...")
    # Maps read name to read length (which should be constant across all
    # alignments of that read). This variable is used both to store this info
    # (which is in turn used for sanity checking) as well as as a
    # crude indication of "have we seen this read yet?"
    readname2len = {}

    # Maps read name to number of match operations to the sequences of this
    # edge or adjacent edges in the graph.
    readname2matchct = defaultdict(int)

    for i, aln in enumerate(bf.fetch(focused_seq), 1):
        check_and_update_alignment(
            aln, readname2len, readname2matchct, focused_seq
        )
    print(f"\t{i} linear alignments to this edge.")

    # Identify adjacent edges to this one in the graph, if present. Allow
    # alignments to these edges to count towards readname2matchct (so we can
    # somewhat avoid penalizing edges that aren't in isolated components).
    # We could also expand this to include _all_ other edges in this edge's
    # component, but that will probs cause problems with hairball
    # component(s) in the graph that span thousands of edges...
    adj_edges = set(graph.neighbors(edge_to_focus_on)) - set(
        [edge_to_focus_on]
    )
    nae = len(adj_edges)

    print(f"\t{nae} adjacent edges in the graph.")

    # To prevent this script from taking a super long amount of time, only
    # allow alignments to adjacent edges if there are < 50 adjacent edges to
    # this edge in the graph. This isn't ideal, I guess, but it's better than
    # the script running forever (and we're already trying to care a lot about
    # the small details here by caring about the graph at all in the alignment
    # filtering process).
    #
    # NOTE: could probs speed this up by looking at SA: tags of reads from
    # minimap2 (although we'd still need to perform fetch to get accurate
    # CIGAR strings of those alignments; see
    # https://github.com/lh3/minimap2/issues/724)
    if nae >= 50:
        print("\tToo many adj. edges; we won't look at their alignments here.")
    elif nae <= 0:
        print("\tNothing to consider in adjacent edges; moving on.")
    else:
        # nae is > 0 and < 50, so we can look at adjacent edges
        print("\tLooking at alignments of shared reads to adjacent edges...")
        # Go through these "allowed" other edges; add to the number
        # of matches in cc for any reads that we see that we've already seen in
        # the edge to focus on. (We implicitly ignore any reads aligned to
        # these edges but not to the edge to focus on.)
        num_other_edge_alns_from_shared_reads = 0
        for other_edge in adj_edges:
            other_edge_seq = f"edge_{other_edge}"
            for aln in bf.fetch(other_edge_seq):
                # If aln.query_name is NOT in readname2len, then this
                # alignment's read wasn't also aligned to the edge to focus on
                # -- in this case we implicitly ignore it, as mentioned above,
                # since we only really care about reads aligned to the edge
                # we're focusing on.
                if aln.query_name in readname2len:
                    check_and_update_alignment(
                        aln, readname2len, readname2matchct, other_edge_seq
                    )
                    num_other_edge_alns_from_shared_reads += 1

        print(
            f"\t{num_other_edge_alns_from_shared_reads} linear alns from "
            "shared reads to adjacent edges."
        )

    # Now that we've considered all relevant edges, we can compute
    # the approximate percentages of each read (aligned to the edge we're
    # focusing on) aligned to all edges in this component.
    #
    # And, finally, we can selectively write these reads' alignments
    # to an output BAM file accordingly.
    #
    # Note that although we're iterating over all alignments, these
    # computations are identical for alignments from the same read (at least in
    # the context of the same focused_seq). So if a read passes the
    # filter for an edge, all its alignments to this edge will be output to the
    # BAM file; and if the read fails the filter for an edge, none of its
    # alignments to this edge will be output to the BAM file (although this
    # does not preclude other alignments from this read to other edges from
    # being output in the context of other edges in this script).
    print(
        "\tComputing percentages and outputting alignments from "
        "passing reads..."
    )
    p = 0
    for aln in bf.fetch(focused_seq):
        if aln.query_name not in readname2len:
            raise ValueError(
                f"We should have seen read {aln.query_name} earlier!"
            )
        read_len = readname2len[aln.query_name]

        if readname2matchct[aln.query_name] == 0:
            # As with a similar error case in check_and_update_alignment(),
            # I *guess* this could happen in practice but it is probably
            # indicative of an error in most cases.
            raise ValueError(
                f"Read {aln.query_name} had no match operations done?"
            )

        read_matchct_in_cc = readname2matchct[aln.query_name]
        perc = read_matchct_in_cc / read_len

        if perc >= MIN_PERCENT_ALIGNED:
            of.write(aln)
            p += 1
    print(
        f"\t{p} / {i} ({((p / i)*100):.2f}%) of alignments in "
        f"edge {edge_to_focus_on} passed the filter."
    )

bf.close()
of.close()

print("Filtered to fully (ish) aligned reads.")
