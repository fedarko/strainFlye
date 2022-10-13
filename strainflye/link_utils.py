# Utilities for strainFlye link.


import os
import pickle
import networkx as nx
from itertools import combinations
from collections import defaultdict
from strainflye import bcf_utils, cli_utils, misc_utils, config
from strainflye.errors import WeirdError, ParameterError


def gen_ddi():
    """Returns a new defaultdict(int).

    Needed because pickle can't handle lambda functions.
    """
    return defaultdict(int)


def get_readname2pos2nt(bam_obj, contig, positions):
    """Returns an object mapping read name -> position -> nucleotide.

    Because we should have already filtered secondary alignments and
    overlapping supplementary alignments from the alignment, we are guaranteed
    that an arbitrary read -- across all its linear alignment(s) to a contig --
    covers each position in this contig at most once. So we don't have to worry
    about this read implying multiple distinct nucleotides at a position.

    Parameters
    ----------
    bam_obj: pysam.AlignmentFile
        Object describing a BAM file.

    contig: str
        Name of a contig for which we will compute the mapping.

    positions: list
        List of (zero-indexed) "positions of interest" in this contig for which
        we will retrieve nucleotide information. (In practice, this will
        probably be a list of mutated positions.)

    Returns
    -------
    readname2pos2nt: defaultdict
        Maps the names of reads aligned to this contig to another defaultdict,
        which in turn maps all (zero-indexed) positions in "positions" spanned
        by this read to the nucleotide aligned to this position at this read.
        Nucleotides are encoded as integers using config.N2I.

        Note that we only consider a read to "span" a position if it has a
        non-degenerate nucleotide (i.e. just one of A, C, G, T) aligned to this
        position using a match/mismatch operation.

        If a read does not span any of the positions of interest, then it will
        not be explicitly included as an outer key in readname2pos2nt (although
        since readname2pos2nt is a defaultdict, trying to access such a read in
        it will just give you {}).

    Raises
    ------
    WeirdError
        If things go wrong during iteration.
    """
    # We have already guaranteed (earlier in run_nt()) that this contig is
    # present in the BAM, so checking for that here would probs be overkill.

    # Part 1: build up readname2pos2nt (exactly what it says on the tin)
    readname2pos2nt = defaultdict(dict)

    # Similar song and dance to smooth_utils.get_smooth_aln_replacements(),
    # but we set matches_only=True (we don't care about deletions seen at a
    # position of interest). It might be worth abstracting this shared code to
    # a utility function, but probs not worth the trouble right now.
    for aln in bam_obj.fetch(contig):
        ap = aln.get_aligned_pairs(matches_only=True)

        # Iterating through the aligned pairs is expensive. Since read
        # lengths are generally in the thousands to tens of thousands of
        # bp (which is much less than the > 1 million bp length of most
        # bacterial genomes), we set things up so that we only iterate
        # through the aligned pairs once. We maintain an integer, pi,
        # that is a poor man's "pointer" to an index in positions.

        pi = 0

        # Go through this aln's aligned pairs. As we see each pair, compare
        # the pair's reference position (refpos) to the pi-th position (herein
        # referred to as "pos").
        #
        # If refpos >  pos, increment pi until refpos <= pos
        #                      (stopping as early as possible).
        # If refpos == pos, we have a match! Update readname2pos2nt.
        # If refpos <  pos, continue to the next pair.

        readname = aln.query_name
        for pair in ap:

            readpos, refpos = pair
            pos = positions[pi]

            no_positions_of_interest_to_right_of_here = False

            # Increment pi until we get to the next position at or after
            # after the reference position for this aligned pair (or until we
            # run out of positions in "positions").
            while refpos > pos:
                pi += 1
                if pi < len(positions):
                    pos = positions[pi]
                else:
                    no_positions_of_interest_to_right_of_here = True
                    break

            # I expect this should happen only for reads aligned near the
            # right end of the genome.
            if no_positions_of_interest_to_right_of_here:
                break

            # If the next position of interest occurs after this aligned pair,
            # continue on to a later pair.
            if refpos < pos:
                continue

            # If we've made it here, refpos == pos!
            # (...unless I messed something up in how I wrote this code.)
            if refpos != pos:
                raise WeirdError(
                    f"refpos = {refpos:,}, but pos = {pos:,}. "
                    "refpos and pos should match. This "
                    "should never happen; please, open an issue on GitHub "
                    "so you can yell at Marcus."
                )

            # Perform a quick sanity check that this read has not already had a
            # match to this position. Should never happen in practice, unless
            # the input BAM includes OSAs / secondary alns. We could make this
            # check much more comprehensive (checking ALL positions in all
            # contigs) but we just limit it here to the positions of interest,
            # because we're already assuming the input alignment doesn't have
            # this problem.
            if pos in readname2pos2nt[readname]:
                raise ParameterError(
                    f"Read {readname} is aligned to the 0-indexed position "
                    f"{pos:,} in contig {contig} multiple times. Make sure "
                    "that you have removed secondary alignments and "
                    "overlapping supplementary alignments from your alignment "
                    'file; you can use "strainFlye align" to compute such an '
                    "alignment."
                )

            # If this read has a degenerate nucleotide aligned to a position of
            # interest, it doesn't count (for our purposes at the moment, at
            # least).
            read_nt = aln.query_sequence[readpos].upper()
            if read_nt not in "ACGT":
                continue

            # Convert the nucleotide at this position on this read to an
            # integer in the range [0, 3] using N2I
            readval = config.N2I[read_nt]

            # Record this specific "allele" for this read.
            readname2pos2nt[readname][pos] = readval

    return readname2pos2nt


def get_pos_nt_info(readname2pos2nt):
    """Converts read-level nucleotide info to nucleotide (co-)occurrence info.

    Parameters
    ----------
    readname2pos2nt: defaultdict
        Maps read name -> (zero-indexed) position -> nucleotide (encoded as an
        integer using config.N2I). We only focus on the "positions of interest"
        described in this dict, which in practice will probably mean mutated
        positions. Can be computed by get_readname2pos2nt().

    Returns
    -------
    pos2nt2ct, pospair2ntpair2ct: defaultdict, defaultdict
        pos2nt2ct maps (one-indexed) position of interest -> nucleotide seen at
        this position (still encoded as an integer) -> count of times we saw
        this nucleotide at this position across all of the reads in
        readname2pos2nt. This corresponds to reads(i, N) as described in the
        strainFlye paper.

        pospair2ntpair2ct maps (sorted so the leftmost position comes first)
        pairs of (one-indexed) positions of interest -> pairs of nucleotides
        seen at these two positions, such that the left nt is for the left
        position in this pair (still encoding nts as integers) --> count of
        times we saw this pair of nucleotides at these positions across all of
        the reads in readname2pos2nt. Corresponds to reads(i, j, Ni, Nj) as
        described in the strainFlye paper.

        Note that the inclusion of a pair of possitions of interest in the
        outer layer of pospair2ntpair2ct implies that these two positions were
        spanned by at least one read. Also, there are 16 possible pairs of
        nucleotides for any pair of positions, since 4^2 = 16 (and we should
        have implicitly ignored degen nucleotides, deletions, etc. during
        get_readname2pos2nt()).

    Raises
    ------
    WeirdError
        If things go wrong with the output of itertools.combinations().

    Notes
    -----
    This converts positions from zero-indexing to one-indexing. Probs would've
    been better to do one-indexing from the start in get_readname2pos2nt(), but
    I'm not sure that'd be worth the trouble.

    The main motivation for converting to one-indexing is that the output of
    "strainFlye link graph" will use one-indexing also, so we may as well be
    consistent with the output of "strainFlye link nt".

    Also: if there is only one "position of interest" included, then
    pospair2ntpair2ct will be empty.

    Example
    -------
    In this example, there are two reads. These reads' "haplotypes" (using
    one-indexed positions) are:

    r1: 1 -> T, 6 -> G
    r2: 1 -> C, 5 -> A

    >>> pos2nt2ct, pospair2ntpair2ct = get_pos_nt_info(
    ... {"r1": {0: 3, 5: 2}, "r2": {0: 1, 4: 0}}
    ... )
    >>> assert pos2nt2ct == {
    ... 1: {1: 1, 3: 1},
    ... 5: {0: 1},
    ... 6: {2: 1},
    ... }
    >>> assert pospair2ntpair2ct == {
    ... (1, 6): {(3, 2): 1},
    ... (1, 5): {(1, 0): 1},
    ... }

    (Minor note about the above doctests: the ...s are needed because otherwise
    doctest complains -- see https://stackoverflow.com/a/46083961.)
    """
    pos2nt2ct = defaultdict(gen_ddi)
    pospair2ntpair2ct = defaultdict(gen_ddi)
    # NOTE: If a read did not span any positions of interest in this contig,
    # then it will not be included in this iteration. This is as expected --
    # such a read has nothing to "contribute" to the information we are
    # computing here.
    for readname in readname2pos2nt:
        # TODO: see if we can avoid sorting here: inefficient
        # when done once for every read, maybe?
        positions_of_interest_covered_in_read = sorted(
            readname2pos2nt[readname].keys()
        )

        # NOTE: it may be possible to include this in the combinations()
        # loop below, but we'd need some snazzy logic to prevent updating
        # the same position multiple times. Also, this would cause problems if
        # there is exactly one position of interest covered in a read. Easiest
        # for my sanity to just be a bit inefficient and make this two separate
        # loops.
        for pos in positions_of_interest_covered_in_read:
            # Convert to one-indexing
            pos2nt2ct[pos + 1][readname2pos2nt[readname][pos]] += 1

        for (i, j) in combinations(positions_of_interest_covered_in_read, 2):

            # We can assume that i and j are sorted because
            # positions_of_interest_covered_in_read is sorted: see
            # https://docs.python.org/3.10/library/itertools.html#itertools.combinations
            # This is guaranteed, but let's be paranoid just in case:
            if j <= i:
                raise WeirdError("combinations() isn't preserving order?")

            # these are integers in the range [0, 3] thanks to config.N2I
            i_nt = readname2pos2nt[readname][i]
            j_nt = readname2pos2nt[readname][j]

            # We know these positions of interest were observed on the same
            # read, and we know the exact nucleotides this read had at both
            # positions -- update this in pospair2ntpair2ct
            #
            # Also, convert to one-indexing for these output files.
            pospair2ntpair2ct[(i + 1, j + 1)][(i_nt, j_nt)] += 1

    return pos2nt2ct, pospair2ntpair2ct


def write_obj_to_pickle(obj, output_dir, contig_name, obj_name):
    """Writes out an object to a pickle file.

    This file will be named [contig_name]_[obj_name].pickle, and will be
    located in the directory [output_dir].

    Parameters
    ----------
    obj: object
        Something to write out to a pickle file.

    output_dir: str
        Directory to which we'll write this pickle file. Should already exist.

    contig_name: str
        Used as the first part of the filename.

    obj_name: str
        Used as the second part of the filename.

    Returns
    -------
    None
    """
    fp = os.path.join(output_dir, f"{contig_name}_{obj_name}.pickle")
    with open(fp, "wb") as dumpster:
        pickle.dump(obj, dumpster)


def run_nt(contigs, bam, bcf, output_dir, verbose, fancylog):
    """Computes (co-)occurrence information for nucleotides at mutations.

    Basically, this is split up into its own command to make it easier for the
    user to construct a graph based on this information multiple times, if
    desired. Parameters meaningful to graph construction (apart, of course,
    from the selection of which positions are mutated or not) are mostly
    contained in the step after this ("strainFlye link graph") -- this
    "strainFlye link nt" step is boring, straightforward, and relatively
    time-consuming, so it makes sense to separate it.

    For each contig, writes out two "pickle" files to the output directory: one
    file named [contig]_[config.POS_FILE_LBL].pickle, and one file named
    [contig]_[config.POSPAIR_FILE_LBL].pickle.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a BAM file mapping reads to contigs.

    bcf: str
        Filepath to a BCF file describing single-nucleotide mutations.

    output_dir: str
        Directory to which we'll write out information for each contig.

    verbose: bool
        Log extra info.

    fancylog: function
        Logging function.

    Returns
    -------
    None

    Raises
    ------
    Various errors are raised by misc_utils.load_triplet() if the input files
    are problematic.
    """
    verboselog = cli_utils.get_verboselog(fancylog, verbose)
    # Load and check the FASTA, BAM, and BCF files.
    #
    # misc_utils.load_triplet() will verify that:
    # - all contigs in the FASTA are present in the BAM and BCF files
    # - there is at least one contig in the FASTA
    # - the BAM and BCF files seem well-formed (e.g. the BCF doesn't describe
    #   indels)
    # This is in addition to writing out some logging messages.
    #
    # (tt and tm will both be None, since we don't care about whether or not
    # this BCF came from strainFlye)
    contig_name2len, bam_obj, bcf_obj, tt, tm = misc_utils.load_triplet(
        contigs, bam, bcf, fancylog
    )
    num_contigs = len(contig_name2len)

    misc_utils.make_output_dir(output_dir)

    fancylog(
        "Going through contigs and computing nucleotide (co-)occurrence "
        "information..."
    )
    for ci, contig in enumerate(contig_name2len, 1):
        clen = contig_name2len[contig]
        cli_utils.proglog(contig, ci, num_contigs, verboselog, contig_len=clen)

        # note that mutated positions here are zero-indexed; this makes it
        # easier to compare stuff with pysam fetch() output below, which also
        # uses zero-indexing.
        mutated_positions = sorted(
            bcf_utils.get_mutated_positions_in_contig(bcf_obj, contig)
        )
        if len(mutated_positions) == 0:
            verboselog(
                f"Contig {contig} has no mutations; ignoring it.", prefix=""
            )
            continue
        else:
            verboselog(
                (
                    f"Contig {contig} has {len(mutated_positions):,} mutated "
                    "position(s). Going through linear alignments to it..."
                ),
                prefix="",
            )

        readname2mutpos2nt = get_readname2pos2nt(
            bam_obj, contig, mutated_positions
        )

        verboselog(
            (
                "Done with that; now computing (co-)occurrence information "
                f"for contig {contig}..."
            ),
            prefix="",
        )

        pos2nt2ct, pospair2ntpair2ct = get_pos_nt_info(readname2mutpos2nt)

        write_obj_to_pickle(pos2nt2ct, output_dir, contig, config.POS_FILE_LBL)
        write_obj_to_pickle(
            pospair2ntpair2ct, output_dir, contig, config.POSPAIR_FILE_LBL
        )

        verboselog(
            (
                "Done with that; wrote out (co-)occurrence information for "
                f"contig {contig}."
            ),
            prefix="",
        )
    fancylog("Done.", prefix="")


def write_linkgraph_to_dot(g, output_dir, contig_name):
    """Writes out a NetworkX link graph to a DOT file.

    This file will be named [contig_name]_linkgraph.gv, and will be
    located in the directory [output_dir].

    Parameters
    ----------
    g: nx.Graph
        Link graph to write out. This should actually be a link graph (i.e.
        containing the extra node/edge attributes like "link" for edges), not
        just an arbitrary NetworkX graph.

    output_dir: str
        Directory to which we'll write this DOT file. Should already exist.

    contig_name: str
        Used in the name of the DOT file to be created.

    Returns
    -------
    None
    """
    fp = os.path.join(output_dir, f"{contig_name}_linkgraph.gv")
    with open(fp, "w") as df:
        # We could name this graph according to the contig, or something, but I
        # don't wanna have to worry about Graphviz breaking due to contigs
        # containing weird characters (using these as part of filenames is
        # enough of a headache already). So let's just omit a name for now.
        # (This is consistent with DOT files from LJA at the moment, at least
        # as far as I can tell)
        df.write("graph {\n")

        # We represent node names in the DOT file by just putting their "real"
        # name (from networkx) in double-quotes. Graphviz *can* support things
        # like integers as node names, but just treating all node names as
        # strings from the start of developing this should save us some
        # headaches. (Knock on wood.)
        for n in g.nodes:
            node_lbl = (
                f"{n[0]:,} ({config.I2N[n[1]]})\\n"
                f"{g.nodes[n]['ct']:,}x "
                f"({(100 * g.nodes[n]['freq']):.2f}%)"
            )
            # we use two spaces of indentation for node lines and edge lines,
            # matching https://www.graphviz.org/doc/info/lang.html, although
            # this doesn't actually matter and the time that i spent writing
            # this comment is probably now wasted to the sands of time wait why
            # am i still writing this on aodsfih dofhid ofuhasd foudhf odufhod
            df.write(f'  "{n}" [label="{node_lbl}"];\n')

        for e in g.edges:
            # In theory, an edge's "link" value could be anywhere in the range
            # (0, 1] (the requirement that --low-link is at least 0 forces link
            # to be positive, if we have created an edge in the first place).
            # So doing interpolation to scale these values into
            # (0, MAX_PENWIDTH] isn't too challenging.
            #
            # However, since a very tiny link value implies a very tiny
            # penwidth, we clamp the min penwidth.
            lw = g.edges[e]["link"] * config.MAX_PENWIDTH
            if lw < config.MIN_PENWIDTH:
                lw = config.MIN_PENWIDTH
            df.write(f'  "{e[0]}" -- "{e[1]}" ' f"[penwidth={lw}];\n")

        df.write("}")


def run_graph(
    nt_dir,
    min_nt_ct,
    min_span,
    low_link,
    output_format,
    output_dir,
    verbose,
    fancylog,
):
    verboselog = cli_utils.get_verboselog(fancylog, verbose)
    havent_created_first_graph_yet = True

    pnf_suffix = f"_{config.POS_FILE_LBL}.pickle"
    pair_suffix = "f_{config.POSPAIR_FILE_LBL}.pickle"

    fancylog(
        "Going through co-occurrence information and creating link graphs..."
    )

    # Sorting os.listdir() makes this deterministic. We use this same trick in
    # "strainFlye smooth assemble", also.
    #
    # (We could also take the FASTA file of contigs as input, but it isn't
    # necessary so let's not bother)
    for fp in sorted(os.listdir(nt_dir)):
        if fp.lower().endswith(pnf_suffix):
            pnf_abs_fp = os.path.join(nt_dir, fp)
            with open(pnf_abs_fp, "rb") as loadster:
                pos2nt2ct = pickle.load(loadster)

            contig = fp[: -len(pnf_suffix)]
            corresponding_pair_abs_fp = os.path.join(
                nt_dir, contig + pair_suffix
            )
            if not os.path.exists(corresponding_pair_abs_fp):
                raise FileNotFoundError(
                    f"Found file {pnf_abs_fp}, but not file "
                    f"{corresponding_pair_abs_fp}."
                )
            with open(corresponding_pair_abs_fp, "rb") as loadster:
                pospair2ntpair2ct = pickle.load(loadster)

            # we don't have enough information about the available contigs to
            # easily use cli_utils.proglog(). We could try to get around this
            # but ehhhh easier to just hack a one-off solution here
            verboselog(
                (
                    f"Found both info files for contig {contig}; adding nodes "
                    "to its link graph..."
                ),
                prefix="",
            )

            # Alright, now create the graph
            g = nx.Graph()

            # Add nodes to the graph
            for pos in pos2nt2ct.keys():

                pos_cov = sum(pos2nt2ct[pos].values())

                # Since this data structure is a defaultdict, this will only
                # iterate over the defined (i.e. seen) nucleotide indices
                # (integers in the range [0, 3]).
                for nt in pos2nt2ct[pos].keys():
                    # Set the "ct" attribute of this allele node to the
                    # number of times this nucleotide was seen at this position
                    # in the reads. This corresponds to reads(i, Ni) for
                    # position i and nucleotide Ni.
                    ct = pos2nt2ct[pos][nt]

                    if ct >= min_nt_ct:
                        # Also, set the "freq" attribute to ct divided by
                        # the total number of matching operations at this
                        # nucleotide -- so we can see what percentage of reads
                        # at this position had a given nucleotide.
                        g.add_node((pos, nt), ct=ct, freq=(ct / pos_cov))

            verboselog(
                (
                    f"The link graph for contig {contig} has {len(g.nodes):,} "
                    "node(s). Adding edges..."
                ),
                prefix="",
            )

            for pospair in pospair2ntpair2ct:
                i = pospair[0]
                j = pospair[1]

                # NOTE: possible to speed this up by bundling this computation
                # into the for loop below, maybe also note that "num spanning
                # reads" only includes reads that meet criteria about not
                # having skips/indels at either position, etc.
                num_spanning_reads = sum(pospair2ntpair2ct[pospair].values())

                if num_spanning_reads >= min_span:
                    for ntpair in pospair2ntpair2ct[pospair]:
                        # these are still ints in the range [0, 3]
                        i_nt = ntpair[0]
                        j_nt = ntpair[1]
                        # if one or both of the nodes failed the reads(i, N)
                        # check above due to min_nt_ct, definitely don't create
                        # an edge adjacent to them!
                        if g.has_node((i, i_nt)) and g.has_node((j, j_nt)):
                            link = pospair2ntpair2ct[pospair][ntpair] / max(
                                pos2nt2ct[i][i_nt], pos2nt2ct[j][j_nt]
                            )
                            if link > low_link:
                                # Yay, add an edge between these alleles!
                                g.add_edge((i, i_nt), (j, j_nt), link=link)

            verboselog(
                (
                    f"The link graph for contig {contig} has {len(g.edges):,} "
                    "edge(s)."
                ),
                prefix="",
            )

            if havent_created_first_graph_yet:
                misc_utils.make_output_dir(output_dir)
                havent_created_first_graph_yet = False

            if output_format == "nx":
                write_obj_to_pickle(g, output_dir, contig, "linkgraph")
            elif output_format == "dot":
                write_linkgraph_to_dot(g, output_dir, contig)
            else:
                raise WeirdError("Unrecognized output format: {output_format}")

            verboselog(
                f"Wrote out the link graph for contig {contig}.",
                prefix="",
            )

    if havent_created_first_graph_yet:
        raise ParameterError(
            f"Didn't find any (co-)occurrence information in {nt_dir}."
        )

    fancylog("Done.", prefix="")
