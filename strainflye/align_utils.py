# Utilities for strainFlye's alignment step.

import re
import subprocess
import pysam
from itertools import combinations
from collections import defaultdict
from . import graph_utils


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
    fancylog(f"Done indexing the {bam_descriptor}.", prefix="")


def get_coords(alnseg):
    """Returns a linear alignment's coordinates.

    Parameters
    ----------
    alnseg: pysam.AlignedSegment

    Returns
    -------
    (s, e): (int, int, str)
        Segment start and end.

        The start and end are both inclusive, to simplify comparison of
        alignment ranges for detecting overlaps.

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
    s = alnseg.reference_start
    e = alnseg.reference_end + 1
    if s > e:
        raise ValueError(
            f"Malformed linear alignment coordinates: start {s}, end {e}"
        )
    return (s, e)


def filter_osa_reads(in_bam, out_bam, fancylog, verbose):
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

    verbose: bool
        Log extra info about individual contigs. This comes out to multiple
        lines per contig, so it may be best to set this to False if the graph
        has thousands of contigs.

    Returns
    -------
    None
    """
    fancylog(
        "Filtering reads with overlapping supplementary alignments (OSAs)..."
    )

    def verboselog(*args, **kwargs):
        if verbose:
            fancylog(*args, **kwargs)

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

    # If literally nothing is aligned to this seq (before filtering), set its
    # entry in this to True -- this way, we can skip some work later
    seq2isempty = defaultdict(bool)

    for si, seq in enumerate(bf.references, 1):
        pct = 100 * (si / bf.nreferences)
        verboselog(
            (
                f"OSA filter pass 1/2: on contig {seq} ({si:,} / "
                f"{bf.nreferences:,}) ({pct:.2f}%)."
            ),
            prefix="",
        )

        # Identify all linear alignments of each read to this sequence
        num_lin_alns = 0
        readname2Coords = defaultdict(list)
        for num_lin_alns, linearaln in enumerate(bf.fetch(seq), 1):
            rn = linearaln.query_name
            alncoords = get_coords(linearaln)
            readname2Coords[rn].append(alncoords)

        # How many (unique) reads are aligned total to this contig?
        n_reads_in_seq = len(readname2Coords)
        if n_reads_in_seq == 0:
            seq2isempty[seq] = True
            verbose(
                f"Nothing is aligned to contig {seq}! Ignoring this contig."
            )
            continue

        # Sanity checking -- should never happen (TM) because we should have
        # already continued if n_reads_in_seq == 0
        if num_lin_alns == 0:
            raise ValueError(
                "0 linear alns, but > 0 aligned reads? Something's wrong."
            )

        # Identify overlapping alignments from the same read
        n_reads_w_osa_in_seq = 0
        for rn in readname2Coords:
            alns = readname2Coords[rn]
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

        verboselog(
            f"There are {num_lin_alns:,} linear alignment(s) (from "
            f"{n_reads_in_seq:,} unique read(s)) to contig {seq}.",
            prefix="",
        )
        # We can compute this percentage without worrying about division by
        # zero because we've already ensured above that n_reads_in_seq != 0.
        rpct = 100 * (n_reads_w_osa_in_seq / n_reads_in_seq)
        verboselog(
            f"{n_reads_w_osa_in_seq:,} / {n_reads_in_seq:,} ({rpct:.2f}%) "
            "of these unique read(s) have OSAs.",
            prefix="",
        )

    fancylog(
        "Done with pass 1 of the OSA filter; moving on to pass 2...",
        prefix="",
    )

    # Now that we've made note of all reads with OSAs across *all* sequences in
    # the alignments, we can make another pass through and output all reads
    # without OSAs into a new BAM file.

    # Output BAM file (filtered to remove reads with overlapping supplementary
    # alignments, aka OSAs)
    of = pysam.AlignmentFile(out_bam, "wb", template=bf)

    # TODO: maybe generalize this iteration code into a generator or something
    # to limit code reuse
    for si, seq in enumerate(bf.references, 1):

        # Ignore already-known-to-be empty sequences
        if seq2isempty[seq]:
            continue

        pct = 100 * (si / bf.nreferences)
        verboselog(
            f"OSA filter pass 2/2: on contig {seq} ({si:,} / "
            f"{bf.nreferences:,}) ({pct:.2f}%).",
            prefix="",
        )

        num_alns_retained = 0
        num_alns_filtered = 0
        for linearaln in bf.fetch(seq):
            rn = linearaln.query_name
            # If this read has OSAs anywhere in the alignment (even if it's to
            # other contigs), don't include it in the output BAM file.
            # Otherwise, *do* include it!
            if rn in reads_with_osa:
                num_alns_filtered += 1
            else:
                of.write(linearaln)
                num_alns_retained += 1

        # Like when we computed rpct above, we know at this point that seq is
        # non-empty (at least pre-OSA-filtering), so num_alns_total must be > 0
        num_alns_total = num_alns_retained + num_alns_filtered
        apct = 100 * (num_alns_retained / num_alns_total)
        verboselog(
            f"{num_alns_retained:,} / {num_alns_total:,} ({apct:.2f}%) "
            f"linear aln(s) retained in contig {seq}.",
            prefix="",
        )

    bf.close()
    of.close()

    fancylog(
        "Done filtering reads with overlapping supplementary alignments.",
        prefix="",
    )


def filter_pm_reads(
    gfa,
    in_bam,
    out_bam,
    fancylog,
    verbose,
    min_percent_aligned=90,
    too_many_adj_contigs=50,
):
    """Filters a BAM file to remove partially-mapped reads.

    Parameters
    ----------
    gfa: str or None
        Location of a GFA file whose segments correspond to contigs in the
        assembly graph. This is used when determining whether or not a read is
        partially-mapped. If this is None, then we just don't consider adjacent
        contigs.

    in_bam: str
        Location of the BAM file on which we will perform filtering. This BAM
        file should already have been sorted and indexed, and should have had
        OSAs filtered as well.

    out_bam: str
        Location to which an output BAM file (containing no partially-mapped
        reads) will be written.

    fancylog: function
        Logging function.

    verbose: bool
        Log extra info about individual contigs.

    min_percent_aligned: float
        Number in the range [0, 100]. A given read is only included in the
        filtered alignment information for some contig if the sum of all
        (mis)match operations in this read in this contig -- and in adjacent
        contigs, if applicable -- is at least this percentage of the read
        length.

    too_many_adj_contigs: int
        If a contig C has at least this many adjacent contigs in the graph,
        then -- when checking if a read is partially aligned to C -- we will
        not bother checking the adjacent contigs of C. Why do we do this?

          - Rationale 1: this saves us time, because there are likely a lot of
            cases in the "hairball" component where a contig has tons of
            adjacencies.

          - Rationale 2: Many of these "tons of adjacencies" are likely
            not very meaningful. The main reason that we check adjacent contigs
            in the first place is because we don't want to unjustly penalize
            contigs that are present in components that happen to be
            incompletely assembled.

        Setting this to 0 or 1 implies that adjacent contigs will never be
        considered at all (this is fine). Setting this to a negative number
        disables the check, meaning that adjacent contigs will ALWAYS be
        considered. (This is not encouraged unless you really know what you're
        doing, since this could drastically increase runtime...)

    Returns
    -------
    None

    References
    ----------
    I realized after the fact that this bears some resemblance to samclip
    (https://github.com/tseemann/samclip), although this differs a bit in the
    sort of alignments this allows to pass the filter and the sort of
    information it takes into account. These are probs ultimately minor
    distinctions, tho.
    """
    fancylog("Filtering partially-mapped reads...")

    def verboselog(*args, **kwargs):
        if verbose:
            fancylog(*args, **kwargs)

    # Sanity check
    if min_percent_aligned < 0 or min_percent_aligned > 100:
        # not gonna bother trying to format this number nicely because it could
        # be, well, anything outside of [0, 100]
        raise ValueError(
            f"min_percent_aligned = {min_percent_aligned} is not in [0, 100]."
        )

    bf = pysam.AlignmentFile(in_bam, "rb")
    bam_contigs = set(bf.references)

    # Super duper ultra paranoid sanity check.
    # This should never happen, but we originally used these two numbers
    # interchangeably and i wanna verify this just in case pysam breaks
    if len(bam_contigs) != bf.nreferences:
        raise ValueError("This BAM file is cursed. Call a priest.")

    graph = None
    if gfa is not None:
        fancylog("Loading assembly graph...", prefix="")
        graph = graph_utils.load_gfa(gfa)

        # Another sanity check!
        gfa_nodes = set(graph.nodes())
        if bam_contigs != gfa_nodes:
            raise ValueError(
                f"Contigs in the BAM file ({in_bam}) and segments in the "
                f"graph ({gfa}) do not match.\nThe BAM file has "
                f"{bf.nreferences:,} contigs and the GFA has "
                f"{len(gfa_nodes):,} segments, for reference."
            )

        fancylog("Loaded assembly graph.", prefix="")
    else:
        fancylog("No assembly graph given.", prefix="")

    # Per the SAM v1 specification, top of page 8:
    # "Sum of lengths of the M/I/S/=/X operations shall equal the length of
    # SEQ."
    #
    # These are the five CIGAR operations that consume character(s) from the
    # query sequence (i.e. a read).
    #
    # We ignore S (soft clipping), since this indicates that a given position
    # in a read is not matched to a seq; and we ignore I (insertion), since
    # this also indicates that a position in a read is not really "matched"
    # anywhere on the seq (although I guess you could argue that including
    # insertions here may be useful; it probably isn't a big deal either way).
    #
    # By just looking at M, X, and = occurrences, we can get a count for each
    # alignment of the number of bases "(mis)matched" to a seq.
    # Since we have already filtered secondary alignments and overlapping
    # supplementary alignments, we can sum these match counts across all
    # alignments of a read and divide by the read length to get the approx
    # percentage (see above for slight caveats) of bases in the read aligned to
    # a contig or group of contigs.
    #
    # (By the way, use a raw string to avoid confusing python / flake8 into
    # thinking that we're trying to add an escape sequence:
    # https://stackoverflow.com/a/52335971)
    matches = re.compile(r"(\d+)[MX=]")

    def check_and_update_alignment(
        aln, readname2len, readname2matchct, contig_name
    ):
        """Updates two dicts in place based on a single linear alignment.

        This is intended to be called multiple times, as we consider lots and
        lots of alignments to a contig (or its adjacent contigs, sometimes) in
        rapid succession.

        Parameters
        ----------
        aln: pysam.AlignedSegment

        readname2len: dict
            Maps read name (str) to the length of this read (int). The first
            dict that will be updated. A read's length is important information
            for us to know when we later define a read as partially-mapped or
            not (since the read length is the denominator in that fraction).

            It's very possible that, when this function is called, it will have
            already been called for another alignment of a given read. In this
            case, readname2len won't be updated, but it will be checked -- if
            the read lengths from the previous and current alignments of this
            read disagree, then this function will raise an error.

        readname2matchct: defaultdict(int)
            Maps read name (str) to the number of (mis)match operations of this
            read observed on this contig (int). Note that this number may
            include operations to adjacent contigs (...if this function is
            called across multiple contigs).

        contig_name: str
            Name of the contig to which aln is aligned. This parameter is just
            used here for sanity checking.

        Raises
        ------
        ValueError
            - If the inferred lengths of a read are observed to be inconsistent
              (see readname2len's description).

            - If contig_name and aln.reference_name do not match.

            - If aln's CIGAR string has no (mis)match operations.
        """
        # Ensure that read length is consistent across all alignments involving
        # this read; also, make a record of previously unseen reads' lengths.
        # (Apparently query length is dependent on the actual alignment of this
        # read, so we use infer_read_length() instead because we care about the
        # actual length of the read.)
        obs_read_len = aln.infer_read_length()
        if aln.query_name in readname2len:
            if readname2len[aln.query_name] != obs_read_len:
                raise ValueError(
                    "Inconsistent read lengths across alignments of read "
                    f"{aln.query_name}: prev aln = "
                    f"{readname2len[aln.query_name]:,}, this aln = "
                    f"{obs_read_len:,}"
                )
        else:
            readname2len[aln.query_name] = obs_read_len

        # Each AlignedSegment returned by fetch(contig) should pertain to that
        # specific contig
        if aln.reference_name != contig_name:
            raise ValueError(
                f"Alignment reference name, {aln.reference_name}, isn't "
                f"{contig_name} as expected"
            )

        # The meat of this -- parse the CIGAR string of this alignment and
        # count all (mis)match operations, updating a defaultdict.
        allmatches = matches.findall(aln.cigarstring)
        if allmatches:
            matchct = sum([int(c) for c in allmatches])
            readname2matchct[aln.query_name] += matchct
        else:
            # Raise an error if an alignment of this read does not involve
            # any (mis)match operations at all. This *could* happen in
            # practice, I guess, but if it does something is likely wrong. If
            # this check needs to be removed in the future, then this block
            # could just be commented out or replaced with a "pass" statement
            # or something.
            raise ValueError(
                "No (mis)match operations (M/X/=) found in an alignment of "
                f"read {aln.query_name}.\nThe CIGAR string for this alignment "
                f"is {aln.cigarstring}, for reference.\nIf you encountered "
                "this on a real dataset, please let the developers know and "
                "we can look into removing this check."
            )

        # The function is done now -- we've updated the two dicts based on
        # this alignment, and all sanity checks have passed.

    of = pysam.AlignmentFile(out_bam, "wb", template=bf)

    # Figure out all reads that are aligned to each contig to focus on
    for ci, contig_to_focus_on in enumerate(bf.references, 1):
        # just for convenience's sake, since we write this out a lot
        cdsc = f"contig {contig_to_focus_on}"
        pct = 100 * (ci / bf.nreferences)
        verboselog(
            f"PM read filter: on {cdsc} ({ci:,} / "
            f"{bf.nreferences:,}) ({pct:.2f}%).",
            prefix="",
        )
        # Maps read name to read length (which should be constant across all
        # alignments of that read). This variable is used both to store this
        # info (which is in turn used for sanity checking) as well as as a
        # crude indication of "have we seen this read yet?"
        readname2len = {}

        # Maps read name to number of match operations to the sequences of this
        # contig or adjacent contigs in the graph.
        readname2matchct = defaultdict(int)

        num_lin_alns = 0
        for num_lin_alns, aln in enumerate(bf.fetch(contig_to_focus_on), 1):
            check_and_update_alignment(
                aln, readname2len, readname2matchct, contig_to_focus_on
            )

        # Analogous to the seq2isempty check in the OSA filter function.
        # Here, we don't need to define such a data structure, because we only
        # apply a single "pass" over our BAM file -- so we can just move on
        # immediately. (Since no reads are aligned to this contig, there are no
        # linear alignments to this contig, and thus nothing from this contig
        # could possibly be included in the output BAM file.)
        if num_lin_alns == 0:
            verboselog(
                f"Nothing is aligned to {cdsc}! Ignoring this contig.",
                prefix="",
            )
            continue
        else:
            verboselog(
                f"{num_lin_alns:,} linear alignment(s) to {cdsc}.", prefix=""
            )

        if graph is not None:
            # Identify adjacent contigs to this one in the graph, if present.
            # Allow alignments to these contigs to count towards
            # readname2matchct (so we can somewhat avoid penalizing contigs
            # that aren't in isolated components). We could also expand this to
            # include _all_ other contigs in this contig's component, if
            # desired (although due to hairballs that will likely cause
            # problems unless we still impose the bound on nae).
            #
            # The removal of set([contig_to_focus_on]) is done because, if a
            # contig has an edge to itself, then it'll be considered by
            # NetworkX as a neighbor of itself (and thus will be included
            # in adj_contigs). We don't want this, which is why we remove it!
            adj_contigs = set(graph.neighbors(contig_to_focus_on)) - set(
                [contig_to_focus_on]
            )
            nae = len(adj_contigs)

            verboselog(
                f"{nae:,} contig(s) in the graph are adjacent to {cdsc}.",
                prefix="",
            )

            # To prevent this script from taking a super long amount of time,
            # only allow alignments to adjacent contigs if there are less than
            # "too_many_adj_contigs" adjacent contigs to this contig in the
            # graph.
            if nae > 0:
                consider_adj = True

                # too_many_adj_contigs being negative implies that we should
                # always consider adj contigs.
                if too_many_adj_contigs >= 0 and nae >= too_many_adj_contigs:
                    verboselog(
                        "Too many adjacent contigs; we won't consider them "
                        "when filtering partially-mapped reads from this "
                        "contig.",
                        prefix="",
                    )
                    consider_adj = False

                if consider_adj:
                    verboselog(
                        "Considering alignments of shared reads to these "
                        "adjacent contig(s)...",
                        prefix="",
                    )

                    # Go through these "allowed" other contigs; add to the
                    # number of matches in cc for any reads that we see that
                    # we've already seen in the contig to focus on. (We
                    # implicitly ignore any reads aligned to these other
                    # contigs but not to the contig to focus on.)
                    num_other_contig_alns_from_shared_reads = 0
                    for other_contig in adj_contigs:
                        for aln in bf.fetch(other_contig):
                            # If aln.query_name is NOT in readname2len, then
                            # this alignment's read wasn't also aligned to the
                            # contig to focus on -- in this case we implicitly
                            # ignore it, as mentioned above, since we only
                            # really care about reads aligned to the contig
                            # we're focusing on.
                            if aln.query_name in readname2len:
                                check_and_update_alignment(
                                    aln,
                                    readname2len,
                                    readname2matchct,
                                    other_contig,
                                )
                                num_other_contig_alns_from_shared_reads += 1

                    verboselog(
                        f"{num_other_contig_alns_from_shared_reads:,} linear "
                        "alignments(s) from shared reads to adjacent "
                        f"contig(s) of {cdsc}.",
                        prefix="",
                    )

        # Now that we've considered all relevant contigs (this contig and its
        # adjacent ones, if applicable and if certain checks passed), we can
        # compute the approximate percentages of each read (aligned to the
        # contig we're focusing on) aligned to all contigs in this component.
        #
        # And, finally, we can selectively write these reads' alignments
        # to an output BAM file accordingly.
        #
        # Note that although we're iterating over all alignments, these
        # computations are identical for alignments from the same read (at
        # least in the context of the same contig_to_focus_on). So if a read
        # passes the filter for a contig, all its alignments to this contig
        # will be output to the BAM file one by one; and if the read fails the
        # filter for a contig, none of its alignments to this contig will be
        # output to the BAM file (although this does not preclude other
        # alignments from this read to other contigs from being output in the
        # context of other contigs in this function).
        passing_aln_ct = 0
        for aln in bf.fetch(contig_to_focus_on):
            if aln.query_name not in readname2len:
                raise ValueError(
                    f"We should have seen read {aln.query_name} earlier!"
                )
            read_len = readname2len[aln.query_name]

            if readname2matchct[aln.query_name] == 0:
                # As with a similar error case in check_and_update_alignment(),
                # I *guess* this could happen in practice but it is almost
                # certainly indicative of an error in most cases.
                #
                # Our earlier error check should prevent us from ever getting
                # to this point, but you never know.
                raise ValueError(
                    f"Read {aln.query_name} had no match operations done? Sus."
                )

            read_matchct_in_cc = readname2matchct[aln.query_name]
            # We say a read passes the filter if
            # (match ct / read len) >= (min_percent_aligned / 100).
            #
            # Or, equivalently, if
            # (100 * match ct) >= (min_percent_aligned * read len).
            #
            # This way we can avoid division. This is analogous to how naive
            # p-mutation calling works in the SheepGut github repo -- see
            # https://github.com/fedarko/sheepgut/blob/ee74a550145c5e0c2fb56b8ed0484b7df6decc1b/notebooks/pileup.py#L301-L314
            lhs = 100 * read_matchct_in_cc
            rhs = min_percent_aligned * read_len

            if lhs >= rhs:
                of.write(aln)
                passing_aln_ct += 1

        # If we've made it this far, we know num_lin_alns != 0. So no need to
        # worry about division by zero.
        passing_pct = 100 * (passing_aln_ct / num_lin_alns)
        verboselog(
            f"{passing_aln_ct:,} / {num_lin_alns:,} ({passing_pct:.2f}%) of "
            f"alignments to {cdsc} passed the filter.",
            prefix="",
        )

    bf.close()
    of.close()

    fancylog("Done filtering out partially-mapped reads.", prefix="")
