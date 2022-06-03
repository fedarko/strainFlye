# Utilities for strainFlye's alignment step.

import os
import re
import subprocess
import pysam
from itertools import combinations
from collections import defaultdict
from . import graph_utils, fasta_utils
from .errors import ParameterError, SequencingDataError


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
            verboselog(
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

        # We already sanity-checked that the GFA and FASTA describe the exact
        # same sequences (at least, going by their names), so no need to do
        # that checking here.
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


def run(reads, contigs, graph, output_dir, fancylog, verbose):
    """Runs the entire alignment-and-filtering process.

    Parameters
    ----------
    reads: tuple
        Collection of paths to FASTA / FASTQ files describing sequencing reads.

    contigs: str
        Filepath to a FASTA file containing contigs.

    graph: str or None
        If str, should be a filepath to a GFA file describing an assembly
        graph. This is only used in the partially-mapped read filter; if this
        is None, this just won't be used there. See filter_pm_reads()' docs for
        more information.

    output_dir: str
        Output directory path. Will be created if it doesn't already exist.

    fancylog: function
        Logging function

    verbose: bool
        Whether or not to log extra info about individual contigs during
        filtering.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If "reads" isn't a tuple, and click betrayed me :(

    SequencingDataError
        If the set of segment names in the "graph" GFA file is different from
        the set of sequence names in the "contigs" FASTA file.

        If something is wrong with the "contigs" FASTA file -- see
        fasta_utils.get_name2len().
    """

    # reads will be a tuple
    if type(reads) != tuple:
        # if we just have a single string (e.g. click changes something in the
        # future about how variadic arguments work) then we can fix this, but
        # for now let's be defensive. (If you encounter this error, go yell at
        # marcus)
        raise ParameterError("Collection of input reads should be a tuple.")

    # Note that get_name2len() does some sanity checking on the FASTA
    fasta_name2len = fasta_utils.get_name2len(contigs)

    # Sanity-check that the GFA segments are identical to the FASTA contigs. If
    # not, we've got problems (in this case, it's probably easiest to just not
    # consider the GFA in the PM read filter).
    # We purposefully perform this check early so we can fail early, if needed.
    # See https://github.com/fedarko/strainFlye/issues/20
    if graph is not None:
        fancylog(
            "Sanity-checking that the assembly graph and contigs describe the "
            "same sequences..."
        )
        # TODO: Loading the entire graph topology is more work than we need to
        # do here; it'd be sufficient to just scan the segment names in the GFA
        # file into a list. Might result in a slight speedup, although we will
        # have to eventually parse this graph again anyway soooo...
        # (I'm not keeping the graph in memory once we parse it, since it'll be
        # quite a while until we get to the PM read filter.)
        graph_obj = graph_utils.load_gfa(graph)
        graph_nodes = set(graph_obj.nodes())
        fasta_nodes = set(fasta_name2len.keys())
        if graph_nodes != fasta_nodes:
            raise SequencingDataError(
                "GFA segment names don't match contig names in the FASTA."
            )
        fancylog("Everything looks good so far.", prefix="")

    # Make the output dir if it doesn't already exist
    os.makedirs(output_dir, exist_ok=True)
    first_output_bam = os.path.join(output_dir, "sorted-unfiltered.bam")

    # There isn't really a need to store the SAM file from minimap2, or the
    # unsorted BAM file from "samtools view". So we use piping.
    # SAMtools stuff based on:
    # https://davetang.org/wiki/tiki-index.php?page=SAMTools#Converting_SAM_directly_to_a_sorted_BAM_file # noqa
    # Python stuff based on https://stackoverflow.com/a/4846923 and
    # https://stackoverflow.com/a/9655939

    # There's probably a way to print stuff after each individual command in
    # the chain finishes, but I don't think that sorta granularity is super
    # necessary right now tbh

    threesteps = "minimap2 --> samtools view --> samtools sort"
    fancylog(f"Running {threesteps}...")

    # NOTE: the -ax asm20 preset is what we use in the paper, but later
    # versions of minimap2 have added in "-ax map-hifi" which is probs a better
    # option in most cases. Shouldn't make too much of a difference; for
    # simplicity's sake we just stick with asm20 here, but we could definitely
    # change this (or add the option to configure it) if desired
    #
    # TODO TODO TODO the I8g thing is gonna cause problems for other computers
    # -- see https://github.com/lh3/minimap2/blob/master/FAQ.md, #3.
    # The way to handle this is make minimap2 CLI params configurable by the
    # user (they can pass a string with params I guess?)
    minimap2_run = subprocess.Popen(
        [
            "minimap2",
            "-ax",
            "asm20",
            "--secondary=no",
            "-I8g",
            "--MD",
            contigs,
            *reads,
        ],
        stdout=subprocess.PIPE,
    )
    sam_to_bam_run = subprocess.Popen(
        ["samtools", "view", "-b", "-"],
        stdin=minimap2_run.stdout,
        stdout=subprocess.PIPE,
    )
    minimap2_run.stdout.close()
    bam_to_sorted_bam_run = subprocess.Popen(
        ["samtools", "sort", "-", "-o", first_output_bam],
        stdin=sam_to_bam_run.stdout,
    )
    sam_to_bam_run.stdout.close()
    bam_to_sorted_bam_run.communicate()

    fancylog(f"Done running {threesteps}.", prefix="")

    index_bam(first_output_bam, "sorted BAM", fancylog)

    osa_filter_bam = os.path.join(output_dir, "sorted-osa-filtered.bam")
    filter_osa_reads(first_output_bam, osa_filter_bam, fancylog, verbose)
    index_bam(osa_filter_bam, "OSA-filtered BAM", fancylog)

    # Now that we've finished the OSA filter step, we can remove its input BAM
    # -- the first one we generated -- to save space. These files are big
    # enough (e.g. upwards of 70 GB for the SheepGut dataset) that keeping all
    # three around at the same time might exceed the space limit on a user's
    # system.
    os.remove(first_output_bam)
    os.remove(first_output_bam + ".bai")

    pm_filter_bam = os.path.join(output_dir, "final.bam")
    filter_pm_reads(graph, osa_filter_bam, pm_filter_bam, fancylog, verbose)
    index_bam(pm_filter_bam, "final BAM", fancylog)

    # Similarly, we can remove the OSA-filtered (but not PM-filtered) BAM now.
    # The PM-filtered BAM represents the "final" BAM produced by the alignment
    # step, and is the one that should be used in downstream analyses.
    os.remove(osa_filter_bam)
    os.remove(osa_filter_bam + ".bai")
