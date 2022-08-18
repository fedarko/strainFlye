# Utilities for strainFlye smooth.

import os
from collections import defaultdict
from strainflye import (
    phasing_utils,
    cli_utils,
    bcf_utils,
    misc_utils,
    fasta_utils,
)
from strainflye.errors import ParameterError, WeirdError


def append_reads(fasta_fp, readname2seq):
    """Appends reads to the end of a FASTA file.

    Parameters
    ----------
    fasta_fp: str
        Filepath to a FASTA file. OK if it does or doesn't exist already.

    readname2seq: dict
        Maps read name to sequence. Both read name and sequence should be
        strings.

    Returns
    -------
    None
    """
    with open(fasta_fp, "a") as of:
        for readname in readname2seq:
            # Write out both the header and sequence for each read
            of.write(f">{readname}\n{readname2seq[readname]}\n")


def get_smooth_aln_replacements(aln, mutated_positions, mp2ra):
    """Returns changes needed to make a smoothed read from a linear alignment.

    Parameters
    ----------
    aln: pysam.AlignedSegment
        Object representing a linear alignment to a contig.

    mutated_positions: list of int
        Sorted list of all (zero-indexed) mutated positions in the contig to
        which aln is aligned.

    mp2ra: dict
        Maps (zero-indexed) mutated positions in the contig (the same contig
        that aln is aligned to) to a tuple of (ref nt, alt nt), as listed in
        the BCF file. Can be generated using
        bcf_utils.get_mutated_position_details_in_contig().

    Returns
    -------
    replacements_to_make: dict or None
        If we need to "ignore" aln for some reason (it has a skip aligned to a
        mutated position, or it has a non-ref or non-alt nucleotide aligned to
        a mutated position), then this will be None to indicate that no
        smoothed read should be generated from aln.

        Otherwise, this will be a dict mapping each zero-indexed mutated
        positions in the contig (that is spanned by aln) to the nucleotides
        that aln has aligned to each position. These represent the replacements
        to make to the "reference" contig sequence to construct a smoothed read
        based on aln.
    """

    replacements_to_make = {}

    # We may choose to ignore this linear alignment, if we think it is
    # error-prone or otherwise not useful. If this gets set to True in
    # the loop below, then we'll notice this and ignore this alignment.
    ignoring_this_aln = False

    # Notably, include skips -- this way, we can figure out if the aln
    # has a deletion at a mutated position, and if so ignore this aln
    ap = aln.get_aligned_pairs(matches_only=False)

    # Iterating through the aligned pairs is expensive. Since HiFi read
    # lengths are generally in the thousands to tens of thousands of
    # bp (which is much less than the > 1 million bp length of any
    # bacterial genome), we set things up so that we only iterate
    # through the aligned pairs once. We maintain an integer, mpi,
    # that is a poor man's "pointer" to an index in mutated_positions.

    mpi = 0

    # Go through this aln's aligned pairs. As we see each pair, compare
    # the pair's reference position (refpos) to the mpi-th mutated
    # position (herein referred to as "mutpos").
    #
    # If refpos >  mutpos, increment mpi until refpos <= mutpos
    #                      (stopping as early as possible).
    # If refpos == mutpos, we have a match! Update
    #                      readname2mutpos2ismutated[mutpos] based on
    #                      comparing the read to the reference at the
    #                      aligned positions.
    # If refpos <  mutpos, continue to the next pair.

    for pair in ap:
        readpos, refpos = pair

        # Since we set matches_only (for get_aligned_pairs()) to False,
        # the alignment could include insertions (which are encoded as
        # the reference pos being set to None). We inherently ignore
        # these insertions as part of the read smoothing process.
        if refpos is None:
            continue

        mutpos = mutated_positions[mpi]

        no_mutations_to_right_of_here = False

        # Increment mpi until we get to the next mutated position at or
        # after the reference pos for this aligned pair (or until we
        # run out of mutated positions).
        while refpos > mutpos:
            mpi += 1
            if mpi < len(mutated_positions):
                mutpos = mutated_positions[mpi]
            else:
                no_mutations_to_right_of_here = True
                break

        # I expect this should happen only for reads aligned near the
        # right end of the genome.
        if no_mutations_to_right_of_here:
            break

        # If the next mutation occurs after this aligned pair, continue
        # on to a later pair.
        if refpos < mutpos:
            continue

        # If we've made it here, refpos == mutpos!
        # (...unless I messed something up in how I designed this.)
        if refpos != mutpos:
            raise WeirdError(
                f"refpos = {refpos}, but mutpos = {mutpos}. "
                "refpos and mutpos should match. This "
                "should never happen; please, open an issue on GitHub "
                "so you can yell at Marcus."
            )

        # Since we set matches_only (for get_aligned_pairs()) to False,
        # there's a chance a read contains a deletion aligned to a
        # mutated position. This is accounted for by this check: in
        # this case, we ignore this alignment.
        if readpos is None:
            ignoring_this_aln = True
            break

        read_nt = aln.query_sequence[readpos]
        # If this read doesn't match the "reference" or "alternate"
        # nucleotide at this position, ignore it.
        #
        # For mutations produced by strainFlye call, "reference" should
        # always be the consensus nucleotide and "alternate" should
        # always be the second-most-common nucleotide at a position.
        #
        # For arbitrary mutations, this isn't necessarily true; we
        # leave things up to the discretion of however these mutations
        # were called.
        if read_nt != mp2ra[mutpos][0] and read_nt != mp2ra[mutpos][1]:
            ignoring_this_aln = True
            break

        # When constructing the smoothed read for this linear alignment,
        # to use the linear alignment's read's nucleotide at this mutated
        # position. (It's entirely possible that this nucleotide might match
        # the contig sequence at this position, in which case this replacement
        # won't actually change the sequence -- but this is OK.)
        relative_pos_on_aln = mutpos - aln.reference_start
        replacements_to_make[relative_pos_on_aln] = read_nt

    if ignoring_this_aln:
        return None
    else:
        return replacements_to_make


def run_apply(
    contigs,
    bam,
    bcf,
    virtual_reads,
    virtual_read_well_covered_perc,
    virtual_read_flank,
    output_dir,
    verbose,
    fancylog,
    sr_chunk_size=1000,
):
    """Generates smoothed and virtual reads.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a BAM file mapping reads to contigs.

    bcf: str
        Filepath to a BCF file describing single-nucleotide mutations.

    virtual_reads: bool
        If True, create virtual reads; otherwise, don't.

    virtual_read_well_covered_perc: float
        Only used if virtual_reads is True. Used to define whether or not a
        position in a contig is "low-coverage," and should thus be spanned by
        virtual reads.

    virtual_read_flank: int
        Only used if virtual_reads is True. A virtual read spanning a
        low-coverage region of length L will be extended by virtual_read_flank
        positions on the left and right sides (clamping to the end of the
        contig if needed).

    output_dir: str
        Directory to which we'll write out reads for each contig.

    verbose: bool
        Log extra info.

    fancylog: function
        Logging function.

    sr_chunk_size: int
        After constructing this many smoothed reads, we'll write out their
        sequences. (Like in graph_utils.gfa_to_fasta(), there's a tradeoff here
        between storing a lot of stuff in memory vs. making a lot of slow write
        operations.)

    Returns
    -------
    None

    Raises
    ------
    Various types of errors if the inputs/parameters are messed up.
    (... Sorry, I'll try to make this description better later.)
    """
    contig_name2len, bam_obj, bcf_obj = phasing_utils.load_triplet(
        contigs, bam, bcf, fancylog
    )
    num_fasta_contigs = len(contig_name2len)
    if virtual_reads:
        rt = "smoothed and virtual reads"
    else:
        rt = "smoothed reads"

    misc_utils.make_output_dir(output_dir)

    fancylog(f"Going through contigs and constructing {rt}...")
    for ci, contig in enumerate(contig_name2len, 1):
        contig_len = contig_name2len[contig]
        if verbose:
            cli_utils.proglog(
                contig, ci, num_fasta_contigs, fancylog, contig_len=contig_len
            )

        mp2ra = bcf_utils.get_mutated_position_details_in_contig(
            bcf_obj, contig
        )

        if len(mp2ra) == 0:
            fancylog(
                f"Contig {contig} has no mutations; ignoring it.", prefix=""
            )
            continue

        pos2srcov = None
        if virtual_reads:
            # Map zero-indexed positions to coverage by smoothed reads
            # (This is only used if we're adding virtual reads)
            pos2srcov = [0] * contig_len

        out_reads_fp = os.path.join(output_dir, f"{contig}.fasta")
        if os.path.exists(out_reads_fp):
            raise FileExistsError(f"File {out_reads_fp} already exists.")
        # (These are zero-indexed.)
        mutated_positions = sorted(mp2ra.keys())
        if verbose:
            fancylog(
                (
                    f"This contig has {len(mutated_positions):,} mutated "
                    "position(s)."
                ),
                prefix="",
            )

        # We'll store the full sequence of the contig here. I imagine that many
        # contigs in most datasets (like in SheepGut) will have zero
        # alignments; so, we defer loading the sequence until we actually know
        # that this contig actually has alignments.
        contig_seq = None

        ######################################################################
        # Task 1: iterate through all alignments to this contig and create
        # smoothed reads.
        ######################################################################
        # Instead of just writing out every smoothed read as soon as we
        # generate it, we build up a "buffer" of these reads and then write a
        # bunch out at once. This way we limit slowdown due to constantly
        # having to open/close files.
        #
        # This object maps the name of a smoothed read to its sequence. The
        # name of a smoothed read includes the name of the original
        # (un-smoothed) read from that was used to create it; it's useful to
        # preserve this information, so we know where smoothed reads came from.
        sr_buffer = {}

        # Just for logging
        num_ignored_alns = 0
        num_sr_generated = 0
        num_alns_total = 0

        # Smoothed read names also include a number, after the original read
        # name. This is because one original read can result in multiple linear
        # alignments to a contig (e.g. a read with supplementary alignments
        # that spans the left and right sides of a circular contig), and each
        # of these linear alignments is considered separately (and could end
        # up being used to create a smoothed read). So, the number
        # distinguishes these linear alignments.
        readname2freq_so_far = defaultdict(int)

        # Go through all linear alignments of each read to this contig...
        for ai, aln in enumerate(bam_obj.fetch(contig), 1):

            readname = aln.query_name

            if aln.is_secondary:
                raise ParameterError(
                    f"Found a secondary alignment to contig {contig} (from a "
                    f"read named {readname}). "
                    "The BAM file should not contain secondary alignments."
                )

            num_alns_total += 1

            # OK, we know this contig actually has (non-secondary) alignments,
            # so load its sequence if we haven't already.
            #
            # NOTE: get_single_seq() is kind of inefficient -- in the
            # worst case, it has to iterate over the entire FASTA file to find
            # a contig's sequence. If this ends up being a bottleneck, we can
            # speed this process up by modifying get_name2len() to figure out
            # which lines a sequence's entry takes up in the FASTA file -- and,
            # from there, use this information to speed up get_single_seq().
            # (But let's not optimize this prematurely...)
            if contig_seq is None:
                contig_seq = fasta_utils.get_single_seq(contigs, contig)

            # Although we disallow secondary alignments, supplementary
            # alignments are ok! We implicitly handle these here.
            #
            # As mentioned above, different alignments of the same read will
            # have different new_readnames, because we're gonna be treating
            # them as distinct "reads" (yes I'm aware that this is a limitation
            # of this method; the paper mentions this). We should have
            # already filtered reference-overlapping supplementary alignments,
            # so these "reads" shouldn't intersect on the contig at least.
            readname2freq_so_far[readname] += 1
            new_readname = f"{readname}_{readname2freq_so_far[readname]}"

            # Should never happen.
            if new_readname in sr_buffer:
                raise WeirdError(
                    f"{new_readname} is already in the smoothed read buffer?"
                )

            # Figure out where on the contig this alignment "hits." These are
            # 0-indexed positions. (reference_end points to the position after
            # the actual final position, since these are designed to be
            # interoperable with Python's half-open intervals.)
            #
            # Of course, there likely will be indels within this range: we're
            # purposefully ignoring those here.
            ref_start = aln.reference_start
            ref_end = aln.reference_end - 1

            # This should never happen (TM)
            if ref_start > ref_end:
                raise WeirdError(
                    f"Ref start {ref_start:,} >= ref end {ref_end:,} for "
                    f"linear alignment {new_readname}?"
                )

            repls = get_smooth_aln_replacements(aln, mutated_positions, mp2ra)

            if repls is None:
                num_ignored_alns += 1
            else:
                # Smoothed read sequence from the start to end of the linear
                # alignment, completely skipping over indels and mutations.
                # Edit this seq so that if the alignment has (mis)matches to
                # any called mutated positions within this sequence, these
                # positions will be updated to match.
                # NOTE that this isn't a string -- it's a skbio.DNA object.
                sr_seq = contig_seq[ref_start : ref_end + 1]
                # Unfortunately, skbio.DNA.replace() doesn't seem to be able to
                # replace multiple positions in multiple ways at the same time,
                # so we've gotta loop through repls ourselves.
                for aln_pos in repls:
                    sr_seq = sr_seq.replace([aln_pos], repls[aln_pos])

                # Prepare this smoothed read to be written out to a FASTA file.
                # See comments above on sr_buffer. (Also, we've already
                # guaranteed new_readname isn't already in sr_buffer, so no
                # need to worry about accidentally overwriting something from
                # earlier.)
                sr_buffer[new_readname] = str(sr_seq)
                num_sr_generated += 1

                if ai % sr_chunk_size == 0:
                    append_reads(out_reads_fp, sr_buffer)
                    # Clear the buffer
                    sr_buffer = {}

                # If we care about coverage info due to creating virtual reads,
                # record which positions this smoothed read covers (of course,
                # the original read may not exactly "cover" these positions due
                # to indels, but the smoothed version will cover them).
                #
                # We purposefully delay performing this update until right now
                # -- this way, if we choose to "ignore" this linear alignment
                # at any point in the process above, we won't increase anything
                # in pos2srcov for this linear alignment.
                if virtual_reads:
                    for pos_in_sr in range(ref_start, ref_end + 1):
                        pos2srcov[pos_in_sr] += 1

        # We're probably going to have left over smoothed reads that we still
        # haven't written out, unless things worked out so that on the final
        # alignment we saw ai was exactly divisible by sr_chunk_size (that's
        # pretty unlikely unless you set sr_chunk_size to a low number).
        # So let's make one last dump of the buffer.
        if len(sr_buffer) > 0:
            append_reads(out_reads_fp, sr_buffer)

        if num_alns_total != num_sr_generated + num_ignored_alns:
            raise WeirdError(
                f"Contig {contig}: # total alns = {num_alns_total:,}, but "
                f"# sr = {num_sr_generated:,} and # ignored = "
                f"{num_ignored_alns:,}."
            )
        if num_sr_generated > 0:
            if verbose:
                fancylog(
                    (
                        f"From the {num_alns_total:,} linear alignment(s) to "
                        f"contig {contig}: constructed {num_sr_generated:,} "
                        f"smoothed read(s) and ignored {num_ignored_alns:,} "
                        "linear alignment(s)."
                    ),
                    prefix="",
                )
        else:
            fancylog(
                (
                    f"Ignored all linear alignments for contig {contig}: "
                    "couldn't generate any smoothed reads. Ignoring this "
                    "contig."
                ),
                prefix="",
            )
            continue

        ######################################################################
        # Task 2: iterate through runs of low-coverage positions and create
        # virtual reads.
        ######################################################################

        # TODO: Unlike the Phasing-LJA.ipynb notebook, we don't have average
        # coverages up front; I guess we could have users pass in the diversity
        # index TSV file which includes these, but that's gross and
        # precludes arbitrary BCF inputs. Compute from scratch? not sure how
        # long that'll take.
        # coverage_sum = 0

    fancylog("Done.", prefix="")


def run_assemble():
    pass
