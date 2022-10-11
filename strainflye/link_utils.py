# Utilities for strainFlye link.


import os
import pickle
from itertools import combinations
from collections import defaultdict
from strainflye import bcf_utils, cli_utils, misc_utils, config
from strainflye.errors import WeirdError


def gen_ddi():
    """Returns a new defaultdict(int).

    Needed because pickle can't handle lambda functions.
    """
    return defaultdict(int)


def run_nt(contigs, bam, bcf, output_dir, verbose, fancylog):
    """Computes (co-)occurrence information for nucleotides at mutations.

    Basically, this is split up into its own command to make it easier for the
    user to construct a graph based on this information multiple times, if
    desired. Parameters meaningful to graph construction (apart, of course,
    from the selection of which positions are mutated or not) are mostly
    contained in the step after this ("strainFlye link graph") -- this
    "strainFlye link nt" step is boring, straightforward, and relatively
    time-consuming, so it makes sense to separate it.

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

    # TODO add logs

    for ci, contig in enumerate(contig_name2len, 1):
        clen = contig_name2len[contig]
        cli_utils.proglog(contig, ci, num_contigs, verboselog, contig_len=clen)

        mp2ra = bcf_utils.get_mutated_position_details_in_contig(
            bcf_obj, contig
        )

        if len(mp2ra) == 0:
            verboselog(
                f"Contig {contig} has no mutations; ignoring it.", prefix=""
            )
            continue

        # (These positions are zero-indexed.)
        mutated_positions = sorted(mp2ra.keys())
        verboselog(
            (
                f"Contig {contig} has {len(mutated_positions):,} mutated "
                "position(s)."
            ),
            prefix="",
        )

        # Part 1: build up readname2mutpos2nt (exactly what it says on the tin)
        readname2mutpos2nt = defaultdict(dict)
        # Similar song and dance to smooth_utils.get_smooth_aln_replacements()
        # (TODO: abstract shared code...?)
        for aln in bam_obj.fetch(contig):
            ap = aln.get_aligned_pairs(matches_only=True)

            # Iterating through the aligned pairs is expensive. Since read
            # lengths are generally in the thousands to tens of thousands of
            # bp (which is much less than the > 1 million bp length of most
            # bacterial genomes), we set things up so that we only iterate
            # through the aligned pairs once. We maintain an integer, mpi,
            # that is a poor man's "pointer" to an index in mutated_positions.

            mpi = 0

            # Go through this aln's aligned pairs. As we see each pair, compare
            # the pair's reference position (refpos) to the mpi-th mutated
            # position (herein referred to as "mutpos").
            #
            # If refpos >  mutpos, increment mpi until refpos <= mutpos
            #                      (stopping as early as possible).
            # If refpos == mutpos, we have a match! Update readname2mutpos2nt.
            # If refpos <  mutpos, continue to the next pair.

            readname = aln.query_name
            for pair in ap:

                readpos, refpos = pair
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
                # (...unless I messed something up in how I wrote this code.)
                if refpos != mutpos:
                    raise WeirdError(
                        f"refpos = {refpos:,}, but mutpos = {mutpos:,}. "
                        "refpos and mutpos should match. This "
                        "should never happen; please, open an issue on GitHub "
                        "so you can yell at Marcus."
                    )

                # (Convert the nucleotide at this position on this read
                # to an integer in the range [0, 3] using N2I)
                readval = config.N2I[aln.query_sequence[readpos]]

                # Record this specific "allele" for this read.
                readname2mutpos2nt[readname][mutpos] = readval

        # Part 2: convert readname2mutpos2nt into the actual nucleotide
        # (co-)occurrence information we're going to output from here

        # Maps mutated position -> nucleotide seen at this position, summed
        # across all reads included here -> freq. This corresponds to
        # Reads(i, N) as described in the paper.
        pos2nt2freq = defaultdict(gen_ddi)

        # This defaultdict has two levels:
        # OUTER: Keys are sorted (in ascending order) 0-indexed pairs (tuples)
        #        of mutated positions. The inclusion of a pair of mutated
        #        positions in this defaultdict implies that these two mutated
        #        positions were spanned by at least one read. The value of each
        #        pair is another defaultdict:
        #
        # INNER: The keys of this inner defaultdict are pairs of integers, each
        #        in the range [0, 3]. These represent the 4 nucleotides
        #        (0 -> A, 1 -> C, 2 -> G, 3 -> T): the first entry represents
        #        the nucleotide seen at the first position in the pair (aka the
        #        position "earlier" in the genome), and the second entry
        #        represents the nucleotide seen at the second position in the
        #        pair (aka the position "later" in the genome). Of course, many
        #        bacterial genomes are circular, so "earlier" and "later" are
        #        kinda arbitrary. Anyway, there are 16 possible pairs in one of
        #        these defaultdicts, since there are 4^2 = 16 different
        #        possible combinations of two nucleotides (ignoring deletions,
        #        degenerate nucleotides, etc.) That said, I expect in practice
        #        only a handful of nucleotide pairs will be present for a given
        #        position pair. The value of each pair in this defaultdict is
        #        an integer representing the frequency with which this pair of
        #        nucleotides was observed on a spanning read at this pair of
        #        positions.
        #
        # So, as an example, if we only have two mutated positions in a genome
        # (at 0-indexed positions 100 and 500), and we saw:
        #
        # - 30    reads with an A at both positions
        # - 1,000 reads with an A at position 100 and a T at position 500
        # - 5     reads with a T at position 100 and an A at position 500
        # - 100   reads with a T at both positions
        # - 3     reads with a C at position 100 and a T at position 500
        # - 1     read  with a G at position 100 and a T at position 500
        #
        # ... then pospair2ntpair2freq would look like
        # {
        #     (100, 500): {
        #         {
        #             (0, 0): 30,
        #             (0, 3): 1000,
        #             (3, 0): 5,
        #             (3, 3): 100,
        #             (1, 3): 3,
        #             (2, 3): 1
        #         }
        #     }
        # }
        pospair2ntpair2freq = defaultdict(gen_ddi)
        for ri, readname in enumerate(readname2mutpos2nt, 1):
            # TODO: see if we can avoid sorting here: inefficient
            # when done once for every read, maybe?
            mutated_positions_covered_in_read = sorted(
                readname2mutpos2nt[readname].keys()
            )

            # NOTE: it may be possible to include this in the combinations()
            # loop below, but we'd need some snazzy logic to prevent updating
            # the same position multiple times. Easiest for my sanity to just
            # be a bit inefficient and make this two separate loops.
            for mutpos in mutated_positions_covered_in_read:
                pos2nt2freq[mutpos][readname2mutpos2nt[readname][mutpos]] += 1

            for (i, j) in combinations(mutated_positions_covered_in_read, 2):

                # We can assume that i and j are sorted because
                # mutated_positions_covered_in_read is sorted: see
                # https://docs.python.org/3.10/library/itertools.html#itertools.combinations
                # This is guaranteed, but let's be paranoid just in case:
                if j <= i:
                    raise WeirdError(
                        "combinations() isn't preserving order as expected?"
                    )

                # these are integers in the range [0, 3] thanks to config.N2I
                i_nt = readname2mutpos2nt[readname][i]
                j_nt = readname2mutpos2nt[readname][j]

                # We know these mutated positions were observed on the same
                # read, and we know the exact nucleotides this read had at both
                # positions -- update this in pospair2ntpair2freq
                pospair2ntpair2freq[(i, j)][(i_nt, j_nt)] += 1

        with open(
            os.path.join(output_dir, f"{contig}_pos2nt2freq.pickle"), "wb"
        ) as dumpster:
            pickle.dump(pos2nt2freq, dumpster)

        with open(
            os.path.join(output_dir, f"{contig}_pospair2ntpair2freq.pickle"),
            "wb",
        ) as dumpster:
            pickle.dump(pospair2ntpair2freq, dumpster)
