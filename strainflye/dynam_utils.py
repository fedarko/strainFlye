# Utilities for strainFlye dynam.


import os
import pysamstats
from math import floor
from statistics import median
from strainflye import misc_utils, fasta_utils, cli_utils, config


def skew(dna):
    """Returns a single numeric value: (G - C) / (G + C) for a DNA sequence.

    Parameters
    ----------
    dna: str
        DNA sequence. This should *not* be a skbio.DNA object, since we are
        going to check individual characters versus "C" or "G".

    Returns
    -------
    skew: float
        (G - C) / (G + C) for dna, where G = number of Gs and C = number of Cs.

        If G + C == 0, this value is undefined -- for convenience's sake, we
        return 0 in this case.

    References
    ----------
    Computed as discussed in Grigoriev 1998, "Analyzing genomes with cumulative
    skew diagrams."
    """
    g_ct = 0
    c_ct = 0
    for nt in dna:
        if nt == "C":
            c_ct += 1
        if nt == "G":
            g_ct += 1

    # Done out of convenience -- I'm not sure what the best practice is here
    if g_ct + c_ct == 0:
        return 0

    return (g_ct - c_ct) / (g_ct + c_ct)


def update_cumulative_binned_skews(cumulative_binned_skews, bin_skew):
    """Updates a list of cumulative binned skews.

    Basically: if this list is empty, we just add bin_skew to it. If this list
    isn't empty, then we compute the *cumulative* skew by adding the rightmost
    entry in it and bin_skew; then we add that to the right of the list.
    """
    if len(cumulative_binned_skews) == 0:
        cumulative_binned_skews.append(bin_skew)
    else:
        cumulative_binned_skews.append(cumulative_binned_skews[-1] + bin_skew)


def contig_covskew(contig, contigs, bam_obj, bin_len, nclb, ncub):
    """Computes normalized binned coverages, binned skews, and bin centers.

    Parameters
    ----------
    contig: str
        Name of a contig for which we'll compute this information.

    contigs: str
        Filepath to a FASTA file describing contigs. We assume that "contig" is
        described in this file.

    bam_obj: pysam.AlignmentFile
        Describes an alignment of reads to contigs. We assume that "contig" is
        one of the references in this alignment.

    bin_len: int
        Bin length.

    nclb: float
        Lower bound to which we'll clamp normalized coverages.

    ncub: float
        Upper bound to which we'll clamp normalized coverages.

    Returns
    -------
    (nb_coverages, cb_skews, left_positions, center_positions): (list, list,
                                                                 list, list)

        Define "b" as the number of bins for this contig.

        nb_coverages: Contains b entries. The bi-th entry describes the
                      normalized coverage for the bi-th bin.

        cb_skews: Contains b entries. The bi-th entry describes the cumulative
                  GC skew for the bi-th bin.

        left_positions: Contains b entries. The bi-th entry describes the
                        leftmost position in this bin, using 1-indexing.
                        (I don't think we need to include the rightmost
                        position, also, because that should be easy to figure
                        out.)

        center_positions: Contains b entries. The bi-th entry describes the
                          center of this bin, computed as the average of the
                          leftmost and rightmost positions included in this
                          bin (using 1-indexing). Useful for plotting.

    Note: How do we compute binned coverages and normalize them?
    ------------------------------------------------------------

    1. Compute the median coverage in each bin.

    2. Compute the median of these medians, agg_total_cov (we call it "M" in
       the paper).

    3. Divide each bin's median coverage by agg_total_cov. This gives us
       "normalized" coverages.

    4. Clamp normalized coverages to the range [nclb, ncub].
    """
    covs = []
    for rec in pysamstats.stat_variation(
        bam_obj,
        chrom=contig,
        fafile=contigs,
        pad=True,
        max_depth=config.MAX_DEPTH_PYSAM,
    ):
        covs.append(rec["A"] + rec["C"] + rec["G"] + rec["T"])

    # NOTE: Using get_single_seq() like this, in the middle of a loop over all
    # contigs, is inefficient. If this becomes a bottleneck, can do indexing or
    # something. See notes on this in smooth_utils.write_smoothed_reads().
    contig_seq = str(fasta_utils.get_single_seq(contigs, contig))

    contig_len = len(contig_seq)

    binned_coverages = []
    cumulative_binned_skews = []
    left_positions = []
    center_positions = []
    left_pos = 1
    while left_pos + bin_len - 1 <= contig_len:
        left_positions.append(left_pos)

        # The -1 is needed to fit things in properly.
        # For example, say our sequence is ABCDEFGHIJKLMNOP (start: 1, end: 16)
        # Bins w/ length 3:                1234567890123456
        #                                  ---===---===---=
        # The [left pos, right pos] intervals (inclusive on both ends) are:
        # [1, 3], [4, 6], [7, 9], [10, 12], [13, 15], [16]
        right_pos = left_pos + bin_len - 1
        center_positions.append((left_pos + right_pos) / 2)

        # We need to convert left_pos and right_pos from 1- to 0-indexing in
        # order to use them as indices in covs. (There's no "- 1" done on
        # right_pos because Python's half-open indexing means that we are
        # already excluding right_pos.)
        bin_covs = covs[left_pos - 1 : right_pos]

        binned_coverages.append(median(bin_covs))

        # Compute binned skew while we're computing binned coverage. In the
        # original analysis notebook these were two separate functions, but
        # it makes sense to do them all at once.
        bin_skew = skew(contig_seq[left_pos - 1 : right_pos])
        update_cumulative_binned_skews(cumulative_binned_skews, bin_skew)

        left_pos = right_pos + 1

    # Unless contig_len was evenly divisible by bin_len, there will be some
    # extra positions not in any bins yet (at the right end of the sequence).
    # Create a new bin to hold these.
    if left_pos <= contig_len:
        left_positions.append(left_pos)
        positions_in_bin = range(left_pos, contig_len + 1)
        center_positions.append(
            (positions_in_bin[0] + positions_in_bin[-1]) / 2
        )

        # compute this special tiny bin's coverage
        bin_covs = covs[left_pos - 1 :]
        binned_coverages.append(median(bin_covs))

        # ... and skew
        bin_skew = skew(contig_seq[left_pos - 1 :])
        update_cumulative_binned_skews(cumulative_binned_skews, bin_skew)

    # We've got binned coverages -- do normalization now.
    agg_total_cov = median(binned_coverages)

    if agg_total_cov == 0:
        # If the median of medians is zero, then we can't do normalization.
        # Set all normalized coverages to "NA" and leave it at that -- we can
        # still give the user skew information, though. (This might actually
        # happen in practice.)
        norm_binned_coverages = ["NA"] * len(binned_coverages)
    else:
        # Thankfully, most cases should ... not have the above case happen.
        # We can actually divide by agg_total_cov.
        norm_binned_coverages = []

        for bc in binned_coverages:
            norm_cov = bc / agg_total_cov

            # Clamp norm_cov to [lower bound, upper bound] if needed
            if norm_cov < nclb:
                norm_cov = nclb
            elif norm_cov > ncub:
                norm_cov = ncub

            norm_binned_coverages.append(norm_cov)

    return (
        norm_binned_coverages,
        cumulative_binned_skews,
        left_positions,
        center_positions,
    )


def log_bin_ct_info(contig, contig_len, bin_len, verboselog):
    """Logs information about the number of bins in a contig."""

    if contig_len >= bin_len:
        # This is the number of "full" bins of exactly bin_len
        num_bins = floor(contig_len / bin_len)
        bin_noun = "bins" if num_bins > 1 else "bin"
        bin_desc = f"{num_bins:,} {bin_noun} of length {bin_len:,} bp"

        # There will probably be extra positions on the right (e.g. if
        # contig_len = 23 and bin_len = 10, then we'll have two "full" bins and
        # three positions remaining outside of these bins). So, let's
        # create an extra bin containing these extra positions.
        rbl = contig_len % bin_len
        if contig_len % bin_len != 0:
            num_bins += 1
            bin_desc += f" and 1 smaller bin of length {rbl:,} bp"
    else:
        # Make a bin of length contig_len including all positions in the contig
        num_bins = 1
        bin_desc = f"1 smaller bin of length {contig_len:,} bp"

    verboselog(f"Creating {bin_desc} for contig {contig}...", prefix="")


def run_covskew(
    contigs, bam, bin_len, norm_cov_epsilon, output_dir, verbose, fancylog
):
    """Computes coverage and skew information for contigs.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a (sorted and indexed) BAM file mapping reads to contigs.

    bin_len: int
        Bin length.

    norm_cov_epsilon: float
        Used to determine the clamp range for normalized coverages: we'll set
        it as [1 - norm_cov_epsilon, 1 + norm_cov_epsilon].

    output_dir: str
        Directory to which we'll write out TSV files for each contig.

    verbose: bool
        Log extra info.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    verboselog = cli_utils.get_verboselog(fancylog, verbose)
    contig_name2len, bam_obj, num_contigs = misc_utils.load_fasta_and_bam(
        contigs, bam, fancylog
    )
    ncl = 1 - norm_cov_epsilon
    ncu = 1 + norm_cov_epsilon
    fancylog(
        "Going through contigs and computing coverage/skew information..."
    )
    misc_utils.make_output_dir(output_dir)
    for ci, contig in enumerate(contig_name2len, 1):

        clen = contig_name2len[contig]
        cli_utils.proglog(contig, ci, num_contigs, verboselog, contig_len=clen)

        if verbose:
            log_bin_ct_info(contig, clen, bin_len, verboselog)

        (
            nb_coverages,
            b_skews,
            left_positions,
            center_positions,
        ) = contig_covskew(contig, contigs, bam_obj, bin_len, ncl, ncu)
        with open(os.path.join(output_dir, f"{contig}_covskew.tsv"), "w") as f:
            f.write(
                "LeftPos_1IndexedInclusive\t"
                "CenterPos\t"
                "NormalizedCoverage\t"
                "CumulativeSkew\n"
            )
            all_lists = zip(
                nb_coverages, b_skews, left_positions, center_positions
            )
            for (cov, skew, left, center) in all_lists:
                f.write(f"{left}\t{center}\t{cov}\t{skew}\n")

    fancylog("Done.", prefix="")
