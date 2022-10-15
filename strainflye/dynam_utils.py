# Utilities for strainFlye dynam.

import pysamstats
from math import floor
from statistics import median
from strainflye import misc_utils, cli_utils, config


def compute_binned_coverages(contig, contigs, bam_obj, bin_len, verboselog):
    """Computes binned coverages (and the centers of these bins) for a contig.

    Parameters
    ----------
    contig: str
        Name of a contig.

    contigs: str
        Filepath to a FASTA file describing contigs. (pysamstats needs this.)

    bam_obj: pysam.AlignmentFile
        Describes an alignment of reads to contigs. We assume that "contig" is
        one of the references in this alignment.

    verboselog: function
        Logging function.

    Returns
    -------
    (binned_coverages, center_positions): (list, list)

        binned_coverages: Contains b entries, where b is the number of bins for
                          this contig. The bi-th entry describes the median
                          coverage of all positions in the bi-th bin.

        center_positions: Contains b entries, where b is the number of bins for
                          this contig. The bi-th entry describes the center of
                          this bin, computed as the average of the leftmost and
                          rightmost positions included in this bin. Useful for
                          plotting.
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

    # we could just pass contig_name2len here, but this is easier (and we've
    # already guaranteed thanks to misc_utils that this matches what
    # contig_name2len says)
    contig_len = bam_obj.get_reference_length(contig)

    binned_coverages = []
    center_positions = []
    left_pos = 1
    while left_pos + bin_len - 1 <= contig_len:

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
        assert len(bin_covs) == bin_len

        binned_coverages.append(median(bin_covs))
        left_pos = right_pos + 1

    # Unless contig_len was evenly divisible by bin_len, there will be some
    # extra positions not in any bins yet (at the right end of the sequence).
    # Create a new bin to hold these.
    if left_pos <= contig_len:
        positions_in_bin = range(left_pos, contig_len + 1)
        center_positions.append(
            (positions_in_bin[0] + positions_in_bin[-1]) / 2
        )
        bin_covs = covs[left_pos - 1 :]
        binned_coverages.append(median(bin_covs))

    return binned_coverages, center_positions


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
        Used to determine the clamp "height" for normalized coverage.

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
    fancylog(
        "Going through contigs and computing coverage/skew information..."
    )
    for ci, contig in enumerate(contig_name2len, 1):

        clen = contig_name2len[contig]
        cli_utils.proglog(contig, ci, num_contigs, verboselog, contig_len=clen)

        if clen >= bin_len:
            # This is the number of "full" bins of exactly bin_len
            num_bins = floor(clen / bin_len)
            bin_noun = "bins" if num_bins > 1 else "bin"
            bin_desc = f"{num_bins:,} {bin_noun} of length {bin_len:,} bp"

            # There will probably be extra positions on the right (e.g. if
            # clen = 23 and bin_len = 10, then we'll have two "full" bins and
            # three positions remaining outside of these bins). So, let's
            # create an extra bin containing these extra positions.
            rbl = clen % bin_len
            if clen % bin_len != 0:
                num_bins += 1
                bin_desc += f" and 1 smaller bin of length {rbl:,} bp"
        else:
            # Make a bin of length clen containing all positions in the contig
            num_bins = 1
            bin_desc = f"1 smaller bin of length {clen:,} bp"

        verboselog(f"Creating {bin_desc} for contig {contig}...", prefix="")
        binned_coverages, center_positions = compute_binned_coverages(
            contig, contigs, bam_obj, bin_len, verboselog
        )
        verboselog(binned_coverages, prefix="")
        verboselog(center_positions, prefix="")
