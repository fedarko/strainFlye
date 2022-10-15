# Utilities for strainFlye dynam.

from strainflye import misc_utils


def run_covskew(contigs, bam, binlen, norm_cov_epsilon, output_dir, fancylog):
    """Computes coverage and skew information for contigs.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a (sorted and indexed) BAM file mapping reads to contigs.

    binlen: int
        Bin length.

    norm_cov_epsilon: float
        Used to determine the clamp "height" for normalized coverage.

    output_dir: str
        Directory to which we'll write out TSV files for each contig.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    contig_name2len, bf, num_fasta_contigs = misc_utils.load_fasta_and_bam(
        contigs, bam, fancylog
    )
