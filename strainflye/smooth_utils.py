# Utilities for strainFlye smooth.

from strainflye import phasing_utils, cli_utils


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
    """
    contig_name2len, bam_obj, bcf_obj = phasing_utils.load_triplet(
        contigs, bam, bcf, fancylog
    )
    num_fasta_contigs = len(contig_name2len)
    if virtual_reads:
        rt = "smoothed and virtual reads"
    else:
        rt = "smoothed reads"

    fancylog(f"Going through contigs and constructing {rt}...")
    for si, seq in enumerate(contig_name2len, 1):
        contig_len = contig_name2len[seq]
        if verbose:
            cli_utils.log_prog(
                seq, contig_len, si, num_fasta_contigs, fancylog
            )
        # Task 1: iterate through all alignments to this contig and create
        # smoothed reads. Also compute average coverages.
        # (Unlike the Phasing-LJA.ipynb notebook, we don't have average
        # coverages up front; I guess we could have users pass in the diversity
        # index TSV file which includes these, but that's gross and
        # precludes arbitrary inputs.)

        # Task 2: iterate through runs of low-coverage positions and create
        # virtual reads.

    # misc_utils.make_output_dir(output_dir)


def run_assemble():
    pass
