# Shared utilities between strainFlye smooth and strainFlye link-graph.


import pysam
from strainflye import fasta_utils, bcf_utils, misc_utils


def load_triplet(contigs, bam, bcf, fancylog):
    """Loads and checks three files: FASTA, BAM, and BCF.

    Mainly, this ensures that the contigs in the FASTA file are all present in
    the BAM and BCF files.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a (sorted and indexed) BAM file mapping reads to contigs.

    bcf: str
        Filepath to an (indexed) BCF file describing single-nucleotide
        mutations in contigs.

    fancylog: function
        Logging function.

    Returns
    -------
    (contig_name2len, bam_obj, bcf_obj): (dict, pysam.AlignmentFile,
                                          pysam.VariantFile)
        contig_name2len: dict mapping contig name to length.
        bam_obj: Object describing the BAM file.
        bcf_obj: Object describing the BCF file.

    Raises
    ------
    This function doesn't raise any errors itself, but it calls various
    functions which can raise errors if the input files are invalid in certain
    ways. See fasta_utils.get_name2len(), misc_utils.verify_contigs_subset(),
    and bcf_utils.parse_bcf() for more details on the sorts of errors that can
    get raised here.
    """
    fancylog("Loading and checking FASTA, BAM, and BCF files...")

    contig_name2len = fasta_utils.get_name2len(contigs)
    fasta_contigs = set(contig_name2len)
    fancylog(
        f"The FASTA file describes {len(fasta_contigs):,} contigs.", prefix=""
    )

    bam_obj = pysam.AlignmentFile(bam, "rb")
    misc_utils.verify_contigs_subset(
        fasta_contigs,
        set(bam_obj.references),
        "the FASTA file",
        "the BAM file",
    )
    fancylog(
        (
            "All FASTA contigs are included in "
            f"the BAM file (this BAM file has {bam_obj.nreferences:,} "
            "references)."
        ),
        prefix="",
    )

    bcf_obj, thresh_type, thresh_min = bcf_utils.parse_bcf(bcf)
    bcf_contigs = set(bcf_obj.header.contigs)
    misc_utils.verify_contigs_subset(
        fasta_contigs,
        bcf_contigs,
        "the FASTA file",
        "the BCF file",
    )
    fancylog(
        (
            "All FASTA contigs are included in "
            "the BCF file (the header of this BCF file describes "
            f"{len(bcf_contigs):,} contigs)."
        ),
        prefix="",
    )
    fancylog("So far, these files seem good.", prefix="")

    return contig_name2len, bam_obj, bcf_obj
