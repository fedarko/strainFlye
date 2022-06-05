# This isn't a comprehensive list of parameter descriptions, at least not yet.
# The main goal is abstracting repeated stuff (e.g. between p-mutation and
# r-mutation calling).

INPUT_CONTIGS_NAIVE_CALL = (
    "FASTA file of contigs in which to na\u00efvely call mutations. "
    "All contigs in this FASTA file should also be contained in the BAM file; "
    "however, if some contigs in the BAM file are not included in this FASTA "
    "file, we won't perform any calling on these absent contigs."
)

INPUT_BAM = "BAM file representing an alignment of reads to contigs."

OUTPUT_VCF_NAIVE_CALL = (
    "Filepath to which an output VCF file (describing the called "
    "mutations) will be written."
)

VERBOSE_CALL = "Display extra details for each contig while running."