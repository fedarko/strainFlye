# This isn't a comprehensive list of parameter descriptions, at least not yet.
# The main goal is abstracting repeated stuff (e.g. between p-mutation and
# r-mutation calling).

INPUT_CONTIGS_NAIVE_CALL = (
    "FASTA file of contigs in which to na\u00efvely call mutations. "
    "All contigs in this FASTA file should also be contained in the BAM file; "
    "however, if some contigs in the BAM file are not included in this FASTA "
    "file, we won't perform any calling on these absent contigs."
)

INPUT_BAM = (
    "Sorted and indexed BAM file representing an alignment of reads to "
    "contigs."
)

OUTPUT_DIR_NAIVE_CALL = (
    "Directory to which an output BCF file (describing the called "
    "mutations), BCF index file, and diversity index TSV file will be "
    "written. Some temporary files will also be written to this directory."
)

VERBOSE_CALL = "Display extra details for each contig while running."

INPUT_CONTIGS = "FASTA file of contigs."

INPUT_BCF_DOWNSTREAM = (
    "Indexed BCF file describing single-nucleotide mutations in a set of "
    "contigs."
)
