# This isn't a comprehensive list of parameter descriptions, at least not yet.
# The main goal is abstracting repeated stuff (e.g. between p-mutation and
# r-mutation calling).

# these two aren't stand-alone descriptions, but are instead eesigned to be
# used in the descriptions of various downstream tasks
FASTA_SUBSET_BAM = (
    "All contigs in this FASTA file should also be contained in the BAM file; "
    "it's ok if the BAM file contains contigs not in this FASTA file (we'll "
    "ignore them)."
)

FASTA_SUBSET_BAM_BCF = (
    "All contigs in this FASTA file should also be contained in the BAM and "
    "BCF files; it's ok if the BAM or BCF files contain contigs not in this "
    "FASTA file (we'll ignore them)."
)

INPUT_CONTIGS_NAIVE_CALL = (
    "FASTA file of contigs in which to na\u00efvely call mutations. "
    f"{FASTA_SUBSET_BAM}"
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

VERBOSE_BASIC = "Display extra details while running."

INPUT_BCF_DOWNSTREAM = (
    "Indexed BCF file describing single-nucleotide mutations in a set of "
    "contigs."
)
