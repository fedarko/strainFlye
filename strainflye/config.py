from strainflye.cli_utils import list2str

# Max depth at any position in the alignment: needed for pysam, see the
# pysamstats docs -- https://github.com/alimanfoo/pysamstats
MAX_DEPTH_PYSAM = 100000000

# Prefix of diversity index columns in diversity index TSV files
DI_PREF = "DivIdx"

# Default parameters for minimap2 for "strainFlye align"
DEFAULT_MM2_PARAMS = "-ax asm20 --secondary=no -I 8g --MD"

# Default parameters for LJA for "strainFlye smooth assemble"
DEFAULT_LJA_PARAMS = "--simpleec --Cov-threshold 10"

# Rationale: this essentially the power set of (CP2, Tv, Nonsyn, Nonsense).
# We ignore:
#
# - combinations where Nonsyn and Nonsense are together (because all nonsense
#   mutations are by definition nonsynonymous).
#
# - CP2 + Tv + Nonsyn (because there is only one synonymous mutation in CP2 in
#   the standard genetic code, TAA <--> TGA, and A <--> G is a transition; so
#   CP2Tv and CP2TvNonsyn would be identical).
DECOY_CONTEXTS = [
    "Full",
    "CP2",
    "Tv",
    "Nonsyn",
    "Nonsense",
    "CP2Tv",
    "CP2Nonsyn",
    "CP2Nonsense",
    "TvNonsyn",
    "TvNonsense",
    "CP2TvNonsense",
]

# Optional value for the --decoy-contexts parameter for FDR estimation:
# if this is given, we'll just use all decoy contexts (that way the user
# doesn't have to write each option out).
DCTX_EVERYTHING = "Everything"

DCTX_ALL = DECOY_CONTEXTS + [DCTX_EVERYTHING]

# Used in the "strainFlye link" module. Probably not an important optimization
# but humor me here, ok?
N2I = {"A": 0, "C": 1, "G": 2, "T": 3}
I2N = "ACGT"

# Used in "strainFlye link", both when naming info files and when later
# searching for these files in a directory (so, please don't change this unless
# you really need to, because that'd make old "strainFlye link nt" outputs
# incompatible with future versions of "strainFlye link graph"...)
POS_FILE_LBL = "pos2nt2ct"
POSPAIR_FILE_LBL = "pospair2ntpair2ct"

# When going through contigs' gene predictions to compute mutation matrices,
# we'll focus on predicted genes that have this as their type in the GFF3 file:
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#
# This may be too conservative, but we'll try to be loud about ignoring these
# features so it should be ok.
CDS_TYPES_LIST = ["CDS", "SO:0000316"]
# Fancy version we can use for logging messages
CDS_TYPES = list2str(CDS_TYPES_LIST)

# Pre-computed information about codons and reverse complements. I figure it
# makes more sense to store this here than it does to re-compute this every
# time the user needs to run "matrix count" or whatever lol
CODONS = [
    "AAA",
    "AAC",
    "AAG",
    "AAT",
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "AGA",
    "AGC",
    "AGG",
    "AGT",
    "ATA",
    "ATC",
    "ATG",
    "ATT",
    "CAA",
    "CAC",
    "CAG",
    "CAT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CGA",
    "CGC",
    "CGG",
    "CGT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GAA",
    "GAC",
    "GAG",
    "GAT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GGA",
    "GGC",
    "GGG",
    "GGT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TAA",
    "TAC",
    "TAG",
    "TAT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "TGA",
    "TGC",
    "TGG",
    "TGT",
    "TTA",
    "TTC",
    "TTG",
    "TTT",
]
CODON2RC = {
    "AAA": "TTT",
    "AAC": "GTT",
    "AAG": "CTT",
    "AAT": "ATT",
    "ACA": "TGT",
    "ACC": "GGT",
    "ACG": "CGT",
    "ACT": "AGT",
    "AGA": "TCT",
    "AGC": "GCT",
    "AGG": "CCT",
    "AGT": "ACT",
    "ATA": "TAT",
    "ATC": "GAT",
    "ATG": "CAT",
    "ATT": "AAT",
    "CAA": "TTG",
    "CAC": "GTG",
    "CAG": "CTG",
    "CAT": "ATG",
    "CCA": "TGG",
    "CCC": "GGG",
    "CCG": "CGG",
    "CCT": "AGG",
    "CGA": "TCG",
    "CGC": "GCG",
    "CGG": "CCG",
    "CGT": "ACG",
    "CTA": "TAG",
    "CTC": "GAG",
    "CTG": "CAG",
    "CTT": "AAG",
    "GAA": "TTC",
    "GAC": "GTC",
    "GAG": "CTC",
    "GAT": "ATC",
    "GCA": "TGC",
    "GCC": "GGC",
    "GCG": "CGC",
    "GCT": "AGC",
    "GGA": "TCC",
    "GGC": "GCC",
    "GGG": "CCC",
    "GGT": "ACC",
    "GTA": "TAC",
    "GTC": "GAC",
    "GTG": "CAC",
    "GTT": "AAC",
    "TAA": "TTA",
    "TAC": "GTA",
    "TAG": "CTA",
    "TAT": "ATA",
    "TCA": "TGA",
    "TCC": "GGA",
    "TCG": "CGA",
    "TCT": "AGA",
    "TGA": "TCA",
    "TGC": "GCA",
    "TGG": "CCA",
    "TGT": "ACA",
    "TTA": "TAA",
    "TTC": "GAA",
    "TTG": "CAA",
    "TTT": "AAA",
}
