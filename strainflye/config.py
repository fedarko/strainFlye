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

MIN_PENWIDTH = 0.01
MAX_PENWIDTH = 5
