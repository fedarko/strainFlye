# Max depth at any position in the alignment: needed for pysam, see the
# pysamstats docs -- https://github.com/alimanfoo/pysamstats
MAX_DEPTH_PYSAM = 100000000

# Prefix of diversity index columns in diversity index TSV files
DI_PREF = "DivIdx"

# Default parameters for LJA for "strainFlye smooth assemble"
DEFAULT_LJA_PARAMS = "--simpleec --Cov-threshold 10"

# Rationale: this essentially the power set of (CP2, Tv, Nonsyn, Nonsense),
# ignoring combinations where Nonsyn and Nonsense are together (because all
# nonsense mutations are by definition nonsynonymous).
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
    "CP2TvNonsyn",
    "CP2TvNonsense",
]
