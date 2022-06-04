# Max depth at any position in the alignment: needed for pysam, see the
# pysamstats docs -- https://github.com/alimanfoo/pysamstats
MAX_DEPTH_PYSAM = 10000000

# We consider a p-mutation with freq(pos) >= this percentage to be
# "high-frequency", aka indisputable. Depending on what we are using mutation
# calling for, we may want to not call high-frequency mutated positions (for
# example, this will throw off the computation of decoy genomes in the context
# of FDR estimation). However, most of the time, we probably want to include
# these as real mutations, since... they are real mutations! (Probably.)
#
# This should be a number in the range (0, 50], analogous to the "p" value that
# is a parameter for p-mutation calling. (I guess you could set it to exactly
# zero, but that would cause literally every position to be a mutation, so...
# look, don't do that, okay?)
HIGH_FREQUENCY_MIN_PCT = 5

# Prefix of diversity index columns in diversity index TSV files
DI_PREF = "DivIdx"
