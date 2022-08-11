# Utilities for strainFlye spot.


import skbio
from strainflye import bcf_utils
from strainflye.errors import ParameterError


def run_hotspot_feature_detection(
    bcf,
    features,
    min_num_mutations,
    min_perc_mutations,
    output_hotspot_features,
    fancylog,
):
    """Identifies "hotspots," based on user-defined threshold(s).

    Writes out information about these hotspots to a TSV file.

    Parameters
    ----------
    bcf: str
        Filepath to a BCF file describing single-nucleotide mutations.

    features: str
        Filepath to a GFF3 file describing "features" in the dataset's contigs.
        Not all contigs in the BCF file need to be represented in the feature
        file, but all contigs in the feature file must be represented in the
        BCF file.

    min_num_mutations: int or None
        Minimum number of mutations that a feature must have to be considered
        a hotspot. Should be > 0.

    min_perc_mutations: float or None
        Minimum percentage of mutations that a feature must have to be
        considered a hotspot. If this is not None, this should be in the range
        (0, 100].

    output_hotspot_features: str
        Filepath to which we'll write a TSV file describing identified hotspot
        features.

    fancylog: function
        Logging function.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        - If min_num_mutations and min_perc_mutations are both None
        - If various things are invalid in the GFF3 file

    References
    ----------
    For details on GFF3 files (used for the features parameter), see
    http://gmod.org/wiki/GFF3.

    GFF3's start and end coordinates are 1-based and inclusive on both the
    start and end. So, for example, a feature in a GFF3 file with
    start = 5 and end = 19 spans exactly 15 characters (in interval notation,
    this is [5, 19]).

    Note that the docs linked above (at least as of writing) are ambiguous
    about whether the *end* coordinates of features are inclusive or exclusive.
    Some of the difficulty I encountered in interpreting this is due to this
    sentence in the GFF3 docs:

      "For zero-length features, such as insertion sites, start equals end and
      the implied site is to the right of the indicated base in the direction
      of the landmark."

    To me, this implied that "one-length" features would thus have end = start
    + 1 -- meaning that, in general, the end coordinate of each feature in a
    GFF3 file should be exclusive. But this is not the case: as is confirmed
    in this excellent article by Daniel Standage,
    https://standage.github.io/on-genomic-interval-notation.html, **there is no
    (easy) distinction in GFF3 files between features of lengths zero and
    one.** We thus make the assumption that all features where start = end are
    of length one, and warn the user using fancylog() if any such features are
    present in the input file. (But, like, come on.)

    ... SO ANYWAY hopefully these paragraphs save you the ~hour I spent trying
    to figure this out lmao. For a moment there I was like "ok but what if i
    just drop out of grad school"
    """
    if min_num_mutations is None and min_perc_mutations is None:
        raise ParameterError(
            "At least one of (--min-num-mutations, --min-perc-mutations) must "
            "be specified."
        )

    fancylog("Loading and checking the BCF file...")
    # Load BCF
    bcf_obj, thresh_type, thresh_min = bcf_utils.parse_bcf(bcf)
    fancylog("Looks good so far.", prefix="")
    bcf_contigs = set(bcf_obj.header.contigs)

    # List of tuples of
    # (contig ID, feature, # mutations in feature, % mutations in feature)
    # ...where "feature" is a skbio.Interval object.
    hotspots = []

    fancylog(
        "Going through features in the GFF3 file and identifying "
        "hotspot features..."
    )
    # Load features. skbio.io.read() with GFF3 files defaults to a generator of
    # tuples, where the first entry is the sequence ID (e.g. "edge_1") and the
    # second entry is a skbio.metadata.IntervalMetadata object describing all
    # features within this sequence.
    contig_and_im_tuples = skbio.io.read(features, format="gff3")
    at_least_one_feature_seen = False
    for contig, im in contig_and_im_tuples:
        at_least_one_feature_seen = True
        if contig not in bcf_contigs:
            raise ParameterError(
                "The GFF3 file describes feature(s) located on contig "
                f"{contig}, but this contig is not described in the BCF file."
            )

        # Using the BCF file, record all mutated positions in this contig
        mutated_positions = bcf_utils.get_mutated_positions_in_contig(
            bcf_obj, contig
        )

        # Also, get its length -- this'll be useful for checking that the
        # feature coordinates are sane later.
        # (NOTE: We implicitly make the assumption here that the input BCF file
        # has lengths given in its header contig tags, even though the VCF 4.3
        # specification -- as of writing -- only says that these tags
        # "typically" include this information. So, if this starts leading to
        # AttributeErrors, we can handle this another way -- likely just by
        # having the users pass contigs to this command, so we can extract
        # lengths from there. But eesh.)
        contig_length = bcf_obj.header.contigs[contig].length

        seen_feature_ids = set()

        # Next, iterate through all features that belong to this contig.
        # Our use of im.query(metadata={}) to do this is kind of a hack: see
        # https://github.com/biocore/scikit-bio/issues/1817 regarding why
        # we can't just say im.query(), at least at the moment.
        #
        # NOTE: alternatively, we could figure out the first and last mutated
        # position in this contig and use that as the "bounds" for im.query()?
        # That should work, since we should be able to rely on the fact that
        # features containing zero mutated positions would never be called as
        # hotspots. But... let's not get too crazy optimizing this for now.
        for feature in im.query(metadata={}):
            # We could definitely add support for multi-boundary (e.g.
            # discontiguous) features if desired, but for the sake of
            # simplicity we don't right now.
            #
            # GFF3 can encode hierarchies of features that happen to be
            # discontiguous -- e.g. a parent gene that contains multiple exon
            # features -- but we would still treat each of those exons as
            # independent features. GFF3 can also have discontiguous features
            # if multiple rows share the same ID. Currently, we will just treat
            # features with different IDs as independent features entirely; and
            # if multiple rows share the same ID, we'll raise an error.
            # See "Nesting Features" and "Discontinuous Features",
            # respectively, at http://gmod.org/wiki/GFF3.
            #
            # (NOTE: as of writing, skbio's GFF3 parser doesn't seem to output
            # any features [Intervals] with multiple bounds at once, so this
            # branch is technically untested. Maybe we could move this to
            # another function that takes as input an Interval object [which
            # would simplify unit-testing], but ... that's probs not necessary
            # right now.)
            if len(feature.bounds) != 1:
                raise ParameterError(
                    f"A feature in the GFF3 file on contig {contig} exists "
                    f"without exactly one set of bounds: {feature}"
                )

            # This should rarely happen. scikit-bio's GFF3 parser makes the
            # implicit assumption that at least one key=value pair is defined
            # in the final "attributes" column in each row -- so, if this
            # column is a "." for any feature, then scikit-bio will fail at
            # this line:
            # https://github.com/biocore/scikit-bio/blob/541498807b67554353fc8aeb65bb66c28966a1f6/skbio/io/format/gff3.py#L446
            # ... However, if a row has attributes defined, but if none of
            # these defined attributes is "ID", then yeah -- we'll run into
            # this problem. And *this* case is something we test against.
            if "ID" not in feature.metadata:
                raise ParameterError(
                    f"A feature in the GFF3 file on contig {contig} exists "
                    f"without a defined ID: {feature}"
                )

            # Enforce that feature IDs are unique with respect to their contig
            # (It's ok if features in different contigs happen to have
            # different IDs, I guess, although this shouldn't happen with
            # Prodigal output). Nonunique feature IDs are *probably* an
            # indication that the GFF3 file describes discontinuous features,
            # which we explicitly do not support yet. So raise an error.
            fid = feature.metadata["ID"]
            if fid in seen_feature_ids:
                raise ParameterError(
                    f"The feature ID {fid} is used in multiple GFF3 rows for "
                    f"contig {contig}. Features of a contig must have unique "
                    "IDs; this command does not support "
                    '"discontinuous features" at the moment.'
                )
            else:
                seen_feature_ids.add(fid)

            # scikit-bio's GFF3 parser ensures that the start coordinate of a
            # feature must be <= the end coordinate. So we can safely create
            # ranges, etc. based on these coordinates. (Also, these bounds are
            # zero-indexed and half-open, like Python intervals, so they can be
            # directly be used as the inputs to range().)
            #
            # (Just for clarity: we access .bounds[0] because there should
            # only be one [start, end] interval per feature -- otherwise our
            # check above on len(feature.bounds) would have thrown an error.)
            fs, fe = feature.bounds[0]
            feature_range = range(fs, fe)

            # See https://standage.github.io/on-genomic-interval-notation.html,
            # and the references above. This is what we in the business refer
            # to as a certified Bioinformatics Moment (TM) (no one actually
            # says this, also hello i'm surprised someone is reading this code)
            if len(feature_range) == 1:
                fancylog(
                    (
                        f"Warning: feature {fid} on contig {contig} has equal "
                        "start and end coordinates. We assume this refers to "
                        "a feature of length 1 spanning this single position, "
                        "rather than a feature of length 0."
                    ),
                    prefix="",
                )

            # This indicates that the GFF3 file is malformed, or maybe designed
            # for a different contig with this same name.
            # Recall that fs and fe are zero-indexed and half-open -- since
            # it's ok for a feature to start at the last position in a contig
            # (if, for example, it has a length of one), the first "invalid"
            # start position is at contig_length
            if fs >= contig_length:
                # Although we use 0-indexing internally, the user's GFF3 file
                # uses 1-indexing. So let's make this easy for them to
                # understand.
                raise ParameterError(
                    f"Feature {fid} on contig {contig} has a (1-indexed) "
                    f"start coordinate of {fs + 1:,}, which is greater than "
                    f"the contig's length of {contig_length:,}."
                )

            # Unlike the above case (feature start past the contig length),
            # having the feature end past the contig length is actually
            # possible in valid GFF3 files. This indicates that this
            # feature is circular, and "loops around" the contig.
            #
            # (fe is zero-indexed and "half-open", so its first "past the
            # contig length" position is at contig_length + 1. If, on the other
            # hand, fe == contig_length, this just indicates that the feature
            # ends exactly at the rightmost position in the contig.)
            #
            # TODO: add support for this eventually? Since this is explicitly
            # allowed in the GFF spec, and could conceivably happen with
            # prokaryotic gene predictions or whatevs.
            if fe > contig_length:
                raise ParameterError(
                    f"Feature {fid} on contig {contig} has a (1-indexed) end "
                    f"coordinate of {fe:,}, which is greater than the "
                    f"contig's length of {contig_length:,}. We do not support "
                    "'circular' features yet."
                )

            # This does the main work -- figure out what positions within the
            # "range" given by this feature are mutated.
            num_mp = len(set(feature_range) & mutated_positions)
            perc_mp = 100 * (num_mp / len(feature_range))

            # We know how many mutations this feature contains, as well as the
            # length of this feature. Determine whether or not this feature is
            # a "hotspot," at least based on how the user has defined the term.

            # Check 1: a feature is a hotspot if it has at least N mutations
            if min_num_mutations is not None:
                if num_mp < min_num_mutations:
                    continue

            # Check 2: a feature is a hotspot if its percentage of mutations is
            # at least P
            if min_perc_mutations is not None:
                if perc_mp < min_perc_mutations:
                    continue

            # If we've made it here, this feature is a hotspot. Yay!
            hotspots.append((contig, feature, num_mp, perc_mp))

    if not at_least_one_feature_seen:
        raise ParameterError(
            "The GFF3 file doesn't seem to describe any features."
        )

    fancylog(
        f"Identified {len(hotspots):,} hotspot feature(s) across all contigs.",
        prefix="",
    )
    fancylog("Writing out this information to a file...")

    with open(output_hotspot_features, "w") as of:
        of.write(
            "Contig\tFeatureID\tFeatureStart_1IndexedInclusive\t"
            "FeatureEnd_1IndexedInclusive\tNumberMutatedPositions\t"
            "PercentMutatedPositions\n"
        )
        for contig, feature, num_mp, perc_mp in hotspots:
            # Note that we output feature boundaries in 1-indexed coordinates,
            # to be consistent with the input GFF3 start and end
            start = feature.bounds[0][0] + 1
            end = feature.bounds[0][1]
            of.write(
                f"{contig}\t{feature.metadata['ID']}\t{start}\t{end}\t"
                f"{num_mp}\t{perc_mp:.2f}%\n"
            )


def run_coldspot_gap_detection(
    bcf,
    min_length,
    circular,
    output_coldspot_gaps,
    fancylog,
):
    """Identifies "coldspot gaps" based on a user-defined threshold.

    Writes out information about these coldspots to a TSV file.

    This is exactly what it says on the tin. There are definitely more
    sophisticated ways to do this, but it's faster for me to just briefly write
    up and test this command than it is to perform a thorough literature
    search. Something something reinventing the wheel.

    Parameters
    ----------
    bcf: str
        Filepath to a BCF file describing single-nucleotide mutations.

    min_length: int
        Minimum length of a gap to be considered a coldspot. Should be > 0.

    circular: bool
        If True, assume that each contig is circular, and examine the presence
        of a potential gap between the rightmost mutation and the leftmost
        mutation in the contig. (If the contig has length N, then we'd treat
        the first position of the contig as N + 1, etc.) If False, don't do
        this. (This is a lazy way to handle this, since not all contigs will be
        circular -- ideally, we'd handle this on a contig-by-contig basis,
        maybe by consulting the assembly graph to see which contigs actually
        are circular.)

    output_coldspot_gaps: str
        Filepath to which we'll write a TSV file describing identified
        coldspots.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    fancylog("Marcus needs to stop being lazy and implement this.")
