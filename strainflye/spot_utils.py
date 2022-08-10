# Utilities for strainFlye spot.


import skbio
from strainflye import bcf_utils
from strainflye.errors import ParameterError


def get_mutated_positions_in_contig(bcf_obj, contig):
    mutated_positions = set()
    for mut in bcf_obj.fetch(contig):
        # "positions" in pysam are 1-based; to be compatible with scikit-bio's
        # 0-based positions (from parsing GFF files), subtract 1
        mutated_positions.add(mut.pos - 1)
    return mutated_positions


def run_hotspot_detection(
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
        a hotspot. Should be greater than 0.

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
        "Going through features in the GFF3 file and identifying hotspots..."
    )
    # Load features. skbio.io.read() with GFF3 files defaults to a generator of
    # tuples, where the first entry is the sequence ID (e.g. "edge_1") and the
    # second entry is a skbio.metadata.IntervalMetadata object describing all
    # features within this sequence.
    contig_and_im_tuples = skbio.io.read(features, format="gff3")
    for contig, im in contig_and_im_tuples:
        if contig not in bcf_contigs:
            raise ParameterError(
                "The GFF3 file describes feature(s) located on contig "
                f"{contig}, but this contig is not described in the BCF file."
            )

        # Using the BCF file, record all mutated positions in this contig
        mutated_positions = get_mutated_positions_in_contig(bcf_obj, contig)

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
            # ranges, etc. based on these coordinates.
            #
            # (Just for clarity: we access .bounds[0] because there should
            # only be one [start, end] interval per feature -- otherwise our
            # check above on len(feature.bounds) would have thrown an error.)
            feature_range = range(*feature.bounds[0])

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

    fancylog(
        f"Identified {len(hotspots):,} hotspots across all contigs.", prefix=""
    )
    fancylog("Writing out this information to a file...")

    with open(output_hotspot_features, "w") as of:
        of.write(
            "Contig\tFeatureID\tFeatureStart_1IndexedInclusive\t"
            "FeatureEnd_1IndexedInclusive\tNumberMutatedPositions\t"
            "PercentMutatedPositions\n"
        )
        for contig, feature, num_mp, perc_mp in hotspots:
            start = feature.bounds[0][0] + 1
            end = feature.bounds[0][1]
            of.write(
                f"{contig}\t{feature.metadata['ID']}\t{start}\t{end}\t"
                f"{num_mp}\t{perc_mp:.2f}%\n"
            )
