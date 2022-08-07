# Utilities for strainFlye spot.


import skbio
from strainflye import bcf_utils
from strainflye.errors import ParameterError


def find_hotspot_features(
    bcf,
    features,
    min_num_mutations,
    min_perc_mutations,
    output_hotspot_features,
):
    if min_num_mutations is None and min_perc_mutations is None:
        raise ParameterError(
            "At least one of (--min-num-mutations, --min-perc-mutations) must "
            "be specified."
        )

    # Load BCF
    bcf_obj, thresh_type, thresh_min = bcf_utils.parse_bcf(bcf)
    bcf_contigs = set(bcf_obj.header.contigs)

    # Load features. skbio.io.read() with GFF3 files defaults to a generator of
    # tuples, where the first entry is the sequence ID (e.g. "edge_1") and the
    # second entry is a skbio.metadata.IntervalMetadata object describing all
    # features within this sequence.
    contig_and_im_tuples = skbio.io.read(features, format="gff3")
    for contig, im in contig_and_im_tuples:
        if contig not in bcf_contigs:
            raise ParameterError(
                "The features file describes feature(s) located on contig "
                f"{contig}, but this contig is not described in the BCF file."
            )

        # TODO: Using the BCF file, record all mutated positions in this contig

        # Next, iterate through all features that belong to this contig.
        # Our use of im.query(metadata={}) to do this is kind of a hack: see
        # https://github.com/biocore/scikit-bio/issues/1817 regarding why
        # we can't just say im.query(), at least at the moment.
        for feature in im.query(metadata={}):
            if len(feature.bounds) < 1:
                raise ParameterError(
                    "A feature exists that doesn't have at least one set of "
                    f"bounds: {feature}"
                )

            # TODO: The minimal required functionality here is just figuring
            # out how many of the mutated positions in this contig lie within
            # this feature's bounds, and then using the min_num / min_perc
            # parameters to classify this feature as a "hotspot" or not.
            #
            # Ideally, we would make this really sophisticated (so that we
            # could do things like use the "Parent" GFF3 attribute to link
            # together exon features from the same gene, or something) ... but
            # I don't think there is demand for that here yet, so for now we
            # will just treat each "feature" in the GFF3 file as separate from
            # all other features.

    # TODO: write out classified hotspots
