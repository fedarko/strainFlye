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
    # second entry is a skbio.metadata.IntervalMetadata object describing this
    # feature.
    contig_and_im_tuples = skbio.io.read(features, format="gff3")
    for contig, im in contig_and_im_tuples:
        if contig not in bcf_contigs:
            raise ParameterError(
                "The features file describes feature(s) located on contig "
                f"{contig}, but this contig is not described in the BCF file."
            )

        # TODO: record all mutated positions in this contig

        # HACK: iterate through all features that belong to this contig.
        # See https://github.com/biocore/scikit-bio/issues/1817 regarding why
        # we can't just say im.query(), at least at the moment.
        for feature in im.query(metadata={}):
            if len(feature.bounds) < 1:
                raise ParameterError(
                    "A feature exists that doesn't have at least one set of "
                    f"bounds: {feature}"
                )

            # TODO: assert that bounds don't overlap with each other, right?
            # or just compute a set of positions correspodning to all positions
            # in this feature, then "&" that with the set of mutated positions.
            # could be sped up but probs good enough for now
