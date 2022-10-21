# Utilities for dealing with GFF3 files.

from strainflye import config
from strainflye.errors import ParameterError


def validate_basic(
    feature,
    contig,
    contig_length,
    seen_feature_ids,
    fancylog,
    zero_indexed_range=True,
):
    """Performs basic validation on a feature in a GFF3 file.

    Parameters
    ----------
    feature: skbio.metadata.Interval
        Feature on which we are performing validation.

    contig: str
        Name of the contig on which this feature is located.

    contig_length: int
        Length of "contig".

    seen_feature_ids: set
        Set of seen feature IDs on this contig. We'll check this to make sure
        that the ID of "feature" has not already been seen in this set; also,
        we'll add the ID of "feature" to this set.

    fancylog: function
        Logging function. As of writing, we just use this to warn about
        ambiguous 1-position features.

    zero_indexed_range: bool
        If True, return a zero-indexed position range for this feature; if
        False, return a one-indexed position range for this feature.

    Returns
    -------
    (fid, feature_range): (str, range)

        fid: ID of this feature.

        feature_range: Inclusive range of all positions this feature spans.

    Raises
    ------
    ParameterError
        If various things seem wrong with this feature.
    """
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

    if not zero_indexed_range:
        feature_range = range(fs + 1, fe + 1)

    return fid, feature_range


def check_cds_attrs(cds_id, contig, cds_length, strand):
    """Sanity-checks attributes of a coding sequence.

    Parameters
    ----------
    cds_id: str
        ID of this coding sequence. This isn't checked or anything, but it is
        used in the resulting error message, if applicable.

    contig: str
        Contig on which this CDS is located. As with cds_id, this is only used
        in the error message, if applicable.

    cds_length: int
        Length of the CDS. We check that this is divisible by 3.

    strand: str
        Strand of the CDS. We check that this is one of {"+", "-"}.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If the strand or length checks fail.
    """
    if strand != "+" and strand != "-":
        raise ParameterError(
            f'Feature {cds_id} on contig {contig} has a strand of "{strand}". '
            'We require here that all CDS features have a strand of "+" or '
            '"-".'
        )

    if cds_length % 3 != 0:
        raise ParameterError(
            f"Feature {cds_id} on contig {contig} has a length of "
            f"{cds_length:,} bp. We require here that all CDS features have "
            "lengths divisible by 3."
        )


def validate_if_cds(feature, contig, verboselog):
    """Checks if a feature is a CDS, and if so performs extra validation.

    "CDS" stands for coding sequence, I think. Or at least it would be pretty
    awkward at this point if it didn't.

    Parameters
    ----------
    feature: skbio.metadata.Interval
        Feature on which we are performing validation.

    contig: str
        Name of the contig on which this feature is located.

    verboselog: function
        Logging function. As of writing, we just use this to let the user know
        that we're ignoring non-CDS features.

    Returns
    -------
    (is_cds, strand): (bool, str)

        is_cds: True if this feature is a CDS, False otherwise.

        strand: If is_cds is True, then this will be one of {"+", "-"}. If
                is_cds is False, then this will be None.

    Raises
    ------
    ParameterError
        If various things seem wrong with this feature.
    """
    fid = feature.metadata["ID"]
    # Ignore features that don't aren't labelled as "CDS"s
    # As far as I can tell, this is the only type (or at least the main
    # one) that is used to identify coding sequences in GFF3 files.
    # "Gene" is a possible type, for example, but it seems mostly to be
    # used as a higher-level type for top-level genes containing
    # multiple CDSs -- for reference, Prodigal's output seems to
    # consist solely of CDSs.
    if feature.metadata["type"].upper() not in config.CDS_TYPES_LIST:
        verboselog(
            (
                f"Feature {fid} on contig {contig} has a type that is not in "
                f"{config.CDS_TYPES}; ignoring this feature."
            ),
            prefix="",
        )
        return False, None

    # The phase can be used to indicate an offset -- i.e. phase = 0
    # indicates that (for + strand genes) the first position is the
    # first position of the first codon, or (for - strand genes) the
    # last position is the last position of the first codon (after
    # complementing). phase = 1 indicates that you should ignore this
    # first position, and that the first codon starts one to the right
    # (+ genes) or one to the left (- genes). Same with phase = 2.
    #
    # At least, this is my interpretation from going through the GFF3
    # docs -- e.g.
    # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    #
    # However: as far as I can tell, Prodigal always outputs "phase 0"
    # CDSs, and I'm not sure how exactly to interpret CDSs with other
    # phases (do these correspond to incomplete genes? why are we
    # including these extra nucleotides?). They raise some questions
    # about if, for example, the length % 3 check makes sense. So, for
    # the sake of simplicity, I'm just gonna mandate "phase 0" CDSs for
    # now.
    if "phase" in feature.metadata:
        phase = feature.metadata["phase"]
        if phase != 0:
            # As of writing, scikit-bio's GFF3 parser actually allows any int
            # (I think) to be the phase for a feature, rather than restricting
            # this to one of {0, 1, 2}. If we ever choose to support phases of
            # 1 or 2, we'd need to add an extra check that would throw an error
            # if the phase isn't one of these three.
            raise ParameterError(
                f"Feature {fid} on contig {contig} has a phase of {phase}. "
                "This command does not support features with non-zero phases, "
                "at least for now."
            )
    else:
        # As of writing, if you set the phase to "." in a GFF3 file (how you'd
        # indicate "no phase" per the GFF3 spec), scikit-bio's GFF3 parser
        # will not assign a phase attribute to this feature. (... In fairness,
        # GFF3's spec mandates that all CDS features have a phase, so this
        # shouldn't happen much.)
        raise ParameterError(
            f"Feature {fid} on contig {contig} does not have a phase "
            "attribute; this is required (more specifically, required to be "
            "0) here for all CDS features."
        )

    strand = feature.metadata["strand"]
    cds_len = len(range(*feature.bounds[0]))
    check_cds_attrs(fid, contig, cds_len, strand)

    return True, strand
