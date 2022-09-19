# Utilities for strainFlye spot.


import skbio
from math import floor
from scipy.special import comb
from strainflye import bcf_utils
from strainflye.errors import ParameterError, WeirdError


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
        - If various things are invalid in the BCF file

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

    bcf_obj, bcf_contigs = bcf_utils.loudly_parse_arbitrary_bcf_and_contigs(
        bcf, fancylog
    )

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
    at_least_one_feature_on_input_contigs_seen = False
    for contig, im in contig_and_im_tuples:

        at_least_one_feature_seen = True

        # Skip features that are not located on contigs described in the BCF
        # file. Not raising an error here is nice, in the case where our BCF
        # file describes contigs that are a subset of those in the GFF3 file --
        # for example, we have gene predictions for *all* contigs in a dataset,
        # but we only bothered to do mutation calling on some of these contigs.
        if contig not in bcf_contigs:
            continue
        at_least_one_feature_on_input_contigs_seen = True

        # Using the BCF file, record all mutated positions in this contig
        mutated_positions = bcf_utils.get_mutated_positions_in_contig(
            bcf_obj, contig
        )

        # Also, get its length -- this'll be useful for checking that the
        # feature coordinates are sane later.
        # (bcf_utils.parse_arbitrary_bcf() should have already raised an error
        # if length was unavailable in the BCF header for any contigs, so this
        # is safe.)
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

    if not at_least_one_feature_on_input_contigs_seen:
        raise ParameterError(
            "None of the feature(s) described in the GFF3 file are located "
            "on contigs that are described in the BCF file."
        )

    fancylog(
        (
            f"Identified {len(hotspots):,} hotspot feature(s) across all "
            f"{len(bcf_contigs):,} contigs in the BCF file."
        ),
        prefix="",
    )
    fancylog("Writing out this information to a TSV file...")

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
    fancylog("Done.", prefix="")


def get_coldspot_gaps_in_contig(muts, contig_length, min_length, circular):
    """Identifies "coldspot gaps" in a certain contig.

    This is a utility function, to make testing coldspot stuff easier.

    Parameters
    ----------
    muts: list
        Sorted list of (1-indexed) mutated positions in a contig. We assume
        that all of these are in the inclusive range [1, contig_length].

    contig_length: int
        Length of the contig.

    min_length: int
        Minimum length of a gap to be considered a coldspot. We assume this is
        greater than 0.

    circular: bool
        If True, assume that each contig is circular, and examine the presence
        of a potential gap between the rightmost mutation and the leftmost
        mutation in the contig. If False, don't do this.

    Returns
    -------
    coldspots: list of (int, int, int)
        Each entry in this list is a 3-tuple that describes a coldspot gap
        identified in this contig. The four entries of each tuple are:

          1. Start position of the gap (1-indexed, inclusive)
          2. End position of the gap (1-indexed, inclusive)
          3. Length of the gap (aka total number of positions in the gap)
    """
    coldspots = []

    # There is almost certainly a big brain way to handle this elegantly
    # but my brain is smooth and my coffee mug is empty
    # These functions take care of modulo arithmetic for us when circular
    # is True
    def r1c(pos):
        n = pos + 1
        if n > contig_length:
            return 1
        else:
            return n

    def l1c(pos):
        p = pos - 1
        if p < 1:
            return contig_length
        else:
            return p

    # Silly corner cases
    if len(muts) == 0:
        if contig_length >= min_length:
            coldspots.append((1, contig_length, contig_length))

    elif len(muts) == 1:
        # We can think of the contig as the following diagram:
        #
        #    L       R
        # -------M-------
        #
        # Where M indicates the position of the one mutated position.
        m = muts[0]
        if circular:
            # The only gap is from M + 1 to M - 1 -- spanning all R and L
            # positions.
            if (contig_length - 1) >= min_length:
                start = r1c(m)
                end = l1c(m)
                coldspots.append((start, end, contig_length - 1))

        else:
            # So, now we gotta treat L and R as two separate gaps.

            # If m == 1, then m - 1 == 0, and we know that min_length > 0.
            # So we'll automatically ignore a possible gap in this case.
            llen = m - 1
            if llen >= min_length:
                coldspots.append((1, llen, llen))

            # Similar idea: if m == contig_length, then (contig_length - m)
            # == 0, and no possible gap is considered.
            rlen = contig_length - m
            if rlen >= min_length:
                coldspots.append((m + 1, contig_length, rlen))

    else:
        # Look for easy, "internal" gaps between mutations
        for mpi in range(1, len(muts)):
            prev_mut = muts[mpi - 1]
            curr_mut = muts[mpi]

            if prev_mut != curr_mut - 1:
                # This is a gap of length (curr_mut - prev_mut - 1), since
                # we exclude both the start and end point of the "gap."
                # For example, if prev_mut = 4 and curr_mut = 6, then we
                # have a gap of length 1 (including just position 5).
                gap_len = curr_mut - prev_mut - 1

                if gap_len >= min_length:
                    gap_start = prev_mut + 1
                    gap_end = curr_mut - 1
                    coldspots.append((gap_start, gap_end, gap_len))

        if circular:
            # Test the loop-around gap from the rightmost mutation to the
            # leftmost # mutation
            loop_gap_len = (contig_length + muts[0]) - muts[-1] - 1
            if loop_gap_len >= min_length:
                coldspots.append((r1c(muts[-1]), l1c(muts[0]), loop_gap_len))
        else:
            # Test the gaps from 1 to muts[0] and
            # from muts[-1] to contig_length. This is analogous to how we
            # handle the 1-mutation case when circular is False above.
            # (Maybe I should reuse that code here... whatevs.)

            # Test the left gap
            llen = muts[0] - 1
            if llen >= min_length:
                coldspots.append((1, llen, llen))

            # Test the right gap
            rlen = contig_length - muts[-1]
            if rlen >= min_length:
                coldspots.append((muts[-1] + 1, contig_length, rlen))
    return coldspots


def longest_success_run_pvalue(m, n, p):
    """Computes Prob(X >= m), where X is the longest run of "successes."

    "Success" is used here in terms of a Bernoulli trial. I'm just gonna
    straight up copy from (Naus 1982):

    "Given a sequence of n Bernoulli trials with probability of success
    on a given trial p, and q = 1 - p, let X denote the length of the
    longest success run."

    We can think of Prob(X >= m) as a p-value: it's the probability that the
    longest run of successes is at least m, under the null hypothesis that
    successes and failures (in our sequence of n Bernoulli trials) are in
    "random order." (This is detailed more in (Bateman 1948).)

    Parameters
    ----------
    m: int
        The longest number of "successes" seen in a row in our sequence of
        n Bernoulli trials.

    n: int
        The number of Bernoulli trials in our sequence.

    p: float
        The probability of success for each Bernoulli trial. Should be in the
        range [0, 1].

    Returns
    -------
    pval: float
        The probability of the longest run of successes in our n Bernoulli
        trials having a length of at least m, under the aforementioned null
        hypothesis.

    Raises
    ------
    WeirdError
        If various things about the input numbers look wrong:
        - m >= n
        - m < 1
        - p < 0 or p > 1

    References
    ----------
    - Bateman, G. (1948). On the Power Function of the Longest Run as a Test
      for Randomness in a Sequence of Alternatives. Biometrika, 35(1/2),
      97-112.

    - Glaz, J., Naus, J. I., Wallenstein, S., Wallenstein, S., & Naus, J. I.
      (2001). Scan Statistics (pp. 243-259). New York: Springer.

    - Naus, J. I. (1982). Approximations for Distributions of Scan Statistics.
      Journal of the American Statistical Association, 77(377), 177-183.

    We use equation (3.1), given in (Naus 1982). This particular equation
    originates from (Bateman 1948). (And I found out about (Naus 1982) from the
    book "Scan Statistics.")
    """
    if m >= n:
        raise WeirdError("n must be greater than m.")

    if m < 1:
        raise WeirdError("m must be at least 1.")

    if p < 0 or p > 1:
        raise WeirdError("p must be in the range [0, 1].")

    q = 1 - p
    # So, we know that floor(n / m) must be at least 1. Since the endpoint
    # of summations are inclusive (just, like, in general --
    # http://www.columbia.edu/itc/sipa/math/summation.html -- I feel like
    # an idiot looking this up but at least I'm a [probably] correct idiot)
    # we know that we will update sum_term at least once, since
    # list(range(1, 2)) == [1].
    sign = 1
    pval = 0
    for j in range(1, floor(n / m) + 1):
        jm = j * m
        pval += (
            sign
            * (p + (((n - jm + 1) * q) / j))
            * comb(n - jm, j - 1)
            * (p**jm)
            * (q ** (j - 1))
        )
        # Corresponds to (-1) ** (j + 1) in the equation.
        sign = -sign

    # Deal with potential precision silliness by clamping to [0, 1]
    # I've seen that sum_term can exceed 1 in some cases (e.g. if
    # num_muts == 5, contig_length == 100, and the max gap length == 2);
    # not sure if it getting less than 0 is possible, but might as well
    # prevent it.
    #
    # (This is probably a harbinger that it would be better to use e.g. log
    # probabilities or something to try to limit precision problems. But
    # these p-values are already a very basic, trivial implementation, and
    # there are like 100 other things on my plate for just this project
    # ._.)
    if pval > 1:
        pval = 1
    elif pval < 0:
        pval = 0

    return pval


def get_coldspot_gap_pvalues(num_muts, contig_length, coldspot_lengths):
    """Computes a p-value for the longest coldspot gap in a contig.

    Parameters
    ----------
    num_muts: int
        Number of mutated positions in a contig. Must be <= contig_length.

    contig_length: int
        Length of the contig.

    coldspot_lengths: list of int
        List of coldspot lengths.

    Returns
    -------
    pvals: list
        Has the same number of entries as coldspot_lengths. The longest
        coldspot in this list (breaking ties arbitrarily) will be assigned
        a p-value, if it is possible to compute this probability for this
        coldspot's length; all other coldspots won't have p-values given.

        What does this p-value represent? Define the "null hypothesis" as the
        case where mutated positions are randomly placed on the contig. The
        p-value reported for a given gap, then, is the probability (given this
        null hypothesis) that the longest gap we would see in a contig is at
        least as long as this gap's length.

    Raises
    ------
    WeirdError
        If various things seem inconsistent with the input parameters.

    References
    ----------
    - Geller, R., Domingo-Calap, P., Cuevas, J. M., Rossolillo, P., Negroni,
      M., & Sanjuán, R. (2015). The external domains of the HIV-1 envelope are
      a mutational cold spot. Nature Communications, 6(1), 1-9.

    The way that we estimate the probability of a given position being
    mutated (under the null hypothesis of randomly distributed mutations) is
    analogous to (Geller, Domingo-Calap, Cuevas et al., 2015) -- see
    https://www.nature.com/articles/ncomms9571#Sec13 (they use a different
    method of defining coldspots with a sliding window approach, but they set
    the probability of mutation in the same way [at least, as far as I can
    tell]).

    Also, see the references for longest_success_run_pvalue() for more details.

    Notes
    -----
    - It's probably possible to give p-values for the second-largest gap size,
      third-largest, etc. But I don't know how to do that and also I don't feel
      like doing it.

    - So, we will report at most one p-value. If the longest gap length is
      shared by multiple gaps, then which of these gets assigned a p-value is
      arbitrary.

    - We implicitly make the assumption that, if the longest gap corresponds to
      the "loop-around" gap obtained by setting circular = True when
      identifying gaps in get_coldspot_gaps_in_contig(), that the contig just
      loops around. I guess this is analogous to saying "let's just shift the
      'endpoints' of the circle so that we can treat the section of it
      containing the longest gap like a line."

    - So you know that null hypothesis we defined above, where mutations are
      placed randomly on a contig? That gets broken basically immediately.
      Like, the codon position figure in our paper shows that there is bias in
      where mutations occur. (So I would not recommend reading too much into
      these p-values...)

    - We don't perform much validation on the number of coldspot gaps -- half
      to make testing easier, and half because this should never get messed up
      in practice. (Also, the number of gaps isn't straightforward to compute
      from the number of mutations and contig length, due to things like
      circular maybe being True.)
    """
    pvals = []
    if num_muts == 0:
        if len(coldspot_lengths) != 1 or coldspot_lengths[0] != contig_length:
            raise WeirdError(
                "A contig with 0 mutations must have exactly 1 coldspot "
                "covering the entire contig."
            )
        # We have no way to estimate the mutation rate of this contig.
        pvals.append("NA")

    elif num_muts == contig_length:
        if len(coldspot_lengths) != 0:
            raise WeirdError(
                "There can't be any gaps if every position is mutated."
            )
    else:
        # Get the index of the largest gap (breaking ties arbitrarily)
        # this is probably a bit inefficient (could bundle the computation of
        # index and len into one iteration through coldspot_lengths) but
        # whatevs
        num_gaps = len(coldspot_lengths)
        max_gap_idx = max(
            range(0, num_gaps), key=lambda i: coldspot_lengths[i]
        )

        # The equation we will use is for the longest run of "successes" in
        # a sequence of n Bernoulli trials. We thus define a "success" as a
        # position *not* being mutated (we could also define this as a failure,
        # it doesn't really matter).
        #
        # If the null hypothesis (mutations occur randomly on the sequence) is
        # true, then all positions have the same probability of being mutated.
        # We can "estimate" the probability of a position being mutated as
        # the number of mutations in the contig divided by the total number of
        # positions in the contig -- this is analogous to how this probability
        # is set in Geller, Domingo-Calap, Cuevas et al., 2015 (see refs
        # above).
        p = 1 - (num_muts / contig_length)

        # Say "NA" for all gaps but the longest one
        # (Since we break ties arbitrarily, this ignores the fact that there
        # may be multiple gaps tied for the "longest.")
        pvals = ["NA"] * num_gaps
        pvals[max_gap_idx] = longest_success_run_pvalue(
            coldspot_lengths[max_gap_idx], contig_length, p
        )
    return pvals


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
    search to find something that matches the analysis we do in the paper.
    Something something reinventing the wheel.

    Parameters
    ----------
    bcf: str
        Filepath to a BCF file describing single-nucleotide mutations.

    min_length: int
        Minimum length of a gap to be considered a coldspot. Should be > 0.

    circular: bool
        If True, assume that each contig is circular, and examine the presence
        of a potential gap between the rightmost mutation and the leftmost
        mutation in the contig. If False, don't do this.
        (This is a lazy way to handle this, since not all contigs will be
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

    Raises
    ------
    ParameterError
        - If min_length < 1
        - If various things are invalid in the BCF file
    """
    # Should never happen due to Click enforcing this, but some of our logic
    # here breaks if min_length == 0 so we manually enforce this to be paranoid
    if min_length < 1:
        raise ParameterError("Minimum coldspot gap length must be at least 1.")

    bcf_obj, bcf_contigs = bcf_utils.loudly_parse_arbitrary_bcf_and_contigs(
        bcf, fancylog
    )

    # Maps contig IDs to a list of 3-tuples of (start, end, length)
    # ... where start/end are 1-indexed and inclusive
    # (i.e. the start should be the first position after a mutation, and the
    # end should be the last position before a mutation, assuming that this gap
    # falls in the interior of a contig and we don't have to worry about
    # looping around stuff)
    #
    # Note for future devs: the reason we store length for these
    # coldspots is that reverse-engineering the length of a coldspot is
    # actually kind of annoying if it loops around the contig (b/c then we
    # gotta know the contig length, and it's a whole ordeal). Easier to just
    # be a little bit memory inefficient and store lengths from the get-go.
    contig2coldspots = {}

    # Maps contig IDs to a list of p-values for each coldspot. Each list should
    # have the same length as the corresponding list in contig2coldspots.
    contig2coldspot_pvals = {}

    total_num_coldspots = 0

    fancylog("Going through contigs and identifying coldspot gaps...")
    for contig in bcf_contigs:

        contig_length = bcf_obj.header.contigs[contig].length

        # (These are 1-indexed positions)
        muts = sorted(
            bcf_utils.get_mutated_positions_in_contig(
                bcf_obj, contig, zero_indexed=False
            )
        )

        # We defer the actual process of finding gaps to a utility function to
        # make testing easier (also to make this code look pretty)
        contig_coldspots = get_coldspot_gaps_in_contig(
            muts, contig_length, min_length, circular
        )
        contig2coldspots[contig] = contig_coldspots

        # The p-value computation is a separate step, just to avoid making me
        # update all of the other tests (it might be slightly faster to bundle
        # this into one loop, but this would take more work)
        coldspot_pvals = get_coldspot_gap_pvalues(
            len(muts), contig_length, [cs[2] for cs in contig_coldspots]
        )
        contig2coldspot_pvals[contig] = coldspot_pvals

        total_num_coldspots += len(contig_coldspots)

    fancylog(
        (
            f"Identified {total_num_coldspots:,} coldspot gap(s) across all "
            f"{len(bcf_contigs):,} contigs in the BCF file."
        ),
        prefix="",
    )
    fancylog("Writing out this information to a TSV file...")

    with open(output_coldspot_gaps, "w") as of:
        of.write(
            "Contig\tStart_1IndexedInclusive\tEnd_1IndexedInclusive\tLength\t"
            "P_Value\n"
        )
        for contig in contig2coldspots.keys():
            both_lists = zip(
                contig2coldspots[contig], contig2coldspot_pvals[contig]
            )
            for (start, end, length), pval in both_lists:
                of.write(f"{contig}\t{start}\t{end}\t{length}\t{pval}\n")
    fancylog("Done.", prefix="")
