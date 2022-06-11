# Utilities for strainFlye fdr.


import re
import pysam
import pandas as pd
from math import floor
from collections import defaultdict
from strainflye import fasta_utils, call_utils
from .errors import ParameterError, SequencingDataError
from .config import DI_PREF


def parse_bcf(bcf):
    """Opens a BCF file, and does sanity checking and p vs. r sniffing on it.

    Thankfully, pysam has the ability to read BCF files, so the main thing
    we do here is checking that the meta-information of the BCF file seems
    seems kosher (i.e. was produced by strainFlye).

    Parameters
    ----------
    bcf: str
        Path to a BCF file produced by one of "strainFlye call"'s subcommands.

    Returns
    -------
    (f, thresh_type, thresh_min): (pysam.VariantFile, str, int)
        f: Object describing the input BCF file.
        thresh_type: either "p" or "r", depending on what type of mutation
                     calling was done to produce this BCF file.
        thresh_min: the minimum value of p or r used in mutation calling to
                    produce this BCF file. Formatted the same as used in the
                    file (so values of p will still be scaled up by 100).

    Raises
    ------
    FileNotFoundError
        If bcf doesn't exist (raised by pysam).

    ValueError
        If bcf doesn't look like a VCF/BCF file (raised by pysam).
        (Note that we could *technically* accept gzipped and indexed VCF files,
        I guess, but I don't want to officially add support for that because
        that sounds like a lot of testing.)

    ParameterError
        If bcf does not have exactly one "strainFlye threshold filter header,"
        which is a term I made up just now. Basically, we rely on there
        existing a single line in the BCF file's header that goes like

            ##FILTER=<ID=strainflye_minT_MMMM,...>

        ... where T corresponds to the type of mutation calling done (p or r)
        and MMMM (variable number of digits) indicates the minimum value of p
        or r used in mutation calling.

        If there isn't exactly one of these lines, then we will be very
        confused! Hence why we raise an error.
    """
    # this will fail with a FileNotFoundError if "bcf" doesn't point to an
    # existing file (although we shouldn't need to worry about that much b/c
    # click should've already checked that this file exists); and it'll fail
    # with a ValueError if it points to a file but this file doesn't look like
    # a VCF/BCF file
    f = pysam.VariantFile(bcf)

    # Now that this at least seems like a VCF/BCF file, make sure it's from
    # strainFlye, and figure out whether it's from p- or r-mutation calling...
    thresh_type = None
    thresh_min = None

    # pysam seems to add an extra filter labelled PASS to the parsed BCF file,
    # for some reason. So let's consider all filters that the file has -- might
    # as well, because it's useful to detect the weird case where there are > 1
    # strainFlye threshold headers (we raise an error about this below)
    for filter_name in f.header.filters:
        filter_match = re.match(r"^strainflye_min([pr])_(\d+)$", filter_name)
        if filter_match is not None:
            # We should only see a filter with this type of ID once. If we see
            # it multiple times, something has gone very wrong.
            if thresh_type is not None or thresh_min is not None:
                raise ParameterError(
                    f"BCF file {bcf} has multiple strainFlye threshold filter "
                    "headers."
                )
            thresh_type = filter_match.group(1)
            thresh_min = int(filter_match.group(2))

    # If we never updated these variables, we never saw a strainFlye filter
    # header -- probably this BCF isn't from strainFlye.
    if thresh_type is None or thresh_min is None:
        raise ParameterError(
            f"BCF file {bcf} doesn't seem to be from strainFlye: no threshold "
            "filter headers."
        )

    # Let's be extra paranoid and verify that this BCF has MDP (coverage based
    # on (mis)matches) and AAD (alternate nucleotide coverage) fields. (It
    # should, because we know at this point that strainFlye generated it, but
    # you never know...)
    info_ids = f.header.info.keys()
    if "MDP" not in info_ids or "AAD" not in info_ids:
        raise ParameterError(
            f"BCF file {bcf} needs to have MDP and AAD info fields."
        )

    return f, thresh_type, thresh_min


def check_decoy_selection(diversity_indices, decoy_contig):
    """Checks that only one of (diversity index file, decoy contig) is given.

    Parameters
    ----------
    diversity_indices: str or None
        If a str, this should be a filepath to a TSV file representing
        diversity index info.

    decoy_contig: str or None
        If a str, this should be the name of a contig described in the
        BCF file.

    Returns
    -------
    selection_type: str
        Will be "DI" if only diversity_indices is not None, and will be "DC" if
        only decoy_contig is not None.

    Raises
    ------
    ParameterError
        - If diversity_indices and decoy_contig are both not None
        - If diversity_indices and decoy_contig are both None
    """
    di = diversity_indices is not None
    dc = decoy_contig is not None

    if di:
        if dc:
            raise ParameterError(
                "Both the diversity indices file and a decoy contig are "
                "specified. These options are mutually exclusive."
            )
        else:
            return "DI"
    else:
        if dc:
            return "DC"
        else:
            raise ParameterError(
                "Either the diversity indices file or a decoy contig must be "
                "specified."
            )


def normalize_series(in_series):
    """Converts a series to values in the range [0, 1].

    Parameters
    ----------
    in_series: pd.Series
        We assume that this does not contain any nan values.

    Returns
    -------
    None or pd.Series
        If the minimum and maximum of in_series are identical, this will just
        return None. (In this case, we can't scale values, because the
        denominator we use when converting a value (max - min) is zero.)

        If the minimum and maximum are not identical (which should usually
        be the case with diversity indices, hopefully...) then this will return
        a pd.Series with the same index as in_series, but with each entry
        scaled to within the range [0, 1] (such that the min value in in_series
        is set to 0, the max is set to 1, and everything else is in between).
    """
    # Small TODO: in theory, it'd be faster to combine the
    # computation of min and max into a single pass over the values
    # (see e.g. https://stackoverflow.com/q/12200580) but this
    # probably won't be a performance bottleneck so I'm not gonna
    # bother for now
    min_val = min(in_series)
    max_val = max(in_series)
    if min_val == max_val:
        return None
    else:
        # Use pandas' vectorization to apply linear interpolation
        # across all diversity indices in this Series
        return (in_series - min_val) / (max_val - min_val)


def autoselect_decoy(diversity_indices, min_len, min_avg_cov, fancylog):
    """Attempts to select a good decoy contig based on diversity index data.

    There are lots of ways to implement this, so here we just stick with
    something simple that combines the diversity index information from
    multiple thresholds:

    1. Filter to all contigs whose lengths and average coverages meet the
       specified thresholds. The number of "passing" contigs is C. If C = 1,
       select this contig; if C = 0, raise an error.

    2. For each of the D diversity index columns provided in the file (where at
       least two contigs have defined diversity indices), compute the minimum
       and maximum diversity index in this column. Assign each contig a score
       for this column in [0, 1] using linear interpolation: the contig with
       the lowest diversity index gets a score of 0, the contig with the
       highest gets a score of 1, and everything else is scaled in between.
       If a contig has an undefined diversity index in such a column, set its
       score for this column to 1.

    3. For each of the C contigs, sum scores across all of the D diversity
       index columns. Select the contig with the lowest score sum. Break ties
       arbitrarily.

    Parameters
    ----------
    diversity_indices: str
        Filepath to a TSV file containing diversity index info, generated by
        one of strainFlye call's subcommands. In addition to diversity index
        values, this also includes length and average coverage information
        for each contig in the file (this will help a lot, since computing
        these values if we don't already have them is time-consuming or at
        the very least annoying).

    min_len: int
        In order for a contig to be selected as the decoy, its length must be
        at least this.

    min_avg_cov: float
        In order for a contig to be selected as the decoy, its average coverage
        must be at least this.

    fancylog: function
        Logging function.

    Returns
    -------
    decoy_contig: str
        The name of the decoy contig we find.

    Raises
    ------
    ParameterError
        If the diversity index file:
        - Describes < 2 contigs (should have already been caught during align,
          but you never know)
        - Doesn't have Length or AverageCoverage columns

    SequencingDataError
        - If none of the contigs in the diversity index file pass the length
          and average coverage thresholds.
        - If none of the diversity index columns has at least two "passing"
          contigs with defined and distinct diversity indices in this column.
    """
    di = pd.read_csv(diversity_indices, sep="\t", index_col=0)

    # We raise an error here that covers the just-one-contig case, too. Because
    # although we could select that contig as a decoy, we wouldn't have any
    # target contigs left!
    if len(di.index) < 2:
        raise ParameterError("Diversity indices file describes < 2 contigs.")

    if "Length" not in di.columns or "AverageCoverage" not in di.columns:
        raise ParameterError(
            "Diversity indices file must include the Length and "
            "AverageCoverage columns."
        )

    # Filter to contigs that pass both the length and coverage thresholds.
    # https://stackoverflow.com/a/13616382
    passing_di = di[
        (di["Length"] >= min_len) & (di["AverageCoverage"] >= min_avg_cov)
    ]
    passing_contigs = passing_di.index
    num_passing_contigs = len(passing_contigs)
    # it isn't clear how much precision min_avg_cov has, so we don't impose any
    # limit on how many digits it goes out to when printing it out. let's let
    # python handle this one
    check_str = (
        f"the min length \u2265 {min_len:,} and min average cov \u2265 "
        f"{min_avg_cov:,}x checks"
    )
    if num_passing_contigs == 0:
        raise SequencingDataError(f"No contigs pass {check_str}.")
    elif num_passing_contigs == 1:
        # Arguably, we could raise an error here -- but we know at this point
        # that there are >= 2 contigs in the file (and this is just the only
        # one that passes the length and coverage thresholds). So we may as
        # well select this contig, albeit after giving a warning.
        fancylog(
            f"Warning: Only one contig passes {check_str}. Selecting it.",
            prefix="",
        )
        return passing_contigs[0]

    # Actually start scoring contigs based on their diversity indices.
    contig2score = defaultdict(int)
    # Diversity index columns where there are at least two "passing" contigs
    # that have defined (non-NA) diversity indices
    good_di_cols = []
    for di_col in passing_di.columns:

        # ignore non-diversity-index columns
        if di_col.startswith(DI_PREF):
            di_vals = passing_di[di_col]

            # Ignore diversity index columns where less than two contigs have
            # defined diversity indices, since these don't mean much for our
            # score computation (at least as currently defined)
            finite_di_vals = di_vals[~di_vals.isna()]
            if len(finite_di_vals.index) >= 2:

                scores = normalize_series(finite_di_vals)
                # normalize_series() will return None if the min and max value
                # in finite_di_vals are identical. In this case, we can't
                # generate meaningful scores, so we just move on.
                if scores is not None:
                    good_di_cols.append(di_col)

                    # Update scores.
                    for passing_contig in passing_di.index:
                        if passing_contig in scores:
                            contig2score[passing_contig] += scores[
                                passing_contig
                            ]
                        else:
                            # Penalize this contig for not having a defined
                            # diversity index in this column: give it the max
                            # possible score
                            contig2score[passing_contig] += 1

    if len(good_di_cols) == 0:
        raise SequencingDataError(
            "No diversity index column has at least two contigs that (1) pass "
            f"{check_str} and (2) have defined and distinct diversity indices "
            "in this column."
        )
    # Find the passing contig with the lowest total score:
    # https://stackoverflow.com/a/3282904
    lowest_score_contig = min(passing_contigs, key=contig2score.get)
    return lowest_score_contig


def compute_number_of_mutations_in_full_contig(
    bcf_obj, thresh_type, thresh_vals, contig
):
    """Counts mutations at certain p or r thresholds in a contig.

    This function is designed to be useful for either decoy or target contigs.

    We perform some sanity checking on thresh_vals, but we'll assume that
    the other parameters are well-formed (e.g. thresh_type is either "p" or
    "r", contig is in bcf_obj, etc.)

    Parameters
    ----------
    bcf_obj: pysam.VariantFile
        Object describing a BCF file produced by strainFlye's naive calling.

    thresh_type: str
        Either "p" or "r", depending on which type of mutations were called in
        bcf_obj.

    thresh_vals: range
        Range of values of p or r (depending on thresh_type) at which to
        count mutations in this contig. Must use a step size of 1.
        If a mutation is a mutation for a value of p or r larger than the
        maximum p or r value here (i.e. it's "indisputable"), it will not be
        counted.

    contig: str
        Name of a contig for which mutation rates will be computed.

    Returns
    -------
    num_muts: list
        List with the same length as thresh_vals. The i-th value in this list
        describes the number of called mutations for the i-th threshold in
        thresh_vals.

    Raises
    ------
    ParameterError
        If thresh_vals doesn't use a step size of 1.
        If len(thresh_vals) is zero.
        If either the stop or start of thresh_vals are zero or below.

        (None of these should happen in practice, but you never know...)
    """
    if thresh_vals.step != 1:
        raise ParameterError("thresh_vals must use a step size of 1.")
    if len(thresh_vals) <= 0:
        raise ParameterError("thresh_vals must have a positive length.")
    if thresh_vals.start <= 0 or thresh_vals.stop <= 0:
        raise ParameterError("thresh_vals' start and stop must be positive.")

    # For each threshold value, keep track of how many mutations we've seen at
    # this threshold.
    num_muts = [0] * len(thresh_vals)

    # We can just infer the "indisputable" mutation value from thresh_vals as
    # the value just above the max thresh_vals entry
    high_val = thresh_vals[-1] + 1
    min_val = thresh_vals[0]

    for mut in bcf_obj.fetch(contig):

        # AAD is technically a tuple since it's defined once for every alt
        # allele, but r/n strainflye call only produces max one alt allele per
        # mutation. So it's a tuple with 1 element (at least for now).
        alt_pos = mut.info.get("AAD")[0]
        cov_pos = mut.info.get("MDP")

        if thresh_type == "p":
            max_passing_val = floor((10000 * alt_pos) / cov_pos)
        else:
            max_passing_val = alt_pos

        # Don't count "indisputable" mutations towards mutation rates
        if max_passing_val >= high_val:
            continue

        # NOTE: This is already more optimized than the analysis
        # notebooks, but I think it could still be made faster. Maybe
        # just increment a single value (corresponding to the max
        # passing p/r), and then do everything at the end after seeing
        # all mutations in one pass? Doesn't seem like a huge bottleneck tho.
        num_vals_to_update = max_passing_val - min_val + 1
        for i in range(num_vals_to_update):
            num_muts[i] += 1
    return num_muts


def compute_full_decoy_contig_mut_rates(
    bcf_obj, thresh_type, thresh_vals, decoy_contig, decoy_contig_len
):
    """Computes mutation rates for the entirety of a decoy contig.

    This is designed for the "Full" option, in which we consider every position
    in the contig as a part of the decoy.

    Parameters
    ----------
    bcf_obj: pysam.VariantFile
        Object describing a BCF file produced by strainFlye's naive calling.

    thresh_type: str
        Either "p" or "r", depending on which type of mutations were called in
        bcf_obj.

    thresh_vals: range
        Range of values of p or r (depending on thresh_type) at which to
        compute mutation rates for this contig. Must use a step size of 1.
        If a mutation is a mutation for a value of p or r larger than the
        maximum p or r value here (i.e. it's "indisputable"), it will not be
        included in the mutation rate computation.

    decoy_contig: str
        Decoy contig name.

    decoy_contig_len: int
        Decoy contig sequence length.

    Returns
    -------
    mut_rates: list
        Mutation rates for each threshold value in thresh_vals.
    """
    num_muts = compute_number_of_mutations_in_full_contig(
        bcf_obj, thresh_type, thresh_vals, decoy_contig
    )
    denominator = 3 * decoy_contig_len
    return [n / denominator for n in num_muts]


def compute_target_contig_fdr_curve_info(
    bcf_obj,
    thresh_type,
    thresh_vals,
    target_contig,
    target_contig_len,
    decoy_mut_rates,
):
    """Computes FDR curve information for a given target contig.

    The intent is to create output that can be dropped straight into a TSV
    file, without too much work on the part of the caller.

    Parameters
    ----------
    bcf_obj: pysam.VariantFile
        Object describing a BCF file produced by strainFlye's naive calling.

    thresh_type: str
        Either "p" or "r", depending on which type of mutations were called in
        bcf_obj.

    thresh_vals: range
        Range of values of p or r (depending on thresh_type) at which to
        count mutations in this contig. Must use a step size of 1.
        If a mutation is a mutation for a value of p or r larger than the
        maximum p or r value here (i.e. it's "indisputable"), it will not be
        included in the FDR curve information computation.

    target_contig: str
        Name of the target contig.

    target_contig_len: int
        Target contig sequence length.

    decoy_mut_rates: list
        List of mutation rates (with the same length as thresh_vals) for the
        decoy contig. The context used to compute these mutation rates (Full,
        CP2, Nonsyn, ...) doesn't matter. All we care about is these rates.

    Returns
    -------
    (fdr_line, num_line): (str, str)
        These are both tab-separated lines (suitable for adding to a TSV file
        where the first column is the target contig name and there are
        len(thresh_vals) additional columns).

        The first entry in fdr_line and num_line is the target contig name.
        The remaining entries in fdr_line describe the estimated FDR for this
        target contig at each threshold value (the x-axis on the FDR curves, as
        we currently draw them). The remaining entries in num_line describe the
        number of mutations per megabase for this target contig at each
        threshold value (the y-axis on the FDR curves, as we currently draw
        them).
    """
    num_muts = compute_number_of_mutations_in_full_contig(
        bcf_obj, thresh_type, thresh_vals, target_contig
    )

    # it's a long story. see docs for compute_num_mutations_per_mb() in
    # https://github.com/fedarko/sheepgut/blob/main/notebooks/DemonstratingTargetDecoyApproach.ipynb
    numpermb_coeff = 1000000 / target_contig_len
    denominator = 3 * target_contig_len
    # The 100 is because we convert FDRs from [0, 1] to percentages
    # (technically, FDRs can exceed 100% if the decoy mutation rate > the
    # target mutation rate, but that shouldn't happen too often)
    fdr_coeff = 100 * denominator

    fdr_line = target_contig
    num_line = target_contig
    for i, n in enumerate(num_muts):
        if n == 0:
            # The FDR is undefined if the target's mutation rate is zero
            fdr_val = "NA"
            num_val = "0"
        else:
            fdr_val = str((fdr_coeff * decoy_mut_rates[i]) / n)
            num_val = str(numpermb_coeff * n)
        fdr_line += f"\t{fdr_val}"
        num_line += f"\t{num_val}"
    return fdr_line + "\n", num_line + "\n"


def compute_decoy_contig_mut_rates(
    contigs,
    bcf_obj,
    thresh_type,
    thresh_vals,
    decoy_contig,
    decoy_context,
):
    """Computes mutation rates for a decoy contig at some threshold values.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs in which mutations were
        naively called. We'll only really use this to extract the decoy
        contig's sequence (it's probably easier to load it here then to rely on
        the caller to load it).

    bcf_obj: pysam.VariantFile
        Object describing a BCF file produced by strainFlye's naive calling.

    thresh_type: str
        Either "p" or "r", depending on which type of mutations were called in
        bcf_obj.

    thresh_vals: range
        Range of values of p or r (depending on thresh_type) at which to
        compute mutation rates. Must use a step size of 1.

    decoy_contig: str
        Name of a contig to compute mutation rates for.

    decoy_context: str
        Context-dependent mutation settings to apply in computing the mutation
        rates. One of the following:
         - "Full": Compute mutation rates across the entire contig.
         - "CP2": Only consider positions that are 1) located in a single
                  predicted protein-coding gene and 2) are located in the
                  second codon position of this gene.
         - "Nonsyn": Compute mutation rates based on treating potential
                     nonsynonymous mutations (for positions located in a single
                     predicted protein-coding gene) as a decoy.
         - "Nonsense": Like Nonsyn, but for nonsense mutations.
         - "CP2Nonsyn": Nonsyn, but only for positions in CP2.
         - "CP2Nonsense": Nonsense, but only for positions in CP2.

    Returns
    -------
    decoy_mutation_rates: list
        Has the same dimensions as thresh_vals. The i-th value of
        decoy_mutation_rates corresponds to the mutation rate computed for this
        decoy contig based on naive calling using the i-th threshold value.

    Raises
    ------
    SequencingDataError
        If decoy_contig isn't in the contigs. (This should never happen if
        this function is called from run_estimate(), which should already have
        ensured that this is the case.)
    """
    decoy_seq = fasta_utils.get_single_seq(contigs, decoy_contig)

    if decoy_context == "Full":
        return compute_full_decoy_contig_mut_rates(
            bcf_obj,
            thresh_type,
            thresh_vals,
            decoy_contig,
            len(decoy_seq),
        )
    else:
        raise NotImplementedError(
            'Only the "Full" context is implemented now.'
        )

    # TODO:
    # - If decoy_context != "Full",
    #   - Predict genes in this sequence using prodigal. Save .sco to tempfile.
    # - Fetch mutations aligned to this contig in the BCF file.
    # -


def run_estimate(
    contigs,
    bcf,
    diversity_indices,
    decoy_contig,
    decoy_context,
    high_p,
    high_r,
    decoy_min_length,
    decoy_min_average_coverage,
    output_fdr_info,
    output_num_info,
    fancylog,
):
    """Runs the pipeline for decoy selection and FDR estimation.

    Notably, both high_p and high_r will be defined regardless of if the BCF
    was generated using p- or r-mutations, because (unlike strainFlye call)
    I don't think it's worth splitting this step up into two sub-commands
    by p- or r-mutations. So we'll just ignore one of these values.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs in which mutations were
        naively called.

    bcf: str
        Filepath to a BCF file generated by one of strainFlye call's
        subcommands.

    diversity_indices: str or None
        If a str, this should be a filepath to a TSV file containing diversity
        index info, also generated by one of strainFlye call's subcommands.
        We'll use this information to automatically select a decoy contig.

    decoy_contig: str or None
        If a str, this should be the name of a contig described in the
        BCF file. We'll use this as a decoy contig.

    decoy_context: str
        Context-dependent mutation settings (e.g. Full, CP2, Nonsyn, ...)
        Thankfully, we know this will be one of a set of allowed choices,
        since we use click.Choice() to screen it at the CLI.

    high_p: int
        "Indisputable" threshold for p-mutations (scaled up by 100).

    high_r: int
        "Indisputable" threshold for r-mutations.

    decoy_min_length: int
        If automatically selecting decoy contigs, we'll only consider contigs
        that are at least this long.

    decoy_min_average_coverage: float
        If automatically selecting decoy contigs, we'll only consider contigs
        with average coverages of at least this.

    output_fdr_info: str
        Filepath to which we'll write a TSV file describing estimated FDRs
        for the target contigs. (x-axis values for the FDR curves, at least as
        plotted in the paper.)

    output_num_info: str
        Filepath to which we'll write a TSV file describing numbers of
        mutations per megabase. (y-axis values for the FDR curves.)

    fancylog: function
        Logging function.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        - If the selected decoy contig isn't present in the contigs file.
        - If the high p (or high r) threshold is <= the minimum threshold value
          used in the BCF file.
        - If the BCF file's contigs do not exactly match the FASTA file's
          contigs.

    - parse_bcf() can also raise various errors if the input BCF is malformed.
    - fasta_utils.get_name2len() can also raise errors if the input FASTA file
      of contigs is malformed.
    """
    fancylog("Loading and checking contig information...")

    # get name -> length mapping for the FASTA file; also sanity check it a bit
    contig_name2len = fasta_utils.get_name2len(contigs)

    # Load the BCF file, also
    bcf_obj, thresh_type, thresh_min = parse_bcf(bcf)

    # Figure out which contigs are considered in the BCF file
    # (thankfully, this header can include contigs with no called mutations,
    # which makes my life easier here)
    bcf_contigs = set(bcf_obj.header.contigs)

    # Ensure that the sets of contigs in the BCF file and FASTA match exactly
    # (In theory, we could allow the BCF to be a subset of the FASTA, but...
    # nah, that's too much work and the user should already have an exactly-
    # matching FASTA file around from when they ran "call".)
    fasta_utils.verify_contigs_subset(
        bcf_contigs,
        set(contig_name2len),
        "the BCF file",
        "the FASTA file",
        exact=True,
    )
    # We *could* try to ensure that the diversity index file's contigs, if
    # a diversity index file is specified, are a subset of the BCF -- but
    # no need to do this extra work right now. the main thing that matters IMO
    # is just checking that the selected decoy contig is in the BCF, which is
    # much easier to do later.
    fancylog(
        (
            "The BCF and FASTA files describe the same "
            f"{len(contig_name2len):,} contigs."
        ),
        prefix="",
    )
    fancylog(
        f"Also, the input BCF file contains {thresh_type}-mutations "
        f"(minimum {thresh_type} = {thresh_min:,}).",
        prefix="",
    )

    # Identify decoy contig
    selection_type = check_decoy_selection(diversity_indices, decoy_contig)
    if selection_type == "DI":
        fancylog("Selecting a decoy contig based on the diversity indices...")
        used_decoy_contig = autoselect_decoy(
            diversity_indices,
            decoy_min_length,
            decoy_min_average_coverage,
            fancylog,
        )
        fancylog(
            f"Selected {used_decoy_contig} as the decoy contig.", prefix=""
        )
    else:
        used_decoy_contig = decoy_contig
        fancylog(f"The specified decoy contig is {used_decoy_contig}.")

    # verify that the decoy contig is actually contained in the FASTA file
    if used_decoy_contig not in contig_name2len:
        raise ParameterError(
            f"Selected decoy contig {used_decoy_contig} is not present in "
            f"{contigs}."
        )
    fancylog(
        (
            "Verified that this decoy contig is present in the BCF and FASTA "
            "files."
        ),
        prefix="",
    )

    # if someone chuckles at this, the project is successful
    # that's how it works
    fancylog("(Sorry for doubting you.)", prefix="")

    # Figure out the exact p or r values we'll iterate through -- these will
    # correspond to the columns in our output TSV of FDR estimation (i.e.
    # for each target contig, we'll produce this many FDR estimates).
    #
    # NOTE 1: For the time being, we just go through in increments of 1 (for
    # p, this is 0.01%; for r, this is just 1 read). We could support other
    # "step" values, but that's a lot of work for probably little benefit
    # (unless this ends up being a bottleneck, idk).
    #
    # NOTE 2: We do not produce an estimate for the exact high_p (or high_r)
    # value. THe maximum threshold is that minus the step value (... which is
    # always 1, at least right now).
    fancylog(f"Determining range of values of {thresh_type} to consider...")
    if thresh_type == "p":
        if high_p <= thresh_min:
            raise ParameterError(
                f"--high-p = {high_p:,} must be larger than the minimum p "
                f"used in the BCF ({thresh_min:,})."
            )
        thresh_high = high_p
    else:
        if high_r <= thresh_min:
            raise ParameterError(
                f"--high-r = {high_r:,} must be larger than the minimum r "
                f"used in the BCF ({thresh_min:,})."
            )
        thresh_high = high_r
    thresh_max = thresh_high - 1
    thresh_vals = range(thresh_min, thresh_high)
    fancylog(
        (
            f"We'll consider {len(thresh_vals):,} values of {thresh_type}: "
            f"from {thresh_min:,} to {thresh_max:,}."
        ),
        prefix="",
    )

    # Extra clarification about indisputable mutations. I want to avoid taking
    # people by surprise with this.
    fancylog(
        (
            f"{thresh_type}-mutations for {thresh_type} \u2265 "
            f'{thresh_high:,} will be considered "indisputable."'
        ),
        prefix="",
    )
    fancylog(
        (
            'These "indisputable" mutations won\'t be included in the FDR '
            "estimation results."
        ),
        prefix="",
    )

    fancylog(
        f"Computing mutation rates for {used_decoy_contig} at these threshold "
        "values..."
    )
    # For each value in thresh_vals, compute the decoy genome's mutation rate.
    decoy_mut_rates = compute_decoy_contig_mut_rates(
        contigs,
        bcf_obj,
        thresh_type,
        thresh_vals,
        used_decoy_contig,
        decoy_context,
    )
    fancylog("Done.", prefix="")

    fancylog(
        "Computing mutation rates and FDR estimates for the "
        f"{len(contig_name2len) - 1:,} target contigs..."
    )

    # Create the header for the TSV files
    tsv_header = "Contig"
    for val in thresh_vals:
        tsv_header += f"\t{thresh_type}{val}"

    # ... and write it out. The header is the same for the FDR estimate file
    # and for the (# of mutations per megabase) file.
    for tsv_fp in (output_fdr_info, output_num_info):
        with open(tsv_fp, "w") as fdr_file:
            fdr_file.write(f"{tsv_header}\n")

    # Compute FDR estimates for each target contig.
    # This is analogous to the "Full" context-dependent option for the decoy
    # genome comptuation, so we can reuse a lot of code from that.
    for target_contig in bcf_contigs - {used_decoy_contig}:
        fdr_line, numpermb_line = compute_target_contig_fdr_curve_info(
            bcf_obj,
            thresh_type,
            thresh_vals,
            target_contig,
            contig_name2len[target_contig],
            decoy_mut_rates=decoy_mut_rates,
        )

        # TODO chunk outputs?
        with open(output_fdr_info, "a") as fdr_file:
            fdr_file.write(fdr_line)
        with open(output_num_info, "a") as num_file:
            num_file.write(numpermb_line)

    fancylog("Done.", prefix="")

    # TODO: Verify that the decoy contig has a nonzero mutation rate?
    # If not, that's problematic, because
    # then we'd estimate the FDR as zero for every target contig. That
    # shouldn't happen most of the time, anyway. Maybe add an option to limit
    # auto-selection to just contigs with mutations? Hm, but that would be
    # annoying to implement, and users can always manually set a decoy contig.)


def get_optimal_threshold_values(fi, fdr):
    """Returns the optimal values of p or r for each contig's mutation calls.

    Parameters
    ----------
    fi: pd.DataFrame
        FDR information produced by "strainFlye fdr estimate". The indices
        (rows) correspond to contigs; the columns correspond to threshold
        values of p or r (sorted in ascending order from left to right). Cells
        indicate the estimated FDR for this contig at this threshold value.

    fdr: float
        FDR at which (non-indisputable) mutation calls for each contig will be
        fixed.

    Returns
    -------
    pd.Series
        Has the same index as fi (so, this has one entry per contig). Each
        contig's entry will be the smallest threshold value (column name) at
        which this contig's FDR is less than or equal to the specified fdr.
        If no such threshold value exists (i.e. all estimated FDRs for this
        contig are greater than the specified fdr), the contig's entry in this
        series will be None.
    """
    # For a given contig's FDR curve, we define four cases regarding how this
    # curve intersects with a fixed FDR (represented as a vertical line).
    # The goal is to, for each contig, identify the "optimal" value of p or r
    # that maximizes the number of mutations / mb we see while keeping the
    # estimated FDR below the fixed value.
    #
    # How do we identify this "optimal" value? By looking at the FDR curve.
    # We define four "cases" of how a FDR curve can look, relative to the fixed
    # FDR (which can be thought of as an infinite vertical line).
    #     |
    # Case 1. The curve crosses the line at exactly one point.
    #     |   This is easy to handle -- select the value of p or r just below
    #     |   the intersection for this contig.
    # ^  _|_____
    # | / |
    # |/  |
    # +---|------->
    #     |
    # Case 2. The curve crosses the line more than one point.
    #     |   (This is an unfortunate possibility, because we do not have the
    #     |   guarantee that FDR estimates increase monotonically with lower
    #     |   values of p or r.) In this case, we select the value of p or r
    #     |   just below the intersection for this contig with the highest
    #     |   # mutations / mb (in the plot below, at the top intersection).
    #     |
    #     |   Notably, this will always correspond to the lowest value of p or
    #     |   r that is <= the fixed FDR: lower values of p or r also "include"
    #     |   the p- or r-mutations called from higher values of p or r, so the
    #     |   lowest "passing" threshold value will also result in the highest
    #     |   number of mutations per megabase in the target contig.
    # ^  _|_____
    # | / |
    # |/  |
    # |\__|__
    # |  _|_/
    # | / |
    # |/  |
    # +---|------->
    #     |
    # Case 3. The curve crosses the line at zero points, and all estimated FDRs
    #     |   are LOWER than the fixed FDR. In this case, just select the
    #     |   lowest value of p or r used.
    #     |
    # ^ | |
    # | | |
    # |/  |
    # +---|------->
    #     |
    # Case 4. The curve crosses the line at zero points, and all estimated FDRs
    #     |   are HIGHER than the fixed FDR. In this case, our hands are tied;
    #     |   don't select any value of p or r. We can't call any
    #     |   non-indisputable mutations for this contig, at least not at this
    #     |   fixed FDR.
    #     |
    # ^   |    /
    # |   |   /
    # |   |  |
    # +---|------->

    # Convert the FDR information to a binary matrix:
    # True  means this FDR is <= the fixed FDR value
    # False means this FDR is >  the fixed FDR value
    lte_fdr = fi <= fdr

    def get_optimal_threshold(fdr_df_row):
        # We implicitly go through from lower to higher threshold values -- so
        # we'll select the lowest threshold value that is <= the fixed FDR.
        # This automatically accounts for Cases 1, 2, and 3 above.
        for ci, val in enumerate(fdr_df_row):
            if val:
                return fdr_df_row.index[ci]
        # If we make it here, then none of the threshold values were <= the
        # fixed FDR. So this is a "Case 4" situation -- return None.
        return None

    # We can select each optimal threshold using DataFrame.apply(). This is
    # faster than naive looping through the DataFrame, but it could still
    # probably be sped up (although I'm not sure how exactly we could use
    # vectorization here). See
    # https://web.archive.org/web/20181106230656/https://engineering.upside.com/a-beginners-guide-to-optimizing-pandas-code-for-speed-c09ef2c6a4d6
    # for details.
    return lte_fdr.apply(get_optimal_threshold, axis=1)


def write_filtered_bcf(
    in_bcf_obj, out_bcf_fp, optimal_thresh_vals, thresh_high
):
    """Writes out a filtered BCF based on optimal p or r values."""
    pass


def load_and_sanity_check_fdr_file(fdr_info, thresh_type):
    """Loads and sanity-checks a FDR TSV file.

    This sanity-checks the structure of the file -- it doesn't say anything
    about the contigs all matching up with a BCF file, for example. The main
    goal is just checking that this looks like the sort of TSV file that the
    estimate command would have generated.

    Parameters
    ----------
    fdr_info: str
        Filepath to a FDR TSV file.

    thresh_type: str
        Either "p" or "r", depending on which type of mutations were called.

    Returns
    -------
    fi: pd.DataFrame
        DataFrame containing the information stored in the TSV file. Indices
        (rows) correspond to contigs; columns correspond to threshold values
        of p or r (sorted in ascending order from left to right). Cells
        indicate the estimated FDR for this contig at this threshold value.

    Raises
    ------
    ParameterError
        If the TSV file seems malformed.
    ParserError
        Can be raised by pd.read_csv() if the file seems malformed.
        We don't attempt to catch this sort of error.
    """
    fi = pd.read_csv(fdr_info, sep="\t", index_col=0)
    # error prefix -- saves us some typing...
    ep = f"TSV file {fdr_info} seems malformed"
    if fi.index.name != "Contig":
        raise ParameterError(f'{ep}: no "Contig" header?')
    if len(fi.index) < 1:
        raise ParameterError(f"{ep}: no contigs described?")
    if len(fi.columns) < 1:
        raise ParameterError(f"{ep}: no threshold values described?")

    prev_col_val = 0
    for col in fi.columns:
        if col[0] != thresh_type:
            raise ParameterError(
                f"{ep}: columns should start with {thresh_type}."
            )
        col_val = int(col[1:])
        if col_val <= prev_col_val:
            raise ParameterError(
                f"{ep}: values of {thresh_type} should increase from left to "
                "right."
            )
        # We could eventually support different step sizes, but not yet.
        # Actually I think this might be overzealous -- I don't think we do
        # anything in "fix" that requires known step sizes -- but whatever, I
        # wanna be careful.
        if col_val != prev_col_val + 1:
            raise ParameterError(
                f"{ep}: values of {thresh_type} should only increase in steps "
                "of 1."
            )
    if (fi < 0).any.any():
        raise ParameterError(f"{ep}: FDR estimates must be nonnegative or NA.")

    return fi


def run_fix(bcf, fdr_info, fdr, output_bcf, fancylog):
    """Runs the pipeline for FDR fixing.

    Parameters
    ----------
    bcf: str
        Filepath to a BCF file generated by one of strainFlye call's
        subcommands.

    fdr_info: str
        Filepath to a TSV file describing estimated FDRs for the target
        contigs.

    fdr: int
        False Discovery Rate (FDR) to fix mutation calls at. Scaled up by 100,
        like values of p are elsewhere -- so fdr = 100 indicates an FDR of 1%.

    output_bcf: str
        Filepath to which we will write an (indexed) BCF file. This will
        contain a subset of the mutations described in the input BCF file.

    fancylog: function
        Logging function.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        - If the set of contigs in the FDR TSV file is not a subset of those
          in the BCF file.
        - If there is not exactly one contig that is present in the BCF file
          but not in the FDR TSV file.

    - parse_bcf() can also raise various errors if the input BCF is malformed.
    """
    fancylog("Loading and checking BCF and TSV files...")
    # Like in run_estimate(): Load the BCF file and figure out what contigs it
    # describes
    bcf_obj, thresh_type, thresh_min = parse_bcf(bcf)
    bcf_contigs = set(bcf_obj.header.contigs)

    # Load the estimated FDR file.
    fi = load_and_sanity_check_fdr_file(fdr_info, thresh_type)
    tsv_desc = "the FDR TSV file"
    # Ensure that the contigs described in these TSV files (tsv_contigs)
    # are all described in the BCF file. We don't check for an exact match,
    # because the decoy contig will be missing.
    tsv_contigs = set(fi.index)
    fasta_utils.verify_contigs_subset(
        tsv_contigs, bcf_contigs, tsv_desc, "the BCF file"
    )
    absent_contigs = bcf_contigs - tsv_contigs
    if len(absent_contigs) != 1:
        # Fun thing about this error message: we can say "contigs" because this
        # error fires if and only if len(absent_contigs) is not 1 :D
        raise ParameterError(
            "Exactly one contig in the BCF file (the decoy contig) should be "
            f"missing from {tsv_desc}. However, {len(absent_contigs):,} "
            "contigs are missing!"
        )

    # If we've made it here, we know there is exactly one contig in the BCF but
    # not in the FDR TSV file. This must be the decoy contig.
    decoy_contig = absent_contigs.pop()
    fancylog(
        f"Looks good so far; decoy contig seems to be {decoy_contig}.",
        prefix="",
    )

    # Figure out the "indisputable" mutation cutoff.
    thresh_max = int(fi.columns[-1][1:])
    thresh_high = thresh_max + 1
    fancylog(
        (
            'Looks like the cutoff for "indisputable" mutations was '
            f"{thresh_type} = {thresh_high:,}."
        ),
        prefix="",
    )
    fancylog(
        (
            "All mutations passing this cutoff will be included in the "
            "output BCF file."
        ),
        prefix="",
    )

    fancylog(f"Finding optimal values of {thresh_type} for each contig...")
    # We can now begin this in earnest. Figure out the "optimal" value of p or
    # r for each contig, based on the estimated FDR information.
    optimal_thresh_vals = get_optimal_threshold_values(fi, fdr)
    fancylog("Done.", prefix="")

    fancylog("Writing a filtered BCF file...")
    # Filter mutations for each contig to those that pass these thresholds (in
    # addition to indisputable mutations, using thresh_high from above).
    write_filtered_bcf(bcf_obj, output_bcf, optimal_thresh_vals, thresh_high)
    fancylog("Done.", prefix="")

    # ... Then, just make sure to index the output BCF file, and we're done!
    call_utils.index_bcf(output_bcf)
