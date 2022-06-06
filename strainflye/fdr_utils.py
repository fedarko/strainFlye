# Utilities for strainFlye fdr.


import re
import pysam
import pandas as pd
from collections import defaultdict
from strainflye import fasta_utils
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


def compute_decoy_mut_rates(
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

    thresh_vals: list
        List of values of p or r (depending on thresh_type) at which to
        compute mutation rates.

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

    Raises
    ------
    ParameterError
        If decoy_contig isn't in the contigs
    """
    # - Load contig sequence from contigs
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
        for the target contigs.

    fancylog: function
        Logging function.

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
        thresh_max = high_p - 1
    else:
        if high_r <= thresh_min:
            raise ParameterError(
                f"--high-r = {high_r:,} must be larger than the minimum r "
                f"used in the BCF ({thresh_min:,})."
            )
        thresh_max = high_r - 1
    # ... yeah, we could have just set thresh_max to high_p or high_r without
    # the - 1 and then used it as the top of this range, which is technically
    # more efficient I guess, but we will need thresh_max (as it currently is,
    # with the - 1) later and this way is clearer to read I guess ._.
    thresh_vals = range(thresh_min, thresh_max + 1)
    fancylog(
        (
            f"We'll consider {len(thresh_vals):,} values of {thresh_type}: "
            f"from {thresh_min:,} to {thresh_max:,}."
        ),
        prefix="",
    )

    # For each value in thresh_vals, compute the decoy genome's mutation
    # rate. Will need to predict genes using prodigal first, if decoy_context
    # isn't "Full".
    #
    # Save these to a list, decoy_mut_rates -- this will have the same
    # dimensions as thresh_vals.
    # decoy_mut_rates = compute_decoy_mut_rates(
    #     contigs,
    #     bcf_obj,
    #     thresh_type,
    #     thresh_vals,
    #     used_decoy_contig,
    #     decoy_context,
    # )

    # TODO: Verify that the decoy contig has a nonzero mutation rate?
    # If not, that's problematic, because
    # then we'd estimate the FDR as zero for every target contig. That
    # shouldn't happen most of the time, anyway. Maybe add an option to limit
    # auto-selection to just contigs with mutations? Hm, but that would be
    # annoying to implement, and users can always manually set a decoy contig.)

    # TODO: For each target genome...
    # - For each value in thresh_vals...
    #   - Compute the mutation rate for this target genome at this
    #     threshold value. The entire target genome, not just the dctx stuff.
    #   - Compute the FDR estimate for this pair of (target, threshold).
    #     Save to a list of target_fdr_ests, which has the same dimensions
    #     as thresh_vals.
    # - Write out a new row to the FDR estimate file describing
    #   target_fdr_ests.


def run_fix(bcf, fdr_info, fdr, output_bcf, fancylog):
    pass
