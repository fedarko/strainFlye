# Utilities for strainFlye fdr.


import re
import pysam
import pandas as pd
from .errors import ParameterError, SequencingDataError


def parse_vcf(vcf):
    """Opens a VCF file, and does sanity checking and p vs. r sniffing on it.

    Thankfully, pysam has the ability to read VCF files, so the main thing
    we do here is checking that the meta-information of the VCF file seems
    seems kosher (i.e. was produced by strainFlye).

    Parameters
    ----------
    vcf: str
        Path to a VCF file produced by one of "strainFlye call"'s subcommands.

    Returns
    -------
    (f, thresh_type, thresh_min): (pysam.VariantFile, str, int)
        f: Object describing the input VCF file.
        thresh_type: either "p" or "r", depending on what type of mutation
                     calling was done to produce this VCF file.
        thresh_min: the minimum value of p or r used in mutation calling to
                    produce this VCF file. Formatted the same as used in the
                    file (so values of p will still be scaled up by 100).

    Raises
    ------
    FileNotFoundError
        If vcf doesn't exist (raised by pysam).

    ValueError
        If vcf doesn't look like a VCF/BCF file (raised by pysam).

    ParameterError
        If vcf does not have exactly one "strainFlye threshold filter header,"
        which is a term I made up just now. Basically, we rely on there
        existing a single line in the VCF file's header that goes like

            ##FILTER=<ID=strainflye_minT_MMMM,...>

        ... where T corresponds to the type of mutation calling done (p or r)
        and MMMM (variable number of digits) indicates the minimum value of p
        or r used in mutation calling.

        If there isn't exactly one of these lines, then we will be very
        confused! Hence why we raise an error.
    """
    # this will fail with a FileNotFoundError if "vcf" doesn't point to an
    # existing file (although we shouldn't need to worry about that much b/c
    # click should've already checked that this file exists); and it'll fail
    # with a ValueError if it points to a file but this file doesn't look like
    # a VCF/BCF file (it shouldn't be BCF unless the user converted it to BCF
    # themselves, I guess, but that's fine)
    f = pysam.VariantFile(vcf)

    # Now that this at least seems like a VCF/BCF file, make sure it's from
    # strainFlye, and figure out whether it's from p- or r-mutation calling...
    thresh_type = None
    thresh_min = None

    # pysam seems to add an extra filter labelled PASS to the parsed VCF file,
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
                    f"VCF file {vcf} has multiple strainFlye threshold filter "
                    "headers."
                )
            thresh_type = filter_match.group(1)
            thresh_min = int(filter_match.group(2))

    # If we never updated these variables, we never saw a strainFlye filter
    # header -- probably this VCF isn't from strainFlye.
    if thresh_type is None or thresh_min is None:
        raise ParameterError(
            f"VCF file {vcf} doesn't seem to be from strainFlye: no threshold "
            "filter headers."
        )

    # Let's be extra paranoid and verify that this VCF has MDP (coverage based
    # on (mis)matches) and AAD (alternate nucleotide coverage) fields. (It
    # should, because we know at this point that strainFlye generated it, but
    # you never know...)
    info_ids = f.header.info.keys()
    if "MDP" not in info_ids or "AAD" not in info_ids:
        raise ParameterError(
            f"VCF file {vcf} needs to have MDP and AAD info fields."
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
        VCF file.

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


def autoselect_decoy(diversity_indices, min_len, min_avg_cov):
    """Attempts to select a good decoy contig based on diversity index data.

    There are lots of ways to implement this, so here we just stick with
    something simple. Filter to all contigs whose lengths and average coverages
    meet some thresholds, then determine the five lowest-diversity-index
    contigs across all diversity index columns provided in the file. Select as
    the decoy the contig that appears most frequently in these lists of five
    contigs, breaking ties based on lowest diversity index.

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
        - If none of the "passing" contigs in the diversity index file have
          defined diversity indices (should never happen unless the
          diversity index values used to generate this file are weird).
    """
    di = pd.read_csv(diversity_indices, sep="\t", index_col=0)
    if len(di.index) < 2:
        raise ParameterError(
            "Diversity indices file describes less than two contigs."
        )
    if "Length" not in di.columns and "AverageCoverage" not in di.columns:
        raise ParameterError(
            "Length and AverageCoverage columns are not contained in the "
            "diversity indices file."
        )
    # Filter to contigs that pass both the length and coverage thresholds.
    # https://stackoverflow.com/a/13616382
    valid_di = di[
        (di["Length"] >= min_len) & (di["AverageCoverage"] >= min_avg_cov)
    ]
    if len(valid_di.index) == 0:
        raise SequencingDataError(
            f"No contigs pass the min length \u2265 {min_len:,} and min "
            f"average cov \u2265 {min_avg_cov:,}x checks."
        )
    # TODO actually do stuff


def run_estimate(
    vcf,
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

    Notably, both high_p and high_r will be defined regardless of if the VCF
    was generated using p- or r-mutations, because (unlike strainFlye call)
    I don't think it's worth splitting this step up into two sub-commands
    by p- or r-mutations. So we'll just ignore one of these values.

    Parameters
    ----------
    vcf: str
        Filepath to a VCF file generated by one of strainFlye call's
        subcommands.

    diversity_indices: str or None
        If a str, this should be a filepath to a TSV file containing diversity
        index info, also generated by one of strainFlye call's subcommands.
        We'll use this information to automatically select a decoy contig.

    decoy_contig: str or None
        If a str, this should be the name of a contig described in the
        VCF file. We'll use this as a decoy contig.

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
    # Are we using p or r? And what's the minimum p or r that was used?
    vcf_obj, thresh_type, thresh_min = parse_vcf(vcf)
    fancylog(
        f"Input VCF file contains {thresh_type}-mutations (minimum "
        f"{thresh_type} = {thresh_min:,})."
    )

    # Identify decoy contig
    selection_type = check_decoy_selection(diversity_indices, decoy_contig)
    if selection_type == "DI":
        fancylog("Selecting a decoy contig based on the diversity indices...")
        # Automatically select a decoy contig from the diversity indices
        used_decoy_contig = autoselect_decoy(
            diversity_indices, decoy_min_length, decoy_min_average_coverage
        )
        fancylog(f"Using {used_decoy_contig} as the decoy contig.", prefix="")
    else:
        used_decoy_contig = decoy_contig
        fancylog(f"The specified decoy contig is {used_decoy_contig}.")

    # Verify that the decoy contig is actually contained in the VCF. (If not,
    # it could still be in the contigs, but it could just not have any called
    # mutations -- but that is problematic, because then we'd estimate the FDR
    # as zero for every target contig. That shouldn't happen most of the time,
    # anyway.)

    # Figure out range of p or r to use. Create a list, threshold_vals.
    # For each value in threshold_vals, compute the decoy genome's mutation
    # rate. Save this to a list, decoy_mut_rates -- this will have the same
    # dimensions as threshold_vals.

    # For each target genome...
    # - For each value in threshold_vals...
    #   - Compute the mutation rate for this target genome at this
    #     threshold value.
    #   - Compute the FDR estimate for this pair of (target, threshold).
    #     Save to a list of target_fdr_ests, which has the same dimensions
    #     as threshold_vals.
    # - Write out a new row to the FDR estimate file describing
    #   target_fdr_ests.


def run_fix(vcf, fdr_info, fdr, output_vcf, fancylog):
    pass
