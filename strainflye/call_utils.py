# Utilities for strainFlye's call step.


import os
import time
import pysam
import pysamstats
from . import cli_utils, config
from .errors import SequencingDataError, ParameterError
from strainflye import __version__


def get_r_increments(min_r, max_r, delta_r, fancylog):
    """Given a min, max, and delta r, computes a list of r values.

    This should be a simple problem if the delta is 1 and the min and max are
    specified sanely, but... yeah, we gotta be defensive.

    Parameters
    ----------
    min_r: int >= 1
        Minimum value of r to use.

    max_r: int >= 1
        Maximum value of r to use.

    delta_r: int >= 1
        Increment r by this much, starting at min_r.

    fancylog: function
        Logging function.

    Returns
    -------
    r_vals: list of ints >= 1
        List of r values beginning at min_r and increasing by delta_r until we
        reach max_r. Note that this will not include max_r if
        (max_r - min_r) % delta_r != 0 -- we'll warn the user via fancylog()
        if this happens.

    Raises
    ------
    ParameterError
        If min_r >= max_r.

    ValueError
        - If something goes seriously wrong when generating r_vals.
          This should never happen (knock on wood).
    """
    if min_r >= max_r:
        raise ParameterError("Minimum r must be less than maximum r.")

    r_vals = list(range(min_r, max_r + 1, delta_r))

    if len(r_vals) == 1:
        fancylog(f"Computing r-mutations for 1 value of r: {r_vals[0]:,}.")
    else:
        fancylog(f"Computing r-mutations for {len(r_vals):,} values of r.")

    # If r is not 1, it's possible for us to "miss" the maximum r value.
    # Warn the user about this, but don't throw an error.
    if r_vals[-1] != max_r:
        divisibility = (max_r - min_r) % delta_r
        # This should never happen, but check it anyway
        if divisibility == 0:
            raise ValueError("Max r was excluded without cause?")
        fancylog(
            f"Warning: --max-r = {max_r:,} will not be included in the "
            "r-values used. This is due to --max-r minus --min-r not "
            f"being divisible by --delta-r: {max_r:,} - {min_r:,} = "
            f"{max_r - min_r:,}, and {max_r - min_r:,} mod {delta_r:,} = "
            f"{divisibility:,}."
        )

    # Extra sanity check (this should definitely not happen)
    if r_vals[0] != min_r:
        raise ValueError("Min r was excluded without cause?")

    return r_vals


def get_p_increments(min_p, max_p, delta_p, fancylog):
    """Given a min, max, and delta p, computes a list of p values.

    We bypass floating-point issues by just keeping the three parameters as
    integers, and then dividing them by 100 later. Makes my life easier,
    although it makes the UI a bit more complicated.

    This also handles conversion back to normal percentages, so the output
    values of p are in the range (0, 50] rather than (0, 5000].

    NOTE: I'm well aware that this function is annoyingly similar to
    get_r_increments(), and the shared code should probably be extracted
    somehow eventually to avoid duplication.

    Parameters
    ----------
    min_p: float in (0, 5000)
        Minimum value of p to use.

    max_p: float in (0, 5000]
        Maximum value of p to use.

    delta_p: float in (0, 5000)
        Increment p by this much, starting at min_p.

    fancylog: function
        Logging function.

    Returns
    -------
    p_vals: list of floats in (0, 50]
        List of p values beginning at min_p and increasing by delta_p until we
        reach max_p. Note that this can not include max_p.

    Raises
    ------
    ParameterError
        If min_p >= max_p.

    ValueError
        - If something goes seriously wrong when generating p_vals.
          This should never happen (knock on wood).
    """
    if min_p >= max_p:
        raise ParameterError("Minimum p must be less than maximum p.")

    p_vals = [p for p in range(min_p, max_p + 1, delta_p)]

    if len(p_vals) == 1:
        fancylog(
            f"Computing p-mutations for 1 value of p: {p_vals[0] / 100:.2f}%."
        )
    else:
        fancylog(f"Computing p-mutations for {len(p_vals):,} values of p.")

    # If p is not 1, it's possible for us to "miss" the maximum p value.
    # Warn the user about this, but don't throw an error.
    #
    # Note that this check is a reason why we delay dividing the p_vals by 100
    # until now (although I guess we could restructure things to do the
    # division immediately; whatever, this probs won't be a bottleneck).
    if p_vals[-1] != max_p:
        divisibility = (max_p - min_p) % delta_p
        # This should never happen, but check it anyway
        if divisibility == 0:
            raise ValueError("Max p was excluded without cause?")
        fancylog(
            f"Warning: --max-p = {max_p:,} (aka {max_p / 100:.2f}%) will "
            "not be included in the values of p used. "
            "This is due to --max-p minus --min-p not "
            f"being divisible by --delta-p: {max_p:,} - {min_p:,} = "
            f"{max_p - min_p:,}, and {max_p - min_p:,} mod {delta_p:,} = "
            f"{divisibility:,}."
        )

    # Extra sanity check (this should definitely not happen)
    if p_vals[0] != min_p:
        raise ValueError("Min p was excluded without cause?")

    return [p / 100 for p in p_vals]


def get_alt_pos_info(rec):
    """Returns info about the second-most-common nucleotide at a position.

    This nucleotide will usually differ from the reference nucleotide, but it
    may be the reference (i.e. at positions where the reference disagrees with
    the alignment's "consensus").

    This breaks ties arbitrarily.

    Parameters
    ----------
    rec: dict
        pysamstats record for a given position in an alignment produced
        by stat_variation().

    Returns
    -------
    (cov, alt nt freq, alt nt, ref nt freq, ref nt): (int, int, str, int, str)
        Describes the first- and second-most-common nucleotides at a position.

        The first entry in this tuple is the (mis)match coverage at this
        position. This is an integer defined as the sum of A, C, G, T
        nucleotides at this position (note that this excludes degenerate
        nucleotides like N -- we could change this in the future if that'd be
        useful, I suppose). Note that this coverage could be zero, if no reads
        are aligned to this specific position.

        The second entry is the raw frequency of the second-most-common
        nucleotide at this position: this will be an integer greater than or
        equal to 0. This is also referred to in the paper, etc. as alt(pos).

        The third entry is just the alternate nucleotide (one of A, C, G, T),
        represented as a string.

        The fourth and fifth entries are analogous to the second and third, but
        this time correspond to the first-most-common nucleotide at this
        position.

    References
    ----------
    I mean, I wrote this code, but it's copied (and modified) from
    https://github.com/fedarko/sheepgut/blob/main/notebooks/pleuk_copied_code.py
    (which is in turn copied from https://github.com/fedarko/pleuk, but that
    project is still in the oven so it's private for now).
    """
    cov = rec["A"] + rec["C"] + rec["G"] + rec["T"]

    ordered_nts = sorted("ACGT", key=rec.get)

    # The literal nucleotide used in the numerator of freq(pos): one of A, C,
    # G, T
    alt_nt = ordered_nts[-2]

    # The raw frequency (in counts) of alt_nt. An integer >= 0.
    alt_freq = rec[alt_nt]

    # Replicate this stuff for the "reference" nucleotide.
    # We'll implicitly treat the reference nucleotide as the most common (aka
    # "consensus") nucleotide at this position, even at "unreasonable"
    # positions where the nucleotide located here on the contig disagrees
    # with the consensus
    ref_nt = ordered_nts[-1]

    ref_nt_freq = rec[ref_nt]

    return (cov, alt_freq, alt_nt, ref_nt_freq, ref_nt)


def p2filter(p):
    return f"p{p:.2f}"


def r2filter(r):
    return f"r{r}"


def run(
    contigs,
    bam,
    output_vcf,
    fancylog,
    verbose,
    p_vals=[],
    r_vals=[],
    min_alt_pos=None,
):
    """Launches the process of naive p- or r-mutation calling.

    In the user-facing code in the CLI, I try to consistently say "naive" with
    the two dots, but here I don't bother because users don't see this part of
    the codebase. Here we can, like, kick off our shoes and lie down on the
    couch and eat Cheetos, because nobody is going to read this docstring (and
    if you are reading it then I'm very surprised and also hi, how's it going?)

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs in which to call mutations.

    bam: str
        Filepath to a BAM file mapping reads to contigs.

    output_vcf: str
        Filepath to which a VCF file describing called mutations will be
        written.

    fancylog: function
        Logging function.

    verbose: bool
        Log extra info about individual contigs.

    p_vals: list of floats in (0, 50]
        List of p-mutation parameters for which to call p-mutations.

    r_vals: list of ints > 0
        List of r-mutation parameters for which to call r-mutations.

    min_alt_pos: int >= 1 or None
        During p-mutation calling, the second-most-common aligned nucleotide's
        frequency must be at least this. Not used when calling r-mutations.

    Returns
    -------
    None

    Raises
    ------
    cli_utils.ParameterError
        If something went wrong with the parameters. We don't check (at least
        within this function) that they fall into their respective allowed
        ranges (since Click should have already enforced that thanks to the
        range parameter settings), but we do check that only p or only r is
        specified -- exactly one of these parameters must be given.
    """
    using_p = len(p_vals) > 0
    using_r = len(r_vals) > 0

    if using_p and using_r:
        raise cli_utils.ParameterError(
            "p and r can't be specified at the same time. Please choose one."
        )
    elif not using_p and not using_r:
        raise cli_utils.ParameterError("Either p or r needs to be specified.")

    filter_header = ""

    if using_p:
        param_name = "p"
        call_str = (
            f"p-mutation calling for p from {p_vals[0]:.2f}% to "
            f"{p_vals[-1]:.2f}%"
        )
        for p in p_vals:
            filter_header += (
                f'##FILTER=<ID={p2filter(p)},Description="p = {p:.2f}%">\n'
            )
    else:
        param_name = "r"
        call_str = (
            f"r-mutation calling for r from {r_vals[0]:,} to {r_vals[-1]:,}"
        )
        for r in r_vals:
            filter_header += (
                f'##FILTER=<ID={r2filter(r)},Description="r = {r:,}">\n'
            )

    fancylog(f"Running {call_str}.", prefix="")

    with open(output_vcf, "w") as vcf_file:
        # Header info gleaned by reading over the VCF 4.2 docs
        # (https://samtools.github.io/hts-specs/VCFv4.2.pdf) and copying how
        # LoFreq organizes their header

        vcf_file.write(
            "##fileformat=VCFv4.3\n"
            f"##fileDate={time.strftime('%Y%m%d')}\n"
            f'##source="strainFlye v{__version__}: {call_str}"\n'
            f"{filter_header}"
            f"##reference={os.path.abspath(contigs)}\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )

    bf = pysam.AlignmentFile(bam, "rb")
    for si, seq in enumerate(bf.references, 1):
        if verbose:
            pct = 100 * (si / bf.nreferences)
            fancylog(
                f"On contig {seq} ({si:,} / {bf.nreferences:,}) ({pct:.2f}%).",
                prefix="",
            )
        num_any_muts = 0
        vcf_text = ""
        for pos, rec in enumerate(
            pysamstats.stat_variation(
                bf,
                chrom=seq,
                fafile=contigs,
                pad=True,
                max_depth=config.MAX_DEPTH_PYSAM,
            ),
            1,
        ):
            # Sanity checks
            if rec["N"] > 0:
                raise SequencingDataError(
                    "Alignments including degenerate nucleotides (e.g. N) are "
                    "not supported"
                )

            rpos = rec["pos"] + 1
            if rpos != pos:
                raise ValueError(
                    f"Found discontinuity in traversal: {pos:,}-th pos, but "
                    f"rec['pos'] + 1 is {rpos:,}"
                )

            cov, alt_freq, alt_nt, ref_freq, ref_nt = get_alt_pos_info(rec)

            if using_p:
                filter_info, any_pass = call_p_mutations(
                    alt_freq, cov, p_vals, min_alt_pos
                )
            else:
                filter_info, any_pass = call_r_mutations(alt_freq, r_vals)

            if any_pass:
                # 1. CHROM = seq (aka contig name)
                #
                # 2. POS = pos (position on this contig)
                #
                # 3. ID = . (this can be a unique identifier for this variant,
                #            but... we don't have that sort of information)
                #
                # 4. REF = ref_nt (consensus nucleotide -- this can disagree
                #                  with the actual nucleotide at this position
                #                  on the contig in the case of "unreasonable"
                #                  positions)
                #
                # 5. ALT = alt_nt (second-most-common nucleotide at this
                #                  position, for which we are calling a
                #                  mutation; note that we currently just call
                #                  at most one mutation per position, although
                #                  this could in theory be generalized to work
                #                  with multiallelic mutations)
                #
                # 6. QUAL = . (In lieu of providing a single probability here,
                #              we use the FDR estimation stuff described in the
                #              paper)
                #
                # 7. FILTER = filter info about p or r thresholds which this
                #             mutation didn't pass
                #
                # 8. INFO = . (Maybe I'll add extra stuff here later)
                vcf_text += (
                    f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t{filter_info}\t."
                    "\n"
                )
                num_any_muts += 1

        if num_any_muts > 0:
            with open(output_vcf, "a") as vcf_file:
                vcf_file.write(vcf_text)

        if verbose:
            fancylog(
                (
                    f"Called {num_any_muts:,} mutation(s) (for any setting of "
                    f"{param_name}) in contig {seq}."
                ),
                prefix="",
            )


def is_position_rare_direct(alt_pos, cov):
    """Determines if a p-mutated position is a "rare" mutation.

    Parameters
    ----------
    alt_pos: int

    cov: int

    Returns
    -------
    bool
        True if this position is "rare" (aka its mutation frequency is less
        than config.HIGH_FREQUENCY_MIN_PCT), False otherwise.
    """
    lhs = 100 * alt_pos
    rhs_upper = config.HIGH_FREQUENCY_MIN_PCT * cov
    return lhs < rhs_upper


def convert_calling_output(filter_info, any_pass):
    if any_pass:
        if len(filter_info) == 0:
            return "PASS", True
        else:
            return filter_info, True
    else:
        return "", False


def call_p_mutations(
    alt_pos, cov, p_vals, min_alt_pos, only_call_if_rare=False
):
    """Calls a p-mutation at a position, given many values of p to check.

    Parameters
    ----------
    alt_pos: int >= 0

    cov: int >= 0

    p_vals: list of float in (0, 50]

    min_alt_pos: int >= 1

    only_call_if_rare: bool

    Returns
    -------
    (filter_info, any_pass): (str, bool)
        If there exist any occurrences of True in results, then any_pass will
        be True. Otherwise, it will be False.
    """
    # This implicitly avoids the messed-up case where alt(pos) and cov are both
    # zero, since min_alt_pos must be at least 1
    if alt_pos >= min_alt_pos:
        # We call a p-mutation if alt(pos) / reads(pos) >= p / 100.
        # Equivalently: we call a p-mutation if 100*alt(pos) >= p*reads(pos).
        #
        # This way, we avoid division; we aren't entirely out of the
        # woods of potential floating-point errors, but this should be
        # more reliable, I think. (And thanks to Python's support for
        # arbitrary-precision integers, we shouldn't need to worry
        # about numbers getting ridiculously large enough to cause
        # overflow problems -- plus, like, the the max value on the
        # right hand side here would be what, p = 50 times a coverage
        # of say 1,000,000x? That's no biggie.)

        lhs = 100 * alt_pos

        filter_info = ""
        any_pass = False

        for p in p_vals:
            rhs = p * cov
            is_mut = False
            if lhs >= rhs:
                # This position counts as a p-mutation, but we may still need
                # to make the "is this a rare mutation?" check.
                if only_call_if_rare:
                    is_mut = is_position_rare_direct(alt_pos, cov)
                else:
                    is_mut = True
            if is_mut:
                any_pass = True
            else:
                if len(filter_info) > 0:
                    filter_info += ";"
                filter_info += p2filter(p)

        return convert_calling_output(filter_info, any_pass)
    else:
        # NOTE: this is silly, just update the docstring to handle this case to
        # avoid doing extra work
        return "", False


def call_r_mutations(alt_pos, r_vals):
    """Calls a r-mutation at a position, given many values of r to check.

    Parameters
    ----------
    alt_pos: int >= 0

    r_vals: list of int >= 1

    Returns
    -------
    (filter_info, any_pass): (str, bool)
        results is a list with the same length as r_vals. Each entry in results
        corresponds to an entry in r_vals: a True indicates that an r-mutation
        was called at this position for this value of r, and a False indicates
        that an r-mutation was not called at this position for this value of r.

        If there exist any occurrences of True in results, then any_pass will
        be True. Otherwise, it will be False.
    """
    filter_info = ""
    any_pass = False
    for r in r_vals:
        if alt_pos >= r:
            any_pass = True
        else:
            if len(filter_info) > 0:
                filter_info += ";"
            filter_info += r2filter(r)
    return convert_calling_output(filter_info, any_pass)
