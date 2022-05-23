# Utilities for strainFlye's call step.


import os
import time
import pysam
import pysamstats
from . import config
from .errors import SequencingDataError, ParameterError
from strainflye import __version__


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


def parse_di_list(di_list_str, param):
    """Parses and checks a list of p or r parameters.

    Useful when taking in a list of these for diversity index computation.
    Probably there's a more elegant way of handling this than having the user
    just give us a string, but... whatevs.

    Parameters
    ----------
    di_list_str: str
        Ideally, a comma-separated list of valid p or r parameters.
        In practice, people are going to throw a lot of stuff at this, so we'll
        want to do some thorough error checking.

    param: str
        Either "p" or "r", indicating how exactly to check these parameters.

    Returns
    -------
    di_list: list of int
        A sorted and cleaned-up list of parameters.

    Raises
    ------
    ParameterError
        - If any of the values are invalid.
        - If param is not "p" or "r" (come on, don't do this to me...)
    """
    # "trust nobody. not even yourself." -- guido van rossum probably idk
    if param != "p" and param != "r":
        raise ParameterError('param must be either "p" or "r".')
    split = di_list_str.split(",")
    di_list = []
    for s in split:
        # We proooobably don't need to do this because int() is ok with
        # surrounding whitespace, but I'm paranoid
        stripped_s = s.strip()
        try:
            val = int(stripped_s)
        except ValueError:
            raise ParameterError(
                f"We couldn't parse \"{stripped_s}\". Doesn't seem to be an "
                "integer?"
            )
        if param == "p":
            if val <= 0 or val > 5000:
                raise ParameterError(
                    f"{val} is not in the range (0, 5000], and is thus not a "
                    "valid value of p."
                )
        else:
            if val <= 0:
                raise ParameterError(
                    f"{val} is not >= 1, and is thus not a valid value of r."
                )

        di_list.append(val)

    if len(set(di_list)) != len(di_list):
        raise ParameterError(
            "The list of diversity index threshold values isn't unique."
        )
    return di_list


def run(
    contigs,
    bam,
    output_vcf,
    fancylog,
    verbose,
    min_p=None,
    min_r=None,
    min_alt_pos=None,
):
    """Launches the process of naive p- or r-mutation calling.

    Also computes diversity index information -- we're already looping through
    each position in each contig, so it's not a huge extra cost.

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

    min_p: int >= 1 or None
        Minimum p parameter, if we're calling p-mutations.
        Should be scaled up by 100, like in the CLI -- so min_p = 1 corresponds
        to 1 / 100 = 0.01%.

    min_r: int >= 1 or None
        Minimum r parameter, if we're calling r-mutations.

    min_alt_pos: int >= 1 or None
        During p-mutation calling, the second-most-common aligned nucleotide's
        frequency must be at least this. Not used when calling r-mutations.

    div_index_p_list: list of int, or None

    div_index_r_list: list of int, or None

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If something went wrong with the parameters. We don't check (at least
        within this function) that they fall into their respective allowed
        ranges (since Click should have already enforced that thanks to the
        range parameter settings), but we do check that only one of
        (min_p, min_r) is specified -- exactly one of these parameters must be
        given.
    """
    using_p = min_p is not None
    using_r = min_r is not None

    if using_p and using_r:
        raise ParameterError(
            "p and r can't be specified at the same time. Please choose one."
        )
    elif not using_p and not using_r:
        raise ParameterError("Either p or r needs to be specified.")

    # TODO? Add some checking to make sure that using_p implies that only
    # div_index_p_list is specified, and same for using_r and div_index_r_list
    # ... not high priority tho

    # The filter_header is important, since we will parse its ID later on to
    # determine what the minimum p or r value was. I'm not sure if there's
    # a better way to encode arbitrary file-level metadata in VCF files -- I
    # guess we could try to parse the "source" header, but I'm not a huge fan
    # of that.
    if using_p:
        param_name = "p"
        min_str = "--min-p = {min_p / 100:.2f}%"
        filter_header = (
            f'##FILTER=<ID=strainflye_minp_{min_p}, Description="min p '
            'threshold (scaled up by 100)">'
        )
    else:
        param_name = "r"
        min_str = "--min-r = {min_r:,}"
        filter_header = (
            f'##FILTER=<ID=minr_{min_r}, Description="min r threshold">'
        )

    call_str = f"{param_name}-mutation calling ({min_str})"
    fancylog(f"Running {call_str}.", prefix="")

    # See the VCF 4.3 docs for details about the Number meanings:
    # - The 1 indicates that we only include one version of this Number per
    #   position (because it doesn't really make sense to, for example,
    #   compute depth multiple times for a position).
    #
    # - The A indicates that we will include one version of this Number per
    #   alternate allele per position (since for now we just include at most
    #   one alt allele per position, there isn't a difference from "1", but
    #   this could change if we start identifying multiallelic mutations).
    #
    # Also, note that MDP and AAD are nonstandard info headers -- I made them
    # up, because I wasn't satisfied with DP and AD in the VCF docs.
    #
    # - MDP is explicitly *just* (mis)matching reads -- it isn't clear if DP
    #   in the VCF docs includes indels also, and I didn't want to risk
    #   misleading anyone.
    #
    # - AAD is like AD, but it just applies for the alternate allele. At least
    #   now, we don't need allele counts for the reference allele as well.
    info_header = (
        "##INFO=<ID=MDP,Number=1,Type=Integer,"
        'Description="(Mis)match read depth">\n'
        "##INFO=<ID=AAD,Number=A,Type=Integer,"
        'Description="Alternate allele read depth">'
    )

    with open(output_vcf, "w") as vcf_file:
        # Header info gleaned by reading over the VCF docs
        # (https://samtools.github.io/hts-specs/VCFv4.3.pdf) and copying how
        # LoFreq organizes their header.
        # ... (It took me an embarrassing amount of time to realize that VCF
        # had already moved on to 4.3.)

        vcf_file.write(
            "##fileformat=VCFv4.3\n"
            f"##fileDate={time.strftime('%Y%m%d')}\n"
            f'##source="strainFlye v{__version__}: {call_str}"\n'
            f"##reference={os.path.abspath(contigs)}\n"
            f"{info_header}\n"
            f"{filter_header}\n"
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
                is_mut = call_p_mutation(alt_freq, cov, min_p, min_alt_pos)
            else:
                is_mut = call_r_mutation(alt_freq, cov, min_r)

            if is_mut:
                info = get_pos_info_str(alt_freq, cov)
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
                # 7. FILTER = . (We only have one "filter" right now: given
                #               the minimum p or r value, is this position a
                #               p- or r-mutation for this minimum value? If so,
                #               this implies that this position is also a
                #               p- or r-mutation for all higher values of p
                #               or r. So no need to encode discrete "filter"
                #               info for higher values of p or r.)
                #
                # 8. INFO = info (Information about coverage and alt(pos):
                #                useful when plotting FDR curves, etc.)
                vcf_text += (
                    f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t{info}\n"
                )
                num_any_muts += 1

        if num_any_muts > 0:
            with open(output_vcf, "a") as vcf_file:
                vcf_file.write(vcf_text)

        if verbose:
            fancylog(
                (
                    f"Called {num_any_muts:,} {param_name}-mutation(s) "
                    f"(using {min_str}) in contig {seq}."
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


def get_pos_info_str(alt_pos, cov):
    return f"MDP={cov};AAD={alt_pos}"


def call_p_mutation(alt_pos, cov, p, min_alt_pos, only_call_if_rare=False):
    """Calls a p-mutation at a position.

    Parameters
    ----------
    alt_pos: int >= 0

    cov: int >= 0

    p: int >= 1 in (0, 50,000]
        (Still scaled up by 100, so p = 1 corresponds to 1 / 100 = 0.01%.)

    min_alt_pos: int >= 1

    only_call_if_rare: bool

    Returns
    -------
    bool
    """
    # This implicitly avoids the messed-up case where alt(pos) and cov are both
    # zero, since min_alt_pos must be at least 1
    if alt_pos >= min_alt_pos:
        # We call a p-mutation if alt(pos) / coverage(pos) >= p / 10,000.
        # (The 10,000 is because p is scaled up by 100, then by 100 again,
        # relative to (0, 50]. It's now on the order of (0, 50,000].)
        #
        # Equivalently: we call a p-mutation if
        # 10,000 * alt(pos) >= p * coverage(pos).
        #
        # This way, we avoid division and mucking around with floats.
        # The numbers might be large ints, but they're not THAT large.

        lhs = 10000 * alt_pos
        rhs = p * cov

        is_mut = False

        if lhs >= rhs:
            # This position counts as a p-mutation, but we may still need
            # to make the "is this a rare mutation?" check.
            if only_call_if_rare:
                is_mut = is_position_rare_direct(alt_pos, cov)
            else:
                is_mut = True

        return is_mut
    return False


def call_r_mutation(alt_pos, r):
    """Calls a r-mutation at a position.

    Parameters
    ----------
    alt_pos: int >= 0

    r: int >= 1

    Returns
    -------
    bool
    """
    return alt_pos >= r
