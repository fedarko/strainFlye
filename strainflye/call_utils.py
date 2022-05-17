# Utilities for strainFlye's call step.


import os
import time
import pysam
import pysamstats
from . import cli_utils, config
from .errors import SequencingDataError
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
    (cov, alt nt freq, alt nt): tuple of (int, int, str)
        Describes the second-most-common nucleotide at a position.

        The first entry in this tuple is the (mis)match coverage at this
        position. This is an integer defined as the sum of A, C, G, T
        nucleotides at this position (note that this excludes degenerate
        nucleotides like N -- we could change this in the future if that'd be
        useful, I suppose). Note that this coverage could be zero, if no reads
        are aligned to this specific position.

        The second entry is the raw frequency of this nucleotide
        at this position: this will be an integer greater than or equal to 0.
        This is also referred to in the paper, etc. as alt(pos).

        The third entry is just the alternate nucleotide (one of A, C, G, T),
        represented as a string. This is returned for reference -- as of
        writing this isn't actually needed for Pleuk itself, but I have other
        code outside of Pleuk that benefits from this!

    References
    ----------
    I mean, I wrote this code, but it's copied from
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
    alt_nt_freq = rec[alt_nt]

    return (cov, alt_nt_freq, alt_nt)


def run(contigs, bam, output_vcf, min_alt_pos, p=None, r=None):
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

    min_alt_pos: int >= 0
        During p-mutation calling, the second-most-common aligned nucleotide's
        frequency must be at least this.

    p: int in [0, 50]
        p-mutation parameter.

    r: int >= 0
        r-mutation parameter.

    Returns
    -------
    call_str: str
        Description of the type of mutation calling that just happened.

    Raises
    ------
    cli_utils.ParameterError
        If something went wrong with the parameters. We don't check that they
        fall into their respective allowed ranges (since Click should have
        already enforced that thanks to the range parameter settings), but we
        do check that only p or only r is specified -- we can't have both
        specified at once, or neither specified.
    """
    using_p = p is not None
    using_r = r is not None

    if using_p:
        if using_r:
            raise cli_utils.ParameterError(
                "p and r can't be specified at the same time. Please choose "
                "one."
            )
        else:
            # p-mutation calling
            return call_p(contigs, bam, output_vcf, min_alt_pos, p)

    elif using_r:
        # r-mutation calling
        return call_r(contigs, bam, output_vcf, r)

    else:
        raise cli_utils.ParameterError("Either p or r needs to be specified.")

    # Header info gleaned by reading over the VCF 4.2 docs
    # (https://samtools.github.io/hts-specs/VCFv4.2.pdf) and copying how LoFreq
    # organizes their header
    vcf = (
        "##fileformat=VCFv4.2\n"
        f"##fileDate={time.strftime('%Y%m%d')}\n"
        f'##source="strainFlye v{__version__}"\n'
        f"##reference={os.path.abspath(contigs)}\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    bf = pysam.AlignmentFile(bam, "rb")
    for si, seq in enumerate(bf.references, 1):
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

            cov, alt_nt_freq, alt_nt = get_alt_pos_info(rec)

            if using_p:
                if call_p_mutation(alt_nt_freq, cov, p, min_alt_pos):
                    pass  # TODO add to VCF
            else:
                if call_r_mutation(alt_nt_freq, r):
                    pass  # TODO add to VCF

    if using_p:
        return f"na\u00efve p-mutation calling at p = {p:.2f}%"
    else:
        return f"na\u00efve r-mutation calling at r = {r:,}"


def is_position_rare_direct(alt_pos, cov):
    lhs = 100 * alt_pos
    rhs_upper = config.HIGH_FREQUENCY_MIN_PCT * cov
    return lhs < rhs_upper


def call_p_mutation(alt_pos, cov, p, min_alt_pos, only_call_if_rare=False):
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

        lhs = 100 * alt_nt_freq
        rhs = p * cov

        if lhs >= rhs:
            # This position counts as a p-mutation, but we may still need to
            # make the "is this a rare mutation?" check.
            if only_call_if_rare:
                return is_position_rare_direct(alt_pos, cov)
            else:
                return True
    else:
        return False


def call_r_mutation(alt_pos, r):
    return alt_pos >= r
