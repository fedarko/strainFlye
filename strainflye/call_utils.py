# Utilities for strainFlye's call step.


import pysam
from . import cli_utils


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


def call_p(contigs, bam, output_vcf, min_alt_pos, p):
    return f"na\u00efve p-mutation calling at p = {p:.2f}%"


def call_r(contigs, bam, output_vcf, r):
    return f"na\u00efve r-mutation calling at r = {r:,}"
