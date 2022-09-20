# Utilities for strainFlye's CLI.

import time
from strainflye import __version__


def fancystart(
    cmd_name, inputs, outputs, prefix="--------\n", extra_info=(), version=True
):
    """Starts logging things for a given strainFlye command.

    Parameters
    ----------
    cmd_name: str
        Name of the command we're running.

    inputs: tuple of 2-tuples of (str, str)
        Each 2-tuple corresponds to an input to the command. For each
        2-tuples, the first element should be a description of this
        input (e.g. "contig file"), and the second element should be a string
        representing the input specified (e.g. the filepath to the
        corresponding FASTA file of contigs).

    outputs: tuple of 2-tuples of (str, str)
        A collection of command output information, formatted analogously to
        how "inputs" is formatted. For example, a variant caller might use
        something like (("BCF file", bcf_filepath)).

    prefix: str
        Prefix to put before every logging message, by default.

    extra_info: tuple of str
        A collection of extra parameter information; each str will be printed
        on its own line (well, maybe multiple lines if these strings contain
        newlines, but probably don't do that) after inputs but before outputs.
        Useful for adding info for flags, etc. An example:
        ("Check for circular coldspot gaps?: No", "Verbose?: Yes").

    version: bool
        If True, show a message before the starting text about the detected
        strainFlye version; if False, don't.

    Returns
    -------
    fancylog: function
        A simple logging function. The only required parameter to this function
        is a string message. Each logging message is prefixed with the command
        name and the time (in seconds) since this command started (or, more
        accurately, since fancystart() was called).
    """
    t0 = time.time()

    # definitely overkill
    fancyprint = lambda text: print(text, flush=True)

    def fancylog(msg, prefix=prefix):
        """Logs a message.

        By default, the prefix (before the command name and time) matches what
        was passed as the prefix parameter to fancystart(), but this can be
        overridden here if desired.

        Parameters
        ----------
        msg: str
            Message to log.

        prefix: str
            Prefix to put before a logging message (and before the command name
            and time).

        Returns
        -------
        None
        """
        t1 = time.time()
        fancyprint(f"{prefix}{cmd_name} @ {t1 - t0:,.2f}s: {msg}")

    if version:
        fancyprint(f'Using strainFlye version "{__version__}".')

    # Report to the user about the inputs, parameters / settings, and outputs.
    # ... This may be over-engineered.
    # Note that this should behave graciously if there are no inputs or
    # outputs, since these will then be empty collections (and nothing will be
    # printed). However, that probably shouldn't happen in practice.

    starting_info = "Starting..."

    for ituple in inputs:
        starting_info += f"\nInput {ituple[0]}: {ituple[1]}"
    for eline in extra_info:
        starting_info += f"\n{eline}"
    for otuple in outputs:
        starting_info += f"\nOutput {otuple[0]}: {otuple[1]}"

    fancylog(starting_info)

    return fancylog


def b2y(b):
    """I knew that Bachelor's degree was good for something."""
    if b:
        return "Yes"
    else:
        return "No"


def proglog(
    contig, contig_1idx, num_contigs, fancylog, contig_len=None, prefix="On "
):
    """Log about progress when iterating through a set of contigs.

    This is just abstracted code that turns up in a few places.
    """
    pct = 100 * (contig_1idx / num_contigs)
    lenstr = " "
    if contig_len is not None:
        lenstr = f" ({contig_len:,} bp) "

    # i'm fancy
    contig_noun = "contigs"
    if num_contigs == 1:
        contig_noun = "contig"

    fancylog(
        (
            f"{prefix}contig {contig}{lenstr}"
            f"({contig_1idx:,} / {num_contigs:,} {contig_noun} = {pct:.2f}%)."
        ),
        prefix="",
    )


def get_verboselog(fancylog, verbose):
    """Returns a new function which calls fancylog only if verbose is True."""

    def verboselog(*args, **kwargs):
        if verbose:
            fancylog(*args, **kwargs)

    return verboselog
