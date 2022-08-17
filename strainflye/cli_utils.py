# Utilities for strainFlye's CLI.

import time


def fancystart(
    cmd_name, inputs, outputs, quiet=False, prefix="--------\n", extra_info=()
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

    quiet: bool
        If this is True, then nothing will be logged from here, and calling
        the returned logging function (fancylog) will not output anything.
        This isn't used anywhere yet, but if folks ask for a "quiet" option we
        can add it using this. (I know, I know, YAGNI...)

    prefix: str
        Prefix to put before every logging message.

    extra_info: tuple of str
        A collection of extra parameter information; each str will be printed
        on its own line (well, maybe multiple lines if these strings contain
        newlines, but probably don't do that) after inputs but before outputs.
        Useful for adding info for flags, etc. An example:
        ("Check for circular coldspot gaps?: No", "Verbose?: Yes").

    Returns
    -------
    fancylog: function
        A simple logging function. The only required parameter to this function
        is a string message. Each logging message is prefixed with the command
        name and the time (in seconds) since this command started (or, more
        accurately, since fancystart() was called).
    """
    t0 = time.time()

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
        if not quiet:
            t1 = time.time()
            print(
                f"{prefix}{cmd_name} @ {t1 - t0:,.2f} sec: {msg}", flush=True
            )

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


def log_prog(
    contig, contig_len, contig_1idx, num_contigs, fancylog, prefix="On "
):
    """Log about progress when iterating through a set of contigs.

    This is just abstracted code that turns up in a few places.
    """
    pct = 100 * (contig_1idx / num_contigs)
    fancylog(
        (
            f"{prefix}contig {contig} ({contig_len:,} bp) "
            f"({contig_1idx:,} / {num_contigs:,} = {pct:.2f}% done)."
        ),
        prefix="",
    )
