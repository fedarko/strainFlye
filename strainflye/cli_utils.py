# Utilities for strainFlye's CLI.

import time


def fancystart(cmd_name, inputs, outputs, verbose, prefix="--------\n"):
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
        something like (("VCF file", vcf_filepath)).

    verbose: bool
        Whether or not to log any messages at all. If this is False, then
        nothing will be logged from here, and calling the returned logging
        function (fancylog) will not output anything.

    prefix: str
        Prefix to put before every logging message.

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
        if verbose:
            t1 = time.time()
            print(f"{prefix}{cmd_name} @ {t1 - t0:.2f} sec: {msg}", flush=True)

    # Report to the user about the inputs and outputs.
    # ... This may be over-engineered.
    # Note that this should behave graciously if there are no inputs or
    # outputs, since these will then be empty collections (and nothing will be
    # printed). However, that probably shouldn't happen in practice.

    starting_info = "Starting..."
    for paramtype in (("Input", inputs), ("Output", outputs)):
        for ptuple in paramtype[1]:
            starting_info += f"\n{paramtype[0]} {ptuple[0]}: {ptuple[1]}"

    fancylog(starting_info)

    return fancylog
