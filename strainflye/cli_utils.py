# Utilities for strainFlye's CLI.

import time


def fancystart(cmd_name, inputs, outputs, verbose):
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

    Returns
    -------
    fancylog: function
        A function that takes one parameter, msg (of type str), and logs this
        message. Each logging message is prefixed with the command name and the
        time (in seconds) since this command started (or, more accurately,
        since fancystart() was called).
    """
    t0 = time.time()

    def fancylog(msg):
        """Logs a message."""
        if verbose:
            t1 = time.time()
            print(f"--------\n{cmd_name} @ {t1 - t0:.2f} sec: {msg}")

    fancylog("Starting running...")
    # Report to the user about the inputs and outputs
    # ... This may be over-engineered
    for paramtype in (("Input", inputs), "Output", outputs):
        for ptuple in paramtype[1]:
            fancylog(f"{paramtype[0]} {ptuple[0]}: {ptuple[1]}")

    return fancylog
