# Utilities for strainFlye fdr.


import re
import pysam
from .errors import ParameterError


def parse_bcf(bcf):
    """Opens a BCF file, and does sanity checking and p vs. r sniffing on it.

    Thankfully, pysam has the ability to read BCF files, so the main thing
    we do here is checking that the meta-information of the BCF file seems
    kosher (i.e. was produced by strainFlye).

    Possible TODO: add another return value, indicating whether or not a header
    indicating FDR fixing was done is present in this BCF. The caller could
    then optionally log a warning if so?

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
