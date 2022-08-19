# General BCF-related utilities.


import os
import re
import subprocess
import pysam
from .errors import ParameterError


def convert_vcf_to_bcf(vcf_fp, bcf_fp, fancylog):
    """Converts a VCF file to a BCF file using "bcftools view".

    Notably, this will delete the VCF file after performing this conversion.

    Parameters
    ----------
    vcf_fp: str
        Path to a VCF file.

    bcf_fp: str
        Path to which a BCF file will be written.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    fancylog("Converting the VCF file to a compressed BCF file...")
    subprocess.run(["bcftools", "view", "-O", "b", vcf_fp, "-o", bcf_fp])
    os.remove(vcf_fp)
    fancylog("Done.", prefix="")


def index_bcf(bcf_fp, fancylog):
    """Indexes a BCF file using "bcftools index".

    This creates a *.bcf.csi file in the same folder as the BCF file.

    Parameters
    ----------
    bcf_fp: str
        Path to a BCF file to be indexed.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    fancylog("Indexing the BCF file...")
    subprocess.run(["bcftools", "index", bcf_fp])
    fancylog("Done indexing the BCF file.", prefix="")


def compress_vcf(vcf_fp, bcf_fp, fancylog):
    """Converts a VCF to a BCF file, and indexes this BCF file.

    Parameters
    ----------
    vcf_fp: str
        Path to a VCF file.

    bcf_fp: str
        Path to which a BCF file will be written.
        We'll also create a .bcf.csi file in the same folder as the BCF file.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    convert_vcf_to_bcf(vcf_fp, bcf_fp, fancylog)
    index_bcf(bcf_fp, fancylog)


def verify_bcf_has_contigs_with_lengths(bcf_obj, bcf_fp):
    if len(bcf_obj.header.contigs) < 1:
        raise ParameterError(
            f"BCF file {bcf_fp} doesn't describe any contigs in its header."
        )

    # The VCF 4.3 specification -- as of writing -- only says that contig tags
    # "typically" include length information. So, let's guarantee the presence
    # of length info here to make life easier for us downstream.
    for contig in bcf_obj.header.contigs:
        contig_entry = bcf_obj.header.contigs[contig]
        # From manual testing, it looks like pysam does actually set the length
        # attribute of contigs without a given length -- but it sets it to
        # None. Just to future-proof this a bit, we check both the case where a
        # contig has no length attribute and the case where a contig has a None
        # length attribute. Both are invalid.
        if not hasattr(contig_entry, "length") or contig_entry.length is None:
            raise ParameterError(
                f"BCF file {bcf_fp} has no length given for contig {contig}."
            )


def verify_no_multiallelic_mutations_in_bcf(bcf_obj, bcf_fp):
    for contig in bcf_obj.header.contigs:
        seen_positions = set()
        for mut in bcf_obj.fetch(contig):
            if mut.pos in seen_positions:
                raise ParameterError(
                    f"BCF file {bcf_fp} has multiple mutations at "
                    f"(1-indexed) position {mut.pos:,} on contig {contig}. "
                    "strainFlye does not currently support BCF "
                    "files containing multi-allelic mutations, sorry."
                )
            seen_positions.add(mut.pos)


def parse_arbitrary_bcf(bcf):
    """Opens a BCF file that may or may not have been produced by strainFlye.

    Unlike parse_sf_bcf(), this function doesn't require that the input
    mutations come from strainFlye call (or strainFlye fdr fix). This will
    still perform some sanity checking on the BCF file (for example, making
    sure that it does not contain any multi-allelic mutations or indels),
    though. The main purpose of this is to enable us to use arbitrary BCF files
    with strainFlye's "downstream" commands (e.g. hotspot/coldspot detection,
    phasing).

    Parameters
    ----------
    bcf: str
        Path to a BCF file.

    Returns
    -------
    f: pysam.VariantFile
        Object describing the BCF file located at the specified filepath.
    """
    f = pysam.VariantFile(bcf)
    verify_bcf_has_contigs_with_lengths(f, bcf)
    verify_no_multiallelic_mutations_in_bcf(f, bcf)
    return f


def parse_sf_bcf(bcf):
    """Opens a BCF file that we expect to have been produced by strainFlye.

    Does sanity checking and p vs. r sniffing on the file.
    Thankfully, pysam has the ability to read BCF files, so the main thing
    we do here is checking that the meta-information of the BCF file seems
    kosher (i.e. was produced by strainFlye).

    Possible TODO: add another return value, indicating whether or not a header
    indicating FDR fixing was done is present in this BCF. The caller could
    then optionally log a warning if this is or isn't the case?

    Parameters
    ----------
    bcf: str
        Path to a BCF file produced by one of "strainFlye call"'s subcommands.

    Returns
    -------
    (f, thresh_type, thresh_min): (pysam.VariantFile, str, int)
        f: Object describing the BCF file located at the specified filepath.
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

        Also, we raise this sort of error if this file doesn't describe any
        contigs (because in that case, you should consult a priest rather than
        strainFlye).
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

    # Similarly, we verify some extra things about the BCF that *should* hold
    # for strainFlye output but you never know
    verify_bcf_has_contigs_with_lengths(f, bcf)
    verify_no_multiallelic_mutations_in_bcf(f, bcf)

    return f, thresh_type, thresh_min


def loudly_parse_arbitrary_bcf_and_contigs(bcf, fancylog):
    """Calls parse_arbitrary_bcf() on a BCF file while logging about it.

    Encapsulated this to its own function because this was the same for both
    "spot" commands.

    Parameters
    ----------
    bcf: str
        Filepath to a BCF file describing single-nucleotide mutations.

    fancylog: function
        Logging function.

    Returns
    -------
    (bcf_obj, bcf_contigs): (pysam.VariantFile, list)
        BCF object and collection of contigs described in its header.
        (Presumably, the order of contigs matches the order of contigs in the
        BCF header, but we probably can't rely on that in 100% of cases.)

    Raises
    ------
    Any of the errors that parse_arbitrary_bcf() would raise; see that
    function's documentation for details.
    """
    fancylog("Loading and checking the BCF file...")
    bcf_obj = parse_arbitrary_bcf(bcf)
    fancylog("Looks good so far.", prefix="")
    # We don't really NEED to convert this to a list but we might as well just
    # for peace of mind
    bcf_contigs = list(bcf_obj.header.contigs)
    return bcf_obj, bcf_contigs


def verify_contig_in_bcf(bcf_obj, contig):
    """Raises an error if a contig is not described in a BCF object's header.

    Parameters
    ----------
    bcf_obj: pysam.VariantFile
        Object describing a BCF file.

    contig: str
        Name of a contig.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If contig is not described in the BCF object's header.
    """
    if contig not in bcf_obj.header.contigs:
        raise ParameterError(
            f"Contig {contig} is not described in the BCF object's header."
        )


def get_mutated_positions_in_contig(bcf_obj, contig, zero_indexed=True):
    """Identifies positions containing mutations in a contig.

    We use the same definition of "mutation" (a single-nucleotide,
    single-allelic variant) defined elsewhere in our paper and code, although
    we don't explicitly check this here -- we just assume that the input BCF
    file only contains these "simple" mutations.

    Parameters
    ----------
    bcf_obj: pysam.VariantFile
        Object describing a BCF file.

    contig: str
        Name of a contig for which we will fetch mutations in bcf_obj.

    zero_indexed: bool
        If True, return zero-indexed positions; otherwise, return one-indexed
        positions.

    Returns
    -------
    mutated_positions: set
        Set of positions (integers), corresponding to the positions on the
        contig at which mutations are defined in the BCF file. If the contig
        has no mutations, then this set will be empty.

    Raises
    ------
    ParameterError
        If contig is not described in the header of bcf_obj.
    """
    verify_contig_in_bcf(bcf_obj, contig)
    # I guess we could return a list instead, but I am 100% not going to make
    # the assumption that pysam.VariantFile.fetch() always returns positions in
    # sorted order, so it's the caller's responsibility to call sorted() on
    # this if they want a sorted list of these positions. Better slow and
    # correct than fast and wrong.
    mutated_positions = set()

    # pysam gives us 1-based positions, so if we want zero-indexed positions we
    # gotta subtract 1 from everything. This is probably an inefficient way to
    # do this but whatever if this is the bottleneck then god help us all
    def index_position(pos):
        if zero_indexed:
            return pos - 1
        else:
            return pos

    for mut in bcf_obj.fetch(contig):
        mutated_positions.add(index_position(mut.pos))

    return mutated_positions


def get_mutated_position_details_in_contig(bcf_obj, contig):
    """Returns a summary of mutated-position details in a contig.

    Parameters
    ----------
    bcf_obj: pysam.VariantFile
        Object describing a BCF file. We make the assumption that this BCF file
        does not contain multi-allelic or non-single-nucleotide mutations (this
        should already have been screened for when parsing the BCF).

    contig: str
        Name of a contig for which we will fetch mutations in bcf_obj.

    Returns
    -------
    mp2ra: dict
        Maps (zero-indexed) mutated positions in the contig to a tuple of
        (ref nt, alt nt), as listed in the BCF file. If the contig has no
        mutations, then this dict will be empty.

    Raises
    ------
    ParameterError
        If contig is not described in the header of bcf_obj.
    """
    verify_contig_in_bcf(bcf_obj, contig)
    mp2ra = {}
    for mut in bcf_obj.fetch(contig):
        # mut.pos is 1-indexed, so gotta subtract 1 to make it 0-indexed
        mp2ra[mut.pos - 1] = (mut.ref, mut.alts[0])
    return mp2ra
