# Miscellaneous utilities.

import os
import json
import pickle
import pysam
import pandas as pd
from strainflye import config, fasta_utils, bcf_utils
from strainflye.errors import ParameterError


def make_output_dir(output_dir):
    """Creates an output directory, if it doesn't already exist."""
    os.makedirs(output_dir, exist_ok=True)


def write_obj_to_pickle(obj, output_dir, contig_name, obj_name):
    """Writes out an object to a pickle file.

    This file will be named [contig_name]_[obj_name].pickle, and will be
    located in the directory [output_dir].

    Parameters
    ----------
    obj: object
        Something to write out to a pickle file.

    output_dir: str
        Directory to which we'll write this pickle file. Should already exist.

    contig_name: str
        Used as the first part of the filename.

    obj_name: str
        Used as the second part of the filename.

    Returns
    -------
    None
    """
    fp = os.path.join(output_dir, f"{contig_name}_{obj_name}.pickle")
    with open(fp, "wb") as dumpster:
        pickle.dump(obj, dumpster)


def write_obj_to_json(obj, output_dir, contig_name, obj_name):
    """Writes out an object to a JSON file.

    This file will be named [contig_name]_[obj_name].json, and will be
    located in the directory [output_dir].

    (So, same schtick as write_obj_to_pickle().)

    Note that, in JSON, keys are strings (see
    https://docs.python.org/3/library/json.html) -- so if you have e.g. int
    keys then these will be represented as strings in the JSON.

    Parameters
    ----------
    obj: object
        Something to write out to a JSON file.

    output_dir: str
        Directory to which we'll write this JSON file. Should already exist.

    contig_name: str
        Used as the first part of the filename.

    obj_name: str
        Used as the second part of the filename.

    Returns
    -------
    None
    """
    fp = os.path.join(output_dir, f"{contig_name}_{obj_name}.json")
    with open(fp, "w") as dumpster:
        json.dump(obj, dumpster)


def load_from_pickle(fp):
    """Loads an object from a pickle file."""
    # nothing special here, just abstracted code i wanted to reuse quickly
    with open(fp, "rb") as loadster:
        return pickle.load(loadster)


def check_executable(fp):
    """Checks that a file is executable. If not, raises a ParameterError.

    Parameters
    ----------
    fp: str
        Filepath to a file whose "executability" (is that a word?) we will
        check.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If fp does not point to an executable file.

    Notes
    -----
    If fp points to a directory, this function will still work -- and I think
    it should always succeed because directories are technically "executable"
    (see https://superuser.com/a/169418). That said, this function is designed
    to test *files* for being executable or not.

    References
    ----------
    Originally, I used click's "executable" parameter (added in click 8.1.0,
    https://click.palletsprojects.com/en/8.1.x/changes/#version-8-1-0) to
    perform this check -- but click 8.1.0 requires Python >= 3.7, and thus
    breaks Python 3.6 support. (So if someone would try to install strainFlye
    into a Python 3.6 environment, they'd install an older version of click,
    and then get an error about our use of the "executable" parameter.)

    To avoid this problem, I added this function, which re-uses the same check
    that click performs when the "executable" parameter is given: see
    https://github.com/pallets/click/blob/a8910b382d37cce14adeb44a73aca1d4e87c2413/src/click/types.py#L903-L910
    """
    if not os.access(fp, os.X_OK):
        raise ParameterError(f"{fp} is not executable.")


def verify_contig_subset(child, parent, child_desc, parent_desc, exact=False):
    """Verifies that one set of contig names is a subset of another set.

    Parameters
    ----------
    child: set
        Set of contig names.

    parent: set
        Set of contig names.

    child_desc: str
        Human-readable description of the child set of contig names.

    parent_desc: str
        Human-readable description of the parent set of contig names.

    exact: bool
        If True, this makes the check stricter -- it now ensures that the two
        sets are identical (so parent must also be a subset of child).

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        - If the child set is not a subset of the parent set.
        - If exact is True, and the parent set and child set are not identical.
        The resulting error message will include child_desc and parent_desc.
    """
    if not child.issubset(parent):
        raise ParameterError(
            f"All contigs in {child_desc} must also be contained in "
            f"{parent_desc}."
        )
    # Equivalently, we could check "exact and not parent.issubset(child)", but
    # I think this way of writing it is clearer
    if exact and parent != child:
        raise ParameterError(
            f"All contigs in {parent_desc} must also be contained in "
            f"{child_desc}."
        )


def verify_contig_lengths(fasta_name2len, bam_obj=None, bcf_obj=None):
    """Verifies that contig lengths from a FASTA file match other sources.

    The only required parameter is fasta_name2len. You can provide either both
    bam_obj and bcf_obj, or just one of them.

    Please note that we only check that the lengths match for contigs in the
    FASTA file (so if the BCF or BAM file contain contigs that aren't in the
    FASTA file, then these "extra" contigs will be implicitly ignored).

    Also, note that we assume that you've already verified that the *names* of
    the contigs match up between files -- if any of the FASTA contigs aren't
    in the BAM or BCF at all, then this will result in uncaught ugly errors.
    (Yeah, it would be most efficient to combine the names check in with
    this, but let's keep this simple for now...)

    Parameters
    ----------
    fasta_name2len: dict
        Maps contig names to lengths. Can be generated by
        fasta_utils.get_name2len().

    bcf_obj: pysam.VariantFile or None
        Object describing a BCF file, produced by pysam. In practice, you
        can get this from bcf_utils.parse_*_bcf(). We make the assumption that
        all contigs defined in the header of this file have the "length"
        attribute, and that this length attribute is not None; this should
        already have been enforced by bcf_utils.parse_*_bcf().

    bam_obj: pysam.AlignmentFile or None
        Object describing a BAM file, produced by pysam.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If, for every contig in fasta_name2len, its length does not exactly
        match its corresponding length in bam_obj or in bcf_obj.

        If both bam_obj and bcf_obj are None.
    """
    chk_bam = bam_obj is not None
    chk_bcf = bcf_obj is not None

    if not chk_bam and not chk_bcf:
        raise ParameterError("Neither bam_obj nor bcf_obj is provided.")

    for contig in fasta_name2len:
        faslen = fasta_name2len[contig]

        # We could make the error message here fancier if needed (e.g. once we
        # see a mismatch in the BAM, then check the BCF also), but let's keep
        # this simple and easy for now
        if chk_bam:
            bamlen = bam_obj.get_reference_length(contig)
            if faslen != bamlen:
                raise ParameterError(
                    f"Contig {contig} has length {faslen:,} in the FASTA "
                    f"file, but length {bamlen:,} in the BAM file."
                )
        if chk_bcf:
            bcflen = bcf_obj.header.contigs[contig].length
            if faslen != bcflen:
                raise ParameterError(
                    f"Contig {contig} has length {faslen:,} in the FASTA "
                    f"file, but length {bcflen:,} in the BCF file."
                )


def load_and_sanity_check_diversity_indices(
    diversity_indices, min_num_contigs=1, min_num_di_columns=1
):
    """Loads and validates a file containing diversity index information.

    Parameters
    ----------
    diversity_indices: str
        Filepath to a TSV file containing diversity index information.
        Generated by "strainFlye call."

    min_num_contigs: int
        If the diversity index file contains less than this many contigs,
        we'll throw an error. (I guess you could set this to zero to
        effectively stop this check, but this should probably always be one,
        right?)

    min_num_di_columns: int
        If the diversity index file contains less than this many columns
        whose names begin with config.DI_PREF, we'll throw an error.

    Returns
    -------
    di: pd.DataFrame
        Represents the input diversity index file. The indices are contigs; the
        columns are, well, columns in the diversity index TSV file.

    Raises
    ------
    ParameterError
        If the diversity index file does not have the expected numbers of
        contig(s) and diversity index column(s), or if it does not have Length
        or AverageCoverage columns.

    Other errors
        Can be raised by pd.read_csv() if the file seems malformed.
        We don't attempt to catch these errors.
    """
    di = pd.read_csv(diversity_indices, sep="\t", index_col=0)

    if len(di.index) < min_num_contigs:
        raise ParameterError(
            f"Diversity indices file describes < {min_num_contigs:,} "
            "contig(s)."
        )

    if "Length" not in di.columns or "AverageCoverage" not in di.columns:
        raise ParameterError(
            'Diversity indices file must include the "Length" and '
            '"AverageCoverage" columns.'
        )

    num_di_cols = 0
    for col in di.columns:
        if col.startswith(config.DI_PREF):
            num_di_cols += 1
    if num_di_cols < min_num_di_columns:
        raise ParameterError(
            f"Diversity indices file describes < {min_num_di_columns:,} "
            "column(s) of diversity indices."
        )
    return di


def load_fasta_and_bam(contigs, bam, fancylog, min_num_contigs=1):
    """Loads and checks a FASTA and BAM file.

    Ensures that all contigs in the FASTA file are also in the BAM file. We
    allow the BAM file to contain "extra" contigs.

    If you need to load a BCF file at the same time, you should use
    load_triplet() instead. This is just designed for situations where you
    don't have a BCF you care about.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a (sorted and indexed) BAM file mapping reads to contigs.

    fancylog: function
        Logging function.

    min_num_contigs: int
        Will be passed to fasta_utils.get_name2len().

    Returns
    -------
    (contig_name2len, bam_obj, num_fasta_contigs): (dict, pysam.AlignmentFile,
                                                    int)

        contig_name2len: dict mapping contig name to length.

        bam_obj: Object describing the BAM file.

        num_fasta_contigs: Number of contigs in contig_name2len.

    Raises
    ------
    This function doesn't raise any errors itself, but it calls various
    functions which can raise errors if the input files are invalid in certain
    ways. See fasta_utils.get_name2len(), verify_contig_subset(), and
    verify_contig_lengths() for details.
    """
    fancylog("Loading and checking FASTA and BAM files...")

    contig_name2len = fasta_utils.get_name2len(
        contigs, min_num_contigs=min_num_contigs
    )

    bam_obj = pysam.AlignmentFile(bam, "rb")

    # Verify that all contigs in the FASTA are also references in the BAM
    # (this will throw an error if not)
    verify_contig_subset(
        set(contig_name2len),
        set(bam_obj.references),
        "the FASTA file",
        "the BAM file",
    )

    # ... and that recorded contig lengths match
    verify_contig_lengths(contig_name2len, bam_obj=bam_obj)

    num_fasta_contigs = len(contig_name2len)

    # unnecessarily fancy, or me not wanting to break existing tests? you
    # decide

    fancylog(
        f"The FASTA file describes {num_fasta_contigs:,} contig(s).", prefix=""
    )
    fancylog(
        (
            "All of these are included in the BAM file (which has "
            f"{bam_obj.nreferences:,} reference(s)), with the same lengths."
        ),
        prefix="",
    )

    return contig_name2len, bam_obj, num_fasta_contigs


def load_triplet(
    contigs,
    bam,
    bcf,
    fancylog,
    bcf_exact=False,
    min_num_contigs=1,
    get_sf_bcf_details=False,
):
    """Loads and checks three files: FASTA, BAM, and BCF.

    Mainly, this ensures that the contigs in the FASTA file are all present in
    the BAM and BCF files. If bcf_exact is True, then this will also make sure
    that there aren't any "extra" contigs in the BCF file.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a (sorted and indexed) BAM file mapping reads to contigs.

    bcf: str
        Filepath to an (indexed) BCF file describing single-nucleotide
        mutations in contigs.

    fancylog: function
        Logging function.

    bcf_exact: bool
        If True, ensure that the sets of contigs in the FASTA and BCF files are
        identical (still allowing there to be "extra" contigs in the BAM file).
        If False, allow there to be "extra" contigs in the BCF file as well.

    min_num_contigs: int
        Will be passed to fasta_utils.get_name2len().

    get_sf_bcf_details: bool
        If True, load the BCF (and the contained threshold type and minimum)
        using bcf_utils.parse_sf_bcf(); if False, just load the BCF using
        bcf_utils.parse_arbitrary_bcf().

        This impacts both how the BCF is validated during loading, as well
        as whether or not the thresh_type and thresh_min return values from
        this function will be None.

    Returns
    -------
    (contig_name2len, bam_obj, bcf_obj, tt, tm): (dict, pysam.AlignmentFile,
                                                  pysam.VariantFile,
                                                  str or None, int or None)

        contig_name2len: dict mapping contig name to length.

        bam_obj: Object describing the BAM file.

        bcf_obj: Object describing the BCF file.

        tt: If get_sf_bcf_details is True, then this will be the threshold type
            of the mutations described in the BCF file (either "p" or "r"); if
            get_sf_bcf_details is False, then this will be None.

        tm: If get_sf_bcf_details is True, then this will be the min threshold
            of the mutations described in the BCF file (see
            bcf_utils.parse_sf_bcf() for details); if get_sf_bcf_Details is
            False, then this will be None.

    Raises
    ------
    This function doesn't raise any errors itself, but it calls various
    functions which can raise errors if the input files are invalid in certain
    ways. See fasta_utils.get_name2len(), verify_contig_subset(),
    verify_contig_lengths(), bcf_utils.parse_sf_bcf(), and
    bcf_utils.parse_arbitrary_bcf() for more details on the sorts of errors
    that can get raised here.
    """
    fancylog("Loading and checking FASTA, BAM, and BCF files...")

    contig_name2len = fasta_utils.get_name2len(
        contigs, min_num_contigs=min_num_contigs
    )
    fasta_contigs = set(contig_name2len)
    fancylog(
        f"The FASTA file describes {len(fasta_contigs):,} contig(s).",
        prefix="",
    )

    bam_obj = pysam.AlignmentFile(bam, "rb")
    verify_contig_subset(
        fasta_contigs,
        set(bam_obj.references),
        "the FASTA file",
        "the BAM file",
        exact=False,
    )
    fancylog(
        (
            "All FASTA contig(s) are included in "
            f"the BAM file (this BAM file has {bam_obj.nreferences:,} "
            "reference(s))."
        ),
        prefix="",
    )

    thresh_type = None
    thresh_min = None
    if get_sf_bcf_details:
        bcf_obj, thresh_type, thresh_min = bcf_utils.parse_sf_bcf(bcf)
    else:
        bcf_obj = bcf_utils.parse_arbitrary_bcf(bcf)

    bcf_contigs = set(bcf_obj.header.contigs)
    verify_contig_subset(
        fasta_contigs,
        bcf_contigs,
        "the FASTA file",
        "the BCF file",
        exact=bcf_exact,
    )
    if not bcf_exact:
        fancylog(
            (
                "All FASTA contig(s) are included in "
                "the BCF file (the header of this BCF file describes "
                f"{len(bcf_contigs):,} contig(s))."
            ),
            prefix="",
        )
    else:
        fancylog(
            "The FASTA file's contig(s) and BCF file's contig(s) match.",
            prefix="",
        )

    if get_sf_bcf_details:
        fancylog(
            f"Also, the input BCF file describes {thresh_type}-mutations "
            f"(minimum {thresh_type} = {thresh_min:,}).",
            prefix="",
        )

    verify_contig_lengths(contig_name2len, bam_obj=bam_obj, bcf_obj=bcf_obj)
    fancylog(
        (
            "The lengths of all contig(s) in the FASTA file match the "
            "corresponding lengths in the BAM and BCF files."
        ),
        prefix="",
    )
    fancylog("So far, these files seem good.", prefix="")
    return contig_name2len, bam_obj, bcf_obj, thresh_type, thresh_min
