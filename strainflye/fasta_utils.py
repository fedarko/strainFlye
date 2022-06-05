# General IO utilities.


import skbio
from strainflye.errors import SequencingDataError, ParameterError


APOLOGY = "This isn't supported at the moment, sorry."


def get_name2len(fasta_fp, min_num_contigs=2):
    """Based on a FASTA file, creates a dict mapping sequence name to length.

    Parameters
    ----------
    fasta_fp: str
        Filepath to a FASTA file.

    min_num_contigs: int
        Minimum number of contigs that must be described in the FASTA file.
        This will raise an error if this condition is not met.

    Returns
    -------
    name2len: dict
        Maps sequence name to length.

    Raises
    ------
    FileNotFoundError
        If the file doesn't exist (skbio handles this).

    SequencingDataError
        If various things are wrong with the FASTA file, as determined by us:
          - Duplicate sequence names
          - Blank sequence name
          - Degenerate nucleotides in sequences
          - Gaps in sequences
          - Less than min_num_contigs contigs are described
    """
    name2len = {}
    # Fails nicely with a FileNotFoundError if this file doesn't exist
    # (... although that shouldn't happen much b/c we generally specify
    # click.Path(exists=True) for contig FASTA inputs)
    #
    # Note the constructor=skbio.DNA thing. This is optional, but we include it
    # -- because otherwise skbio.io.read() will yield generic "Sequence"
    # objects, which lack some useful functionalities (e.g. has_degenerates())
    for seq in skbio.io.read(fasta_fp, format="fasta", constructor=skbio.DNA):

        seq_name = seq.metadata["id"]

        # Fail if a given sequence with this name already exists in name2len,
        # since this means that there exist duplicate names in the FASTA file.
        # skbio's FASTA reader allows this, but I don't because accounting for
        # this sounds like a pain.
        #
        # For what it's worth: the average time complexity of checking if a
        # name is in the dict is O(1), per https://stackoverflow.com/a/1963514.
        # So this shouldn't be a big bottleneck, although I guess ideally we'd
        # do this check more efficiently (e.g. subclassing dict to check this
        # on adding an item? -- like https://stackoverflow.com/a/4999284)
        if seq_name in name2len:
            raise SequencingDataError(
                f"Duplicate sequence name in {fasta_fp}: {seq_name}"
            )

        # Also, if this sequence's name is blank (which is possible), raise an
        # error. Probably something is horribly wrong in this case.
        if len(seq_name) == 0:
            raise SequencingDataError(
                "There exists a blank (or just whitespace) sequence name in "
                f"{fasta_fp}"
            )

        # Lastly, skbio is kind enough to check for degenerate nucleotides and
        # gaps early, so we can fail loudly about these here.
        # (Supporting contigs with degenerates / gaps is possible, but it'll
        # require some more thinking and testing. I'd prefer to be too strict
        # than too lax.)
        if seq.has_degenerates():
            raise SequencingDataError(
                f"Sequence {seq_name} in {fasta_fp} has at least one "
                f"degenerate nucleotide. {APOLOGY}"
            )

        if seq.has_gaps():
            raise SequencingDataError(
                f"Sequence {seq_name} in {fasta_fp} has at least one gap. "
                f"{APOLOGY}"
            )

        # If we've made it here, this sequence seems kosher
        name2len[seq_name] = len(seq)

    if len(name2len) < min_num_contigs:
        raise SequencingDataError(
            f"Less than {min_num_contigs:,} contigs are given in {fasta_fp}."
        )

    return name2len


def verify_contigs_subset(child, parent, child_desc, parent_desc):
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

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If the child set is not a subset of the parent set.
        The resulting error message will include child_desc and parent_desc.
    """
    if not child.issubset(parent):
        raise ParameterError(
            f"All contigs in {child_desc} must also be contained in "
            f"{parent_desc}."
        )
