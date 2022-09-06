# General IO utilities.


import skbio
from strainflye.errors import SequencingDataError


APOLOGY = "This isn't supported at the moment, sorry."


def get_name2len(fasta_fp, min_num_contigs=1):
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
        If the file doesn't exist (scikit-bio handles this).

    SequencingDataError
        If various things are wrong with the FASTA file, as determined by us:
          - Duplicate sequence names
          - Blank sequence name
          - Degenerate nucleotides in sequences
          - Gaps in sequences
          - Less than min_num_contigs contigs are described

    FASTAFormatError
        Can be raised by scikit-bio while trying to parse the FASTA file.
        Notably, if any of the sequences are empty (i.e. the sequence line is
        empty and/or just whitespace), then scikit-bio will throw this -- this
        case is explicitly tested, too (see
        test_get_name2len_empty_contig_sequence()).
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
            f"Less than {min_num_contigs:,} contig(s) are given in {fasta_fp}."
        )

    return name2len


def get_single_seq(fasta_fp, contig_name):
    """Retrieves a single sequence from a FASTA file.

    In the worst case, this requires iteration over the entire file at once, so
    please don't call it once for every contig in a for loop or something.
    The goal is to just retrieve a single sequence relatively quickly, without
    having to store every sequence in the file in memory or something.

    Parameters
    ----------
    fasta_fp: str
        Filepath to a FASTA file.

    contig_name: str
        Name of a contig. This contig's sequence will be retrieved from the
        FASTA file.

    Returns
    -------
    contig_seq: skbio.DNA
        Sequence of this contig.

    Raises
    ------
    SequencingDataError
        If this contig isn't present in the FASTA file.
    """
    for seq in skbio.io.read(fasta_fp, format="fasta", constructor=skbio.DNA):
        if seq.metadata["id"] == contig_name:
            return seq
    raise SequencingDataError(f"Contig {contig_name} is not in {fasta_fp}.")
