# General IO utilities.


import skbio
from strainflye.errors import SequencingDataError


APOLOGY = "This isn't supported at the moment, sorry."


def get_name2len(fasta_fp):
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

    return name2len
