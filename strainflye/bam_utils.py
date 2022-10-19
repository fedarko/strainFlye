# Utilities for dealing with BAM files.

import os
import subprocess
from strainflye.errors import WeirdError


def index_bam(in_bam, bam_descriptor, fancylog):
    """Indexes a BAM file using samtools.

    This creates a .bam.bai file in the same location as the BAM file.

    Parameters
    ----------
    in_bam: str
        Location of the BAM file to be indexed.

    bam_descriptor: str
        Description of the in_bam file, to be used in the logging message.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    fancylog(f"Indexing the {bam_descriptor}...")
    subprocess.run(["samtools", "index", in_bam], check=True)
    fancylog(f"Done indexing the {bam_descriptor}.", prefix="")


def rm_bam(bam_fp, bam_descriptor, do_removal, fancylog):
    """Removes a temporary BAM file and its index.

    Parameters
    ----------
    bam_fp: str
        Location of the BAM file to be removed. We assume another file in the
        same directory (with the same name, and ending in .bai) exists; we'll
        remove this "index" file as well.

    bam_descriptor: str
        Description of the bam_fp file, to be used in the logging message.

    do_removal: bool
        If True, actually remove the BAM file and its index; if False, don't do
        anything.

    fancylog: function
        Logging function.

    Returns
    -------
    None
    """
    if do_removal:
        fancylog(
            f"Removing the {bam_descriptor} and its index to save space..."
        )
        os.remove(bam_fp)
        os.remove(bam_fp + ".bai")
        fancylog("Done removing these files.", prefix="")


def get_coords(aln, zero_indexed=True):
    """Returns and checks a linear alignment's coordinates on the reference.

    Parameters
    ----------
    aln: pysam.AlignedSegment
        Representation of a linear alignment. We don't say "read" because a
        single read can have multiple linear alignments, since we allow
        supplementary alignments.

    zero_indexed: bool
        If True, return 0-indexed coordinates; if False, return 1-indexed
        coordinates.

    Returns
    -------
    (ref_start, ref_end): (int, int)
        Segment start and end. These are inclusive coordinates, so this linear
        alignment should span the range [ref_start, ref_end].

    Raises
    ------
    WeirdError
        If either the start or end are None (which is technically possible
        due to weird phenomena like unmapped reads, per the pysam docs).

        If the start is greater than the end (both in inclusive coordinates).
        We allow start == end in the rare case that this alignment is exactly
        one nucleotide long.
    """
    # These are 0-indexed positions -- so e really points to the rightmost
    # position plus one.
    s = aln.reference_start
    e = aln.reference_end

    # However, they could be None -- check for this and error out if so.
    if s is None or e is None:
        raise WeirdError(f"Alignment {aln.query_name} is unmapped?")

    # Now that we know that these are not None, we can assume they're numbers.
    if zero_indexed:
        # Subtract 1 from the end to make this range inclusive.
        e -= 1
    else:
        # If we are converting from 0- to 1-indexing, then e is already "good"
        # -- we just need to add 1 to s.
        s += 1

    # In any case, now that these coordinates define an inclusive range, we can
    # verify that they represent a valid start and end.
    if s > e:
        idx = "0" if zero_indexed else "1"
        raise WeirdError(
            f"Malformed linear alignment coordinates for {aln.query_name}: "
            f"start {s:,} is > end {e:,} (both {idx}-indexed and inclusive)"
        )
    return (s, e)
