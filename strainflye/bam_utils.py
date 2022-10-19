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


def get_coords(aln):
    """Returns a linear alignment's coordinates.

    Parameters
    ----------
    alnseg: pysam.AlignedSegment

    Returns
    -------
    (s, e): (int, int, str)
        Segment start and end.

        The start and end are both inclusive, to simplify comparison of
        alignment ranges for detecting overlaps.

    Raises
    ------
    WeirdError
        If the segment's start is greater than its end (both in inclusive
        coordinates). (If this ends up being a problem in practice, maybe
        because there of reverse-mapped reads or something (???), then this
        could probs be modified to just reverse the start and end in this
        case.)
    """
    # We add 1 to the end since this is a half-open interval -- we want
    # the coordinates we use for computing overlap to be completely
    # inclusive intervals. TODO NOPE FIX
    s = aln.reference_start
    e = aln.reference_end + 1
    if s > e:
        raise WeirdError(
            f"Malformed linear alignment coordinates: start {s}, end {e}"
        )
    return (s, e)
