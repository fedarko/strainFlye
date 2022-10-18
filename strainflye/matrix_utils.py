# Utilities for strainFlye matrix.


import skbio
import pandas as pd
from collections import defaultdict
from strainflye import cli_utils, misc_utils, gff_utils, config


def get_contig_cds_info(im, contig, contig_name2len, fancylog, verboselog):
    num_features = im.num_interval_features
    feature_noun = "feature" if num_features == 1 else "features"
    fs = f"{num_features:,} {feature_noun}"

    # Ignore "extra" contigs in the GFF3 but not the FASTA
    if contig not in contig_name2len:
        verboselog(
            (
                f"Found {fs} belonging to sequence {contig} in the GFF3 "
                "file. Ignoring them, since this sequence isn't in the "
                "FASTA file."
            ),
            prefix="",
        )
        return None, None

    verboselog(
        f"Found {fs} belonging to contig {contig}; inspecting...",
        prefix="",
    )

    # Set of feature IDs seen in this contig -- will be updated and checked
    # to make sure that no features in this contig have duplicate IDs
    seen_feature_ids = set()

    # Build up lists of information about CDS features: IDs, coordinates,
    # and strands. We'll then convert this information into a DataFrame
    # later.
    feature_ids = []
    feature_left_ends = []
    feature_right_ends = []
    feature_strands = []

    fid2codon2alignedcodons = {}

    # This "hack" to go through all features in "im" is taken from
    # spot_utils.run_hotspot_feature_detection() -- see that function for
    # context.
    for feature in im.query(metadata={}):
        fid, feature_range = gff_utils.validate_basic(
            feature,
            contig,
            contig_name2len[contig],
            seen_feature_ids,
            fancylog,
            zero_indexed_range=False,
        )
        is_cds, strand = gff_utils.validate_if_cds(feature, contig, verboselog)
        if is_cds:
            feature_ids.append(fid)
            feature_left_ends.append(feature_range[0])
            feature_right_ends.append(feature_range[-1])
            feature_strands.append(strand)
            fid2codon2alignedcodons[fid] = {}
            for cp_left in range(feature_range[0], feature_range[-1] + 1, 3):
                fid2codon2alignedcodons[fid][cp_left] = defaultdict(int)

    num_cds = len(feature_ids)
    feature_noun = "feature" if num_cds == 1 else "features"
    fs = (
        f"Found {num_cds:,} {feature_noun} with a type in "
        f"{config.CDS_TYPES} in contig {contig}."
    )
    # This case is a bit weird, so it merits being loud about
    if num_cds == 0:
        fancylog(f"{fs} Ignoring this contig.", prefix="")
        return None, None
    verboselog(f"{fs} Going through alignments...", prefix="")

    # Construct a DataFrame of gene coordinate information. This makes it
    # easy to use vectorization to figure out which alignments overlap
    # which genes. I think we could probs also do something similar with
    # just numpy arrays, but it seems tricky to relate feature bounds back
    # to the corresponding IDs/strands/etc., and I already have code from
    # the analysis notebooks that assumes that genes are available in a
    # DataFrame...
    cds_df = pd.DataFrame(
        {
            "LeftEnd": feature_left_ends,
            "RightEnd": feature_right_ends,
            "Strand": feature_strands,
        },
        index=feature_ids,
    )
    return cds_df, fid2codon2alignedcodons


def run_count(contigs, bam, genes, output_dir, verbose, fancylog):
    """Counts the 3-mers aligned to each predicted gene in each contig.

    For each contig, writes out a "pickle" file to the output directory.

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs.

    bam: str
        Filepath to a BAM file mapping reads to contigs.

    genes: str
        Filepath to a GFF3 file describing gene coordinates in contigs.

    output_dir: str
        Directory to which we'll write out information for each contig.

    verbose: bool
        Log extra info.

    fancylog: function
        Logging function.

    Returns
    -------
    None

    Raises
    ------
    Various errors are raised by misc_utils.load_triplet() if the input files
    are problematic.
    """
    verboselog = cli_utils.get_verboselog(fancylog, verbose)

    # Load and check the FASTA and BAM files.
    contig_name2len, bam_obj, num_contigs = misc_utils.load_fasta_and_bam(
        contigs, bam, fancylog
    )

    fancylog("Counting aligned 3-mers to coding sequences in contigs...")

    # See spot_utils.run_hotspot_feature_detection() -- same idea here.
    contig_and_im_tuples = skbio.io.read(genes, format="gff3")
    for contig, im in contig_and_im_tuples:

        cds_df, fid2codon2alignedcodons = get_contig_cds_info(
            im, contig, contig_name2len, fancylog, verboselog
        )

        # This DataFrame will be None if we couldn't get any CDS information
        # for this contig -- maybe this contig was an "extra" contig in the
        # GFF3 but not in the FASTA, or maybe this contig was in the FASTA but
        # just didn't have any CDS features in the GFF3 file. In either case,
        # we won't create any output for this contig in particular.
        if cds_df is None:
            continue

        # ... But if we've made it here, then we know that this contig has at
        # least one CDS feature, so we can count 3-mers here.

    fancylog("Done.", prefix="")
