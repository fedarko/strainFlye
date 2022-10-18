# Utilities for strainFlye matrix.


import skbio
from collections import defaultdict
from strainflye import cli_utils, misc_utils, gff_utils, config


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
    for ci, (contig, im) in enumerate(contig_and_im_tuples, 1):

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
            continue

        verboselog(
            f"Found {fs} belonging to contig {contig}; inspecting...",
            prefix="",
        )

        seen_feature_ids = set()
        # TODO: It would probably make sense to construct a DataFrame from
        # these genes. This would let us use vectorization to compute
        # genes_overlapping_aln.
        #
        # I think we could probs also do this with just numpy arrays, but it
        # seems tricky to relate feature bounds back to the corresponding
        # IDs/strands/etc., and I already have code from the analysis notebooks
        # that uses DataFrames lol
        fid2codon2alignedcodons = {}
        fid2range = {}
        fid2strand = {}
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
            is_cds, strand = gff_utils.validate_if_cds(
                feature, contig, verboselog
            )
            if is_cds:
                fid2range[fid] = feature_range
                fid2strand[fid] = strand
                fid2codon2alignedcodons[fid] = {}
                for cp_left in range(
                    feature_range[0], feature_range[-1] + 1, 3
                ):
                    fid2codon2alignedcodons[fid][cp_left] = defaultdict(int)

        num_cds = len(fid2range)
        feature_noun = "feature" if num_cds == 1 else "features"
        fs = (
            f"Found {num_cds:,} {feature_noun} with a type in "
            f"{config.CDS_TYPES} in contig {contig}."
        )
        # This case is a bit weird, so it merits being loud about
        if num_cds == 0:
            fancylog(f"{fs} Ignoring this contig.", prefix="")
            break
        verboselog(f"{fs} Going through alignments...", prefix="")

        # OK, now we can go through the alignment file and figure out which
        # alignments span which codons. This'll let us "finalize"
        # fid2codon2alignedcodons.
        # After that, we can output that to a file for this contig.

    fancylog("Done.", prefix="")
