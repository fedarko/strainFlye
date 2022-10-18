# Utilities for strainFlye matrix.


import skbio
from strainflye import cli_utils, misc_utils, config


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

    fancylog("Going through contigs and counting aligned 3-mers...")

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
            f"Found {fs} belonging to contig {contig}...",
            prefix="",
        )

        num_gene_features_seen = 0
        # gene2codon2alignedcodons = {}
        # This "hack" to go through all features in "im" is taken from
        # spot_utils.run_hotspot_feature_detection() -- see that function for
        # context.
        for feature in im.query(metadata={}):
            # Ignore features that don't aren't labelled as "CDS"s
            # As far as I can tell, this is the only type (or at least the main
            # one) that is used to identify coding sequences in GFF3 files.
            # "Gene" is a possible type, for example, but it seems mostly to be
            # used as a higher-level type for top-level genes containing
            # multiple CDSs -- for reference, Prodigal's output seems to
            # consist solely of CDSs.
            if feature.metadata["type"].upper() not in config.GFF_CDS_TYPES:
                verboselog(
                    (
                        f"Ignoring feature {feature.metadata['ID']}, since "
                        "its type isn't one of "
                        f"{cli_utils.list2str(config.GFF_CDS_TYPES)}."
                    ),
                    prefix="",
                )
                continue
            num_gene_features_seen += 1

            # TODO: read phase info -- if needed, we can just mandate that it
            # be zero and raise an error if not. but i think accounting for it
            # wouldn't be too bad
            #
            # Then, figure out strand.
            #
            # Then we can figure out codons in this gene. Update
            # gene2codon2alignedcodons.

        # OK, now we can go through the alignment file and figure out which
        # alignments span which codons. This'll let us "finalize"
        # gene2codon2alignedcodons.
        #
        # After that, we can output that to a file for this contig.

    fancylog("Done.", prefix="")
