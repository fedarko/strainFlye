# Utilities for strainFlye matrix.


import skbio
import pandas as pd
from collections import defaultdict
from strainflye import cli_utils, misc_utils, bam_utils, gff_utils, config


def get_contig_cds_info(im, contig, contig_name2len, fancylog, verboselog):
    """Returns information about the CDS feature(s) in a contig.

    Parameters
    ----------
    im: skbio.metadata.IntervalMetadata
        Object describing the feature(s) in a GFF3 file for this contig.

    contig: str
        Name of the contig that these feature(s) belong to, according to
        the GFF3 file from which this information was parsed.

    contig_name2len: dict
        Maps contig names to lengths; this should have been obtained from
        parsing a FASTA file. There is no guarantee that "contig",
        specifically, is present within this dict as a key.

    fancylog: function
        Logging function. We'll use this to log particularly unusual things
        that deserve the user's attention even if --verbose was not specified.

    verboselog: function
        Logging function. We'll use this to log about minor details.

    Returns
    -------
    (cds_df, fid2codon2alignedcodons): (pd.DataFrame, dict)

        If no CDS features belong to "contig" -- and/or if "contig" is not in
        contig_name2len -- then both cds_df and fid2codon2alignedcodons will be
        None. Otherwise:

        cds_df describes all of the CDS features in "im". Having this
        information in a DataFrame, in particular, is useful because we can
        perform vectorized operations on it to do things like figure out which
        alignments overlap which genes. (I think we could probs also do
        something similar with just numpy arrays, but it seems tricky to relate
        feature bounds back to the corresponding IDs/strands/etc., and I
        already have code from the analysis notebooks that assumes that genes
        are available in a DataFrame...) Anyway, this DataFrame has one row per
        feature. Each row has an index (the feature ID) and LeftEnd, RightEnd,
        and Strand values.

        fid2codon2alignedcodons maps CDS IDs (the indices in cds_df) --> the
        leftmost position of all codons in this CDS --> an empty
        defaultdict(int). This innermost defaultdict(int) can be used to count
        how many times a given 3-mer is aligned to this codon. Note that we say
        "leftmost" position regardless of strand -- we'll figure out
        reverse-complementing stuff later.

    Raises
    ------
    ParameterError
        Raised by gff_utils.validate_basic() or gff_utils.validate_if_cds() if
        various things about the features in "im" are invalid. See those
        functions for details.
    """
    num_features = im.num_interval_features
    feature_noun = "feature" if num_features == 1 else "features"
    fs = f"{num_features:,} {feature_noun}"

    # Ignore "extra" contigs in the GFF3 but not the FASTA
    if contig not in contig_name2len:
        verboselog(
            (
                f"Found {fs} belonging to sequence {contig} in the GFF3 "
                f"file. Ignoring this sequence and its {feature_noun}, "
                f"since {contig} isn't in the FASTA file."
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
        # If this feature is a CDS feature, do extra validation on it -- make
        # sure it has a strand of + or -, that it has a length divisible by 3,
        # etc.
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
    # This case is a bit weird, so it merits being loud about. Note that us
    # having reached this point implies that this contig *did* have other
    # feature(s) in this GFF3 file, so this case won't trigger for contigs that
    # simply just have no features in the GFF3 file -- it'll only trigger for
    # contigs that have features, but no CDS features in particular.
    if num_cds == 0:
        fancylog(f"{fs} Ignoring this contig.", prefix="")
        return None, None
    verboselog(f"{fs} Going through alignments...", prefix="")

    cds_df = pd.DataFrame(
        {
            "LeftEnd": feature_left_ends,
            "RightEnd": feature_right_ends,
            "Strand": feature_strands,
        },
        index=feature_ids,
    )
    return cds_df, fid2codon2alignedcodons


def count_aligned_3mers(bam_obj, contig, cds_df, fid2codon2alignedcodons):
    """Counts 3-mers aligned to codons in CDS features in a contig.

    This function doesn't return anything, but it will fill in
    fid2codon2alignedcodons with this information.

    Parameters
    ----------
    bam_obj: pysam.AlignmentFile
        Object describing a BAM file mapping reads to contigs.

    contig: str
        Name of a contig. We assume that this contig already exists as a
        reference in bam_obj.

    cds_df: pd.DataFrame
        Contains information about coordinates and strands of coding sequences
        within this contig. Can be produced by get_contig_cds_info().

    fid2codon2alignedcodons: dict
        The other output from get_contig_cds_info(). This maps CDS IDs (the
        indices in cds_df) --> the leftmost position of all codons in this CDS
        --> an empty defaultdict(int). These empty defaultdicts will be updated
        when this function is run.

    Returns
    -------
    None
    """
    for aln in bam_obj.fetch(contig):
        # Find all genes that this aln intersects in this contig

        # Get one-indexed and inclusive coordinates
        ref_start, ref_end = bam_utils.get_coords(aln, zero_indexed=False)

        # Use vectorization to find genes overlapping this aln: see
        # https://stackoverflow.com/a/17071908
        # for details on why parentheses, etc., and
        # https://engineering.upside.com/a-beginners-guide-to-optimizing-pandas-code-for-speed-c09ef2c6a4d6
        # for justification on why this is useful (tldr: makes code go fast)
        genes_overlapping_aln = list(
            cds_df.loc[
                (cds_df["RightEnd"] >= ref_start)
                & (cds_df["LeftEnd"] <= ref_end)
            ].itertuples()
        )
        # Note about the above thing: you may be shaking your fist and saying
        # "wait itertuples is slow!" And yeah, kinda. But for whatever reason
        # I've tried multiple times to keep genes_overlapping_aln as a
        # DataFrame (and then later vectorize stuff like checking that a given
        # aligned pair covers a codon within the genes, etc) and the overhead
        # costs seem to slow things down. I am sure it's possible to speed
        # things up more, but right now things seem good enough.

        # If no genes overlap this aln, we are free to move on to the next aln.
        if len(genes_overlapping_aln) > 0:

            # Computing this is relatively slow, which is why we jump through
            # so many hoops before we do this. Each entry in
            # get_aligned_pairs() is a tuple with 2 elements:
            # the first is the query/read pos and the second is the reference
            # pos. TODO: would it be possible to only do this for certain
            # positions we care about? get_aligned_pairs() returns a lot of
            # stuff we don't need, e.g. regions of the aln that don't intersect
            # with any genes.
            ap = aln.get_aligned_pairs(matches_only=True)

            # Doesn't look like getting this in advance saves much time, but I
            # don't think it hurts.
            read_seq = aln.query_sequence

            # We only consider the leftmost position of each codon, so we don't
            # need to bother checking the last two pairs of positions (since
            # neither could be the leftmost position of a codon that this aln
            # fully covers).
            for api, pair1 in enumerate(ap[:-2]):

                # Convert to 1-indexed position for ease of comparison with
                # gene coordinates
                pair1_refpos = pair1[1] + 1

                havent_checked_next_pairs = True
                for gene_data in genes_overlapping_aln:
                    gl = gene_data.LeftEnd
                    gr = gene_data.RightEnd

                    # Check that this pair is located within this gene and is
                    # the leftmost position of a codon in the gene. (Note that
                    # this check works for both + or - strand genes. Whether
                    # the leftmost position is the "start" [i.e. CP 1] or "end"
                    # [i.e. CP 3] of the gene changes with the strand of the
                    # gene, but we'll account for that later on when we
                    # reverse-complement the codon if needed.)
                    if (
                        pair1_refpos >= gl
                        and pair1_refpos <= gr - 2
                        and ((pair1_refpos - gl) % 3 == 0)
                    ):

                        # Nice! Looks like this aln fully covers this codon.

                        # If we haven't yet, check that this aln doesn't skip
                        # over parts of the codon, or stuff like that. The
                        # reason this check is located *here* (and not before
                        # we loop over the genes) is that it seems like this
                        # is a faster strategy: only run these checks once we
                        # KNOW that this pair looks like it fully covers a
                        # codon, since many pairs might not meet that criteria.
                        #
                        # (And by recording that we've run this check once, in
                        # havent_checked_next_pairs, we can save the time cost
                        # of running the check multiple times.)
                        #
                        # "I feel like an insane person trying to optimize this
                        # so much lmao" -- Me when I originally wrote this
                        # code in like April 2021 (it is now October 2022 and I
                        # am getting a little tired of this project and also my
                        # soul hurts)
                        if havent_checked_next_pairs:
                            # Check that the pairs are all consecutive (i.e. no
                            # "jumps" in the read, and no "jumps" in the
                            # reference). Since we don't consider the last two
                            # pairs in ap, pair2 and pair3 should always be
                            # available.
                            pair2 = ap[api + 1]
                            pair3 = ap[api + 2]

                            # Ensure that the read positions are consecutive
                            p10 = pair1[0]
                            p20 = pair2[0]
                            p30 = pair3[0]
                            readpos_consec = ((p10 + 1) == p20) and (
                                (p20 + 1) == p30
                            )
                            if not readpos_consec:
                                break

                            # Ensure that the ref. positions are consecutive
                            p11 = pair1[1]
                            p21 = pair2[1]
                            p31 = pair3[1]
                            refpos_consec = ((p11 + 1) == p21) and (
                                (p21 + 1) == p31
                            )
                            if not refpos_consec:
                                break

                            havent_checked_next_pairs = False

                        # Figure out what the read actually *says* in the
                        # alignment here. (It'll probably be a complete match
                        # most of the time, but there will be some occasional
                        # mismatches -- and seeing those is ... the whole point
                        # of this notebook.)

                        # We make sure to index the read by read coords, not
                        # reference coords!
                        aligned_codon = read_seq[p10 : p10 + 3]

                        # Finally, update information about codon counts.
                        gi = gene_data.Index
                        if gene_data.Strand == "-":
                            true_aligned_codon = config.CODON2RC[aligned_codon]
                        else:
                            true_aligned_codon = aligned_codon

                        fid2codon2alignedcodons[gi][pair1_refpos][
                            true_aligned_codon
                        ] += 1


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

    misc_utils.make_output_dir(output_dir)

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

        # ... So if we've made it here, then we know that this contig has at
        # least one CDS feature, and we can count 3-mers in these feature(s).
        count_aligned_3mers(bam_obj, contig, cds_df, fid2codon2alignedcodons)

        # Now that we've finished this counting operation, write out the 3-mer
        # count information to a file.
        misc_utils.write_obj_to_pickle(
            fid2codon2alignedcodons, output_dir, contig, "3mers"
        )

    fancylog("Done.", prefix="")
