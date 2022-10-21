# Utilities for strainFlye matrix.


import os
import skbio
import pandas as pd
from collections import defaultdict
from strainflye import (
    cli_utils,
    misc_utils,
    bam_utils,
    gff_utils,
    fasta_utils,
    config,
)
from strainflye.errors import WeirdError


class CodonCounter(object):
    """Represents counts of 3-mers aligned to codons in a contig.

    It's worth noting that this is a fairly "naive" implementation of this. We
    store codons as strings, for example. It should be possible to convert
    codons to 0-based indices into config.CODONS, and thus compress this, but
    I am erring on the side of simplicity for now.
    """

    def __init__(self, contig, contig_seq):
        """Initializes the object.

        Parameters
        ----------
        contig: str
            Name of the contig.

        contig_seq: skbio.DNA
            Sequence of the contig. Including this information here is one of
            the main reasons for having this be its own class -- this is a lot
            more convenient than making the user pass the contigs FASTA file to
            "strainFlye matrix fill" again, double-checking that things are
            consistent, etc. (We could try to save space by only storing the
            sequence for codons in CDSs, but that's an optimization for another
            day. Also, for densely packed overlapping genes on a contig, that
            might actually increase storage space without using some tricks.)
        """
        self.contig = contig
        self.contig_seq = contig_seq
        # The main object we use to store counts
        self.cds2left2counter = {}
        self.cds2strand = {}

    def __str__(self):
        num_cds = len(self.cds2left2counter)
        clen = len(self.contig_seq)
        return f"CodonCounter({self.contig}, {clen:,} bp, {num_cds:,} CDSs)"

    def add_cds(self, cds_id, cds_left, cds_right, strand):
        # this is overkill because we should've already validated this feature
        # using gff_utils, but better safe than sorry
        if (cds_right - cds_left + 1) % 3 != 0:
            raise WeirdError(
                f"{cds_id} has left = {cds_left:,}, right = {cds_right:,}: "
                "not divisible by 3?"
            )
        if strand != "+" and strand != "-":
            raise WeirdError(f"{cds_id} has strand = {strand}?")

        self.cds2left2counter[cds_id] = {}
        for cp_left in range(cds_left, cds_right + 1, 3):
            self.cds2left2counter[cds_id][cp_left] = defaultdict(int)

        self.cds2strand[cds_id] = strand

    def add_count(self, cds_id, cp_left, raw_aligned_codon):
        # We'll have to reverse-complement the codons aligned to - strand genes
        # eventually, so we might as well do it here.
        if self.cds2strand[cds_id] == "-":
            true_aligned_codon = config.CODON2RC[raw_aligned_codon]
        else:
            true_aligned_codon = raw_aligned_codon

        self.cds2left2counter[cds_id][cp_left][true_aligned_codon] += 1


def get_contig_cds_info(im, contig, contig_seq, fancylog, verboselog):
    """Returns information about the CDS feature(s) in a contig.

    Parameters
    ----------
    im: skbio.metadata.IntervalMetadata
        Object describing the feature(s) in a GFF3 file for this contig.

    contig: str
        Name of the contig that these feature(s) belong to, according to
        the GFF3 file from which this information was parsed.

    contig_seq: skbio.DNA
        Sequence of this contig.

    fancylog: function
        Logging function. We'll use this to log particularly unusual things
        that deserve the user's attention even if --verbose was not specified.

    verboselog: function
        Logging function. We'll use this to log about minor details.

    Returns
    -------
    (cds_df, cc): (pd.DataFrame, CodonCounter)

        If no CDS features belong to "contig", then both cds_df and cc will be
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

        cc is an object that can be used to count aligned 3-mers to the codons
        in the CDSs in a contig. (Although at this point it is "empty," because
        we have not considered the alignment yet.)

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

    # I assume that, if a contig is described in the GFF3 file (to the point
    # where scikit-bio's parser returns a pair of (contig, IntervalMetadata),
    # then there is at least one feature for this contig in the GFF3 file. (If
    # a contig has zero features total but it still shows up here for some
    # bizarre reason, then the zero-CDS-features case should still be hit.)
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

    cc = CodonCounter(contig, contig_seq)

    # This "hack" to go through all features in "im" is taken from
    # spot_utils.run_hotspot_feature_detection() -- see that function for
    # context.
    for feature in im.query(metadata={}):
        fid, feature_range = gff_utils.validate_basic(
            feature,
            contig,
            len(contig_seq),
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
            left = feature_range[0]
            right = feature_range[-1]
            feature_left_ends.append(left)
            feature_right_ends.append(right)
            feature_strands.append(strand)
            cc.add_cds(fid, left, right, strand)

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
    return cds_df, cc


def count_aligned_3mers(bam_obj, contig, cds_df, cc):
    """Counts 3-mers aligned to codons in CDS features in a contig.

    This function doesn't return anything, but it will fill in cc with this
    information.

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

    cc: CodonCounter
        The other output from get_contig_cds_info(). This stores count
        information of 3-mers aligned to the codons of the CDSs of this contig.
        We'll update this object.

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
                        # of this analysis.)

                        # We make sure to index the read by read coords, not
                        # reference coords!
                        aligned_codon = read_seq[p10 : p10 + 3]
                        # Finally, count this aligned codon! (or 3-mer,
                        # depending on your terminology preferences)
                        cc.add_count(
                            gene_data.Index, pair1_refpos, aligned_codon
                        )


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
        Directory to which we'll write out count information for each contig.

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

        # Ignore "extra" contigs in the GFF3 but not the FASTA
        if contig not in contig_name2len:
            verboselog(
                (
                    f"Contig {contig} is in the GFF3 file but not the FASTA "
                    "file; ignoring it."
                ),
                prefix="",
            )
            continue

        # NOTE: Using get_single_seq() like this is inefficient, but it gets
        # the job done; see notes in smooth_utils.write_smoothed_reads().
        contig_seq = fasta_utils.get_single_seq(contigs, contig)

        cds_df, cc = get_contig_cds_info(
            im, contig, contig_seq, fancylog, verboselog
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
        count_aligned_3mers(bam_obj, contig, cds_df, cc)

        # Now that we've finished counting, write out the 3-mer count
        # information to a file.
        misc_utils.write_obj_to_pickle(
            cc, output_dir, contig, config.CT_FILE_LBL
        )
        verboselog(
            f"Wrote out 3-mer count info for contig {contig}.",
            prefix="",
        )

    fancylog("Done.", prefix="")


def init_matrix_structs():
    # 64x63 dict: each key is a codon string, and each value is
    # another dict with all the other codons
    codon2codon2ct = {c: defaultdict(int) for c in config.CODONS}

    # 21x20 dict: each key is a proteinogenic amino acid (A, C, D,
    # E, F, ...), limited to just stuff in the standard genetic code
    # (i.e. ignoring selenocystine and pyrrolsine) but including "*",
    # representing a stop codon.
    aa2aa2ct = {a: defaultdict(int) for a in config.AAS}

    # 64-key dict: maps each codon to an integer indicating how
    # frequently this codon occurs in all CDSs in the "reference"
    # contig (i.e. not counting mutations into this codon).
    codon2ct = defaultdict(int)

    # 21-key dict: maps amino acid/stop codon to integer indicating
    # frequency across all CDSs in this contig.
    aa2ct = defaultdict(int)

    return codon2codon2ct, aa2aa2ct, codon2ct, aa2ct


def run_fill(
    ct_dir, p, min_alt_pos, r, output_format, output_dir, verbose, fancylog
):
    """Calls codon mutations and creates mutation matrices.

    For each contig, creates a subdirectory within output_dir containing
    codon2codon2ct, aa2aa2ct, codon2ct, and aa2ct files. The x2x2ct files
    represent mutation matrices, and the x2ct files can be useful in
    visualizing these.

    Parameters
    ----------
    count_dir: str
        Directory containing 3-mer count information from run_count().

    p: int >= 1 or None
        If specified, is a value of p used for calling p-mutations.

    min_alt_pos: int
        During p-mutation calling, the second-most-common aligned nucleotide's
        frequency must be at least this to call a p-mutation at a position. Not
        used when calling r-mutations.

    r: int >= 1 or None
        If specified, is a value of r used for calling r-mutations.

    output_format: str
        Should be one of {"tsv", "json"}.

    output_dir: str
        Directory to which we'll write out matrix information for each contig.

    verbose: bool
        Log extra info.

    fancylog: function
        Logging function.

    Returns
    -------
    None

    Raises
    ------
    TBD
    """
    verboselog = cli_utils.get_verboselog(fancylog, verbose)

    ctf_suffix = f"_{config.CT_FILE_LBL}.pickle"

    fancylog("Going through aligned 3-mer counts and creating matrices...")
    for fp in sorted(os.listdir(ct_dir)):
        if fp.lower().endswith(ctf_suffix):

            # Load the CodonCounter object for this contig
            cc = misc_utils.load_from_pickle(os.path.join(ct_dir, fp))
            verboselog(cc)

            c2c2ct, c2ct, a2a2ct, a2ct = init_matrix_structs()
