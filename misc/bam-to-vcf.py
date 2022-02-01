#! /usr/bin/env python3
#
# NOTE: This is a very early version of this code that just does naive variant
# calling without any false discovery rate control. This code is an
# amalgamation of multiple analysis scripts I've written, and will be replaced
# with an actual easier-to-use pipeline soon.
#
# Converts a BAM file to VCF using naive variant calling.
#
# Dependencies: Python 3.6+ (uses f-strings), pysam, pysamstats, scikit-bio

import os
import time
import skbio
import pysam
import pysamstats

# Script configuration settings:
###############################################################################
# 1. SEQS: Which MAGs / contigs should we use for variant calling?
#          These should have FASTA files located in FASTA_DIR. Each name should
#          match the name of a sequence in the BAM file.
SEQS = ["edge_6104", "edge_1671", "edge_2358"]

# 2. FASTA_DIR: Directory containing one FASTA file per sequence described in
#               SEQ. The name of each FASTA file should match the name listed
#               in SEQS (i.e. for the current data, the FASTA_DIR directory
#               contains three FASTA files -- named "edge_6104.fasta",
#               "edge_1671.fasta", etc.)
FASTA_DIR = "seqs/"

# 3. BAM_FILE: Filepath pointing to a BAM file. All of the sequences listed in
#              SEQS should be included as "reference" sequences here.
BAM_FILE = "main-workflow/output/fully-filtered-and-sorted-aln.bam"

# 4. P_THRESHOLD: Number in the range (0, 50] describing the percentage we will
#       use to naively call SNPs. For example, p = 40 means that we'll only
#       call SNPs if freq(pos) (described in the paper) is >= 40%.
#
#       Right now, you give this percentage and we use it to
#       output the SNP calls, but soon this will be reversed (you give a FDR
#       and we use that to determine the value of P used to call SNPs).
#
#       A "good" setting for this value will depend on how high-coverage your
#       MAGs are (again, this will eventually be set automatically); a good
#       rule-of-thumb lower limit might be P = 100 * (5 / average coverage).
P_THRESHOLD = 1
assert P_THRESHOLD <= 50 and P_THRESHOLD > 0

# 5. MIN_READ_NUMBER: Used (along with p) to determine the minimum coverage of
#                     a position we can call SNPs at. The default of 5 is
#                     probably good for most datasets.
#
MIN_READ_NUMBER = 5
assert MIN_READ_NUMBER >= 0

# 6. MIN_ALT_POS: minimum count of alternate nucleotides required to call a
#                 SNP, regardless of the percentages. See the
#                 "if alt(pos) > 1 [...]" section in the paper for explanation.
MIN_ALT_POS = 2

# 7. VCF_FILE: Output filepath for SNP calls
VCF_FILE = "main-workflow/output/strainflye-naive-calls.vcf"

###############################################################################
# Other settings that probably don't need changing
###############################################################################
MIN_SUFFICIENT_COV = (100 * MIN_READ_NUMBER) / P_THRESHOLD

# How often (i.e. after processing how many positions of a given sequence)
# should we update the user about what's happening? (This won't impact the
# actual VCF file, it'll just impact how fast text is shown on the terminal
# while running this script.)
UPDATE_FREQ = 5000

# Max depth at any position in the alignment: needed for pysam, see pysamstats
# docs
MAX_DEPTH_PYSAM = 1000000
###############################################################################

print(f"Using seqs {SEQS}.")

bf = pysam.AlignmentFile(BAM_FILE)

# Header info gleaned by reading over the VCF 4.2 docs
# (https://samtools.github.io/hts-specs/VCFv4.2.pdf) and copying how LoFreq
# organizes their header
vcf = (
    "##fileformat=VCFv4.2\n"
    f"##fileDate={time.strftime('%Y%m%d')}\n"
    # TODO include version number in the source, as shown in VCF 4.2 docs?
    f'##source="strainFlye (naive calling: p = {P_THRESHOLD}, min alt pos = '
    f'{MIN_ALT_POS}, min read number = {MIN_READ_NUMBER})"\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)
sep = "-" * 79
print(
    f"{sep}\nFor reference, the VCF header of the resulting file will be:\n"
    f"{vcf}{sep}"
)

for seq in SEQS:
    fastapath = os.path.join(FASTA_DIR, f"{seq}.fasta")
    fasta = str(skbio.DNA.read(fastapath))
    seqlen = len(fasta)
    print(f"Seq {seq} has length {seqlen:,} bp.")

    num_variants_seen_so_far = 0

    # We use start=0 and end=(sequence length) because the start and end params
    # of pysamstats.stat_variation() are 0-indexed (although the normal
    # samtools pileup format's coordinates are 1-indexed, and although our
    # output will be 1-indexed).
    # We use pad=True so that even uncovered positions are included -- for
    # example, edge 1671 has some regions in the middle that are uncovered
    # after read filtering, at least as of writing.
    for pos, rec in enumerate(
        pysamstats.stat_variation(
            bf,
            chrom=seq,
            fafile=fastapath,
            start=0,
            end=seqlen,
            truncate=True,
            max_depth=MAX_DEPTH_PYSAM,
            pad=True,
        ),
        1,
    ):
        if rec["N"] > 0:
            # we can definitely handle this case, but it'll take a bit of
            # thinking; there are multiple ways we could do this, and for now
            # we avoid choosing
            raise ValueError("There shouldn't be any Ns in this data?")

        rpos = rec["pos"] + 1
        if rpos != pos:
            raise ValueError(
                f"Found discontinuity in traversal: {pos}-th pos, but "
                f"rec['pos'] + 1 is {rpos}"
            )

        ref_nt = rec["ref"]
        # Paranoid checking
        if ref_nt != fasta[pos - 1]:
            raise ValueError(
                "Looks like FASTA and pysamstats disagree on ref nt at pos "
                f"{pos - 1} (0-indexed), a.k.a. {pos} (1-indexed). "
                f"FASTA says it's {fasta[pos - 1]}; pysamstats says it's "
                f"{ref_nt}."
            )

        # from notebooks/pleuk_copied_code.py
        cov = rec["A"] + rec["C"] + rec["G"] + rec["T"]

        is_variant = False
        if cov >= MIN_SUFFICIENT_COV:
            ordered_nts = sorted("ACGT", key=rec.get)

            # The "reference" nucleotide we'll list in the VCF. This is the
            # maximum-frequency nucleotide in the pileup at this position
            # (breaking ties arbitrarily). This overrides the actual reference
            # in the contig sequence, for weird positions where the contig
            # sequence isn't also the alignment consensus.
            ref_nt = ordered_nts[-1]

            # The literal nucleotide used in the numerator of freq(pos): one
            # of A, C, G, T
            alt_nt = ordered_nts[-2]

            # The raw frequency (in counts) of alt_nt. An integer >= 0.
            alt_nt_freq = rec[alt_nt]

            if alt_nt_freq >= MIN_ALT_POS:
                # We call a p-mutation if alt(pos) / reads(pos) >= p / 100.
                # Equivalently: we call a p-mutation if
                #                         100*alt(pos) >= p*reads(pos).
                # More descriptions given in pileup.py; long story short,
                # this avoids reliance on division, and should help prevent
                # subtle floating-point errors.
                lhs = 100 * alt_nt_freq
                rhs = P_THRESHOLD * cov
                is_variant = lhs >= rhs

        if is_variant:
            vcf += f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t.\n"
            num_variants_seen_so_far += 1

        # Print occasional status updates for my sanity
        if pos % UPDATE_FREQ == 0:
            pct = 100 * (pos / seqlen)
            print(
                f"{seq}: Seen {pos:,} / {seqlen:,} ({pct:.2f}%) "
                f"positions ({num_variants_seen_so_far:,} variants so far)..."
            )

with open(VCF_FILE, "w") as vf:
    vf.write(vcf)

print(f"Wrote out the VCF file to {VCF_FILE}.")
