# Utilities for strainFlye's call step.


import os
import time
import pysamstats
from . import config
from .errors import ParameterError, WeirdError
from strainflye import (
    __version__,
    misc_utils,
    cli_utils,
    bcf_utils,
)


def get_alt_pos_info(rec):
    """Returns info about the second-most-common nucleotide at a position.

    This nucleotide will usually differ from the reference nucleotide, but it
    may be the reference (i.e. at positions where the reference disagrees with
    the alignment's "consensus").

    This breaks ties arbitrarily.

    Parameters
    ----------
    rec: dict
        pysamstats record for a given position in an alignment produced
        by stat_variation().

    Returns
    -------
    (cov, alt nt freq, alt nt, ref nt freq, ref nt): (int, int, str, int, str)
        Describes the first- and second-most-common nucleotides at a position.

        The first entry in this tuple is the (mis)match coverage at this
        position. This is an integer defined as the sum of A, C, G, T
        nucleotides at this position (note that this excludes degenerate
        nucleotides like N -- we could change this in the future if that'd be
        useful, I suppose). Note that this coverage could be zero, if no reads
        are aligned to this specific position.

        The second entry is the raw frequency of the second-most-common
        nucleotide at this position: this will be an integer greater than or
        equal to 0. This is also referred to in the paper, etc. as alt(pos).

        The third entry is just the alternate nucleotide (one of A, C, G, T),
        represented as a string.

        The fourth and fifth entries are analogous to the second and third, but
        this time correspond to the first-most-common nucleotide at this
        position.

    References
    ----------
    I mean, I wrote this code, but it's copied (and modified) from
    https://github.com/fedarko/sheepgut/blob/main/notebooks/pleuk_copied_code.py
    (which is in turn copied from https://github.com/fedarko/pleuk, but that
    project is still in the oven so it's private for now).
    """
    cov = rec["A"] + rec["C"] + rec["G"] + rec["T"]

    ordered_nts = sorted("ACGT", key=rec.get)

    # The literal nucleotide used in the numerator of freq(pos): one of A, C,
    # G, T
    alt_nt = ordered_nts[-2]

    # The raw frequency (in counts) of alt_nt. An integer >= 0.
    alt_freq = rec[alt_nt]

    # Replicate this stuff for the "reference" nucleotide.
    # We'll implicitly treat the reference nucleotide as the most common (aka
    # "consensus") nucleotide at this position, even at "unreasonable"
    # positions where the nucleotide located here on the contig disagrees
    # with the consensus
    ref_nt = ordered_nts[-1]

    ref_nt_freq = rec[ref_nt]

    return (cov, alt_freq, alt_nt, ref_nt_freq, ref_nt)


def parse_di_list(di_list_str, param):
    """Parses, checks, and sorts a list of p or r parameters.

    Useful when taking in a list of these for diversity index computation.
    Probably there's a more elegant way of handling this than having the user
    just give us a string, but... whatevs.

    Parameters
    ----------
    di_list_str: str
        Ideally, a comma-separated list of valid p or r parameters.
        In practice, people are going to throw a lot of stuff at this, so we'll
        want to do some thorough error checking.

    param: str
        Either "p" or "r", indicating how exactly to check these parameters.

    Returns
    -------
    di_list: list of int
        A sorted (in ascending order) and cleaned-up list of parameters.

    Raises
    ------
    ParameterError
        - If param is not "p" or "r" (come on, don't do this to me...)
        - If any of the values are invalid (we do different checks depending
          on if param is "p" or "r")
        - If any of the values are repeated
    """
    # "trust nobody. not even yourself." -- guido van rossum probably idk
    if param != "p" and param != "r":
        raise ParameterError('param must be either "p" or "r".')
    split = di_list_str.split(",")
    di_list = []
    for s in split:
        # We proooobably don't need to do this because int() is ok with
        # surrounding whitespace, but I'm paranoid.
        # (For reference, if the user passes the list of numbers in quotes
        # like --div-index-r-list "1, 2, 3, 4, 5" then spaces are possible.
        # However, if quotes aren't used, then -- at least in the fish terminal
        # -- there won't be any spaces. We can handle either scenario.)
        stripped_s = s.strip()
        try:
            val = int(stripped_s)
        except ValueError:
            raise ParameterError(
                f"We couldn't parse \"{stripped_s}\". Doesn't seem to be an "
                "integer?"
            )
        if param == "p":
            if val <= 0 or val > 5000:
                raise ParameterError(
                    f"{val} is not in the range (0, 5000], and is thus not a "
                    "valid value of p."
                )
        else:
            if val <= 0:
                raise ParameterError(
                    f"{val} is not >= 1, and is thus not a valid value of r."
                )

        di_list.append(val)

    if len(set(di_list)) != len(di_list):
        raise ParameterError(
            "The list of diversity index threshold values isn't unique."
        )

    # sort these ints in ascending order
    return sorted(di_list)


def get_min_sufficient_coverages_p(p_vals, min_read_number):
    """Computes the "minimum sufficient coverage" for value(s) of p.

    Parameters
    ----------
    p_vals: list of int in the range (0, 5000]
        Values of p.

    min_read_number: int >= 1
        Parameter of this computation.

    Returns
    -------
    list of float
        This has the exact same dimensions as p_vals. The i-th entry in this
        list corresponds to the minimum sufficient coverage for the i-th value
        of p in p_vals.
    """
    numerator = 10000 * min_read_number
    return [numerator / p for p in p_vals]


def get_min_sufficient_coverages_r(r_vals, min_cov_factor):
    """Computes the "minimum sufficient coverage" for value(s) of r.

    Parameters
    ----------
    r_vals: list of int >= 1
        Values of r.

    min_cov_factor: float
        Parameter of this computation.

    Returns
    -------
    list of float
        This has the exact same dimensions as r_vals. The i-th entry in this
        list corresponds to the minimum sufficient coverage for the i-th value
        of r in r_vals.
    """
    return [min_cov_factor * r for r in r_vals]


def run(
    contigs,
    bam,
    output_dir,
    fancylog,
    verbose,
    min_p=None,
    min_r=None,
    min_alt_pos=None,
    div_index_p_list=None,
    div_index_r_list=None,
    min_read_number=None,
    min_cov_factor=None,
):
    """Launches the process of naive p- or r-mutation calling.

    Also computes diversity index information -- we're already looping through
    each position in each contig, so it's not a huge extra cost.

    In the user-facing code in the CLI, I try to consistently say "naive" with
    the two dots, but here I don't bother because users don't see this part of
    the codebase. Here we can, like, kick off our shoes and lie down on the
    couch and eat Cheetos, because nobody is going to read this docstring (and
    if you are reading it then I'm very surprised and also hi, how's it going?)

    Parameters
    ----------
    contigs: str
        Filepath to a FASTA file containing contigs in which to call mutations.

    bam: str
        Filepath to a (sorted and indexed) BAM file mapping reads to contigs.

    output_dir: str
        Directory to which BCF / BCF index / TSV diversity index outputs will
        be written. We'll also write out temporary files -- VCF, bgzipped VCF,
        etc. (The temporary files are the motivation for using an output
        directory rather than letting users name the output files explicitly.)

    fancylog: function
        Logging function.

    verbose: bool
        Log extra info about individual contigs.

    min_p: int >= 1 or None
        Minimum p parameter, if we're calling p-mutations.
        Should be scaled up by 100, like in the CLI -- so min_p = 1 corresponds
        to 1 / 100 = 0.01%.

    min_r: int >= 1 or None
        Minimum r parameter, if we're calling r-mutations.

    min_alt_pos: int >= 1 or None
        During p-mutation calling, the second-most-common aligned nucleotide's
        frequency must be at least this. Not used when calling r-mutations.

    div_index_p_list: list of int, or None
        List of values of p for which we'll compute the diversity index.

    div_index_r_list: list of int, or None
        List of values of r for which we'll compute the diversity index.

    min_read_number: int >= 1 or None
        Parameter used in determining "minimum sufficient coverage" when
        computing diversity indices based on p-mutations.

    min_cov_factor: float >= 1 or None
        Parameter used in determining "minimum sufficient coverage" when
        computing diversity indices based on r-mutations.

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        If something went wrong with the parameters. We don't check here
        that they fall into their respective allowed ranges (since Click
        should have already enforced that thanks to its range parameter
        settings for min_p, min_r, and min_alt_pos, and parse_di_list()
        should've taken care of that for the div index lists). However, we do
        check here that only one of (min_p, min_r) is specified.
    """
    # Sanity-check whether we're using p or r. Since (currently) p-mutation and
    # r-mutation calling are separate commands, we should never see these
    # errors in practice, but you never know.
    using_p = min_p is not None
    using_r = min_r is not None

    if using_p and using_r:
        raise ParameterError(
            "p and r can't be specified at the same time. Please choose one."
        )
    if not using_p and not using_r:
        raise ParameterError("Either p or r needs to be specified.")

    # load_fasta_and_bam() will throw an error if there are < 2 contigs in the
    # FASTA. The naive calling stuff requires the use of a decoy contig, so
    # this is an ok restriction to enforce here -- I don't think this stuff
    # would be super useful with just one contig (although if you really do
    # only have one contig and want to use this script then let me know)
    contig_name2len, bf, num_fasta_contigs = misc_utils.load_fasta_and_bam(
        contigs, bam, fancylog, min_num_contigs=2
    )

    # Create the output directory (if it doesn't already exist) and determine
    # the output filepaths within this directory.
    # I'm waffling over if it makes sense to try to name these files based on
    # the threshold type and minimum used to create them (e.g.
    # "naive-p0.15.vcf"), but that will complicate using this in a pipeline.
    # Let's leave it as is for now.
    misc_utils.make_output_dir(output_dir)
    output_diversity_indices = os.path.join(
        output_dir, "diversity-indices.tsv"
    )
    output_vcf = os.path.join(output_dir, "naive-calls.vcf")

    fancylog("Creating and writing diversity index and VCF file headers...")
    # TODO? Add some checking to make sure that using_p implies that only
    # div_index_p_list is specified, and same for using_r and div_index_r_list
    # ... not high priority tho, since the user can't cause these problems from
    # the command line

    # Set up the header for the diversity index TSV file
    di_header = "Contig\tAverageCoverage\tLength"

    # ... And the header for the mutation calls' VCF file, at the same time.
    # The filter_header is important, since we will parse its ID later on to
    # determine what the minimum p or r value was. I'm not sure if there's
    # a better way to encode arbitrary file-level metadata in VCF files -- I
    # guess we could try to parse the "source" header, but I'm not a huge fan
    # of that.
    if using_p:
        param_name = "p"
        min_str = f"--min-p = {min_p / 100:.2f}%"
        filter_header = (
            f'##FILTER=<ID=strainflye_minp_{min_p}, Description="min p '
            'threshold (scaled up by 100)">\n'
        )
        min_suff_coverages = get_min_sufficient_coverages_p(
            div_index_p_list, min_read_number=min_read_number
        )
        di_list = div_index_p_list
        for p, msc in zip(di_list, min_suff_coverages):
            di_header += f"\t{config.DI_PREF}(p={p},minSuffCov={msc})"

    else:
        param_name = "r"
        min_str = f"--min-r = {min_r:,}"
        filter_header = (
            f'##FILTER=<ID=strainflye_minr_{min_r}, Description="min r '
            'threshold">\n'
        )
        min_suff_coverages = get_min_sufficient_coverages_r(
            div_index_r_list, min_cov_factor=min_cov_factor
        )
        di_list = div_index_r_list
        for r, msc in zip(di_list, min_suff_coverages):
            di_header += f"\t{config.DI_PREF}(r={r},minSuffCov={msc})"

    # Write out DI file header
    with open(output_diversity_indices, "w") as di_file:
        di_file.write(f"{di_header}\n")

    # Write out VCF file header.
    #
    # See https://samtools.github.io/hts-specs/VCFv4.3.pdf for details on the
    # VCF file format.
    #
    # We add two nonstandard info fields, in addition to those discussed in
    # the VCF docs: MDP and AAD.
    #
    # About these fields
    # ------------------
    # - MDP is explicitly *just* (mis)matching reads -- it isn't clear if DP
    #   in the VCF docs includes indels also, and I didn't want to risk
    #   misleading anyone. So I made this up.
    #
    # - AAD is like AD, but it just applies for the alternate allele. At least
    #   now, we don't need allele counts for the reference allele as well.
    #
    # About "Number"
    # --------------
    # - MDP: The 1 indicates that we only include one version of this Number
    #   per position (because it doesn't really make sense to, for example,
    #   compute depth multiple times for a position).
    #
    # - AAD: The A indicates that we will include one version of this Number
    #   per alternate allele per position (since for now we just include at
    #   most one alt allele per position, there isn't a difference from "1",
    #   but this could change if we start identifying multiallelic mutations).
    info_header = (
        "##INFO=<ID=MDP,Number=1,Type=Integer,"
        'Description="(Mis)match read depth">\n'
        "##INFO=<ID=AAD,Number=A,Type=Integer,"
        'Description="Alternate allele read depth">\n'
    )

    # Now, create a header for the VCF file listing all contigs -- this is
    # technically optional, but 1) it's recommended practice and 2) it'll help
    # a lot with parsing this VCF downstream and 3) it's needed for BCF files
    # anyway
    contig_header = ""
    for c in contig_name2len:
        contig_header += f"##contig=<ID={c},length={contig_name2len[c]}>\n"

    call_str = f"{param_name}-mutation calling ({min_str})"
    with open(output_vcf, "w") as vcf_file:
        # Header info gleaned by reading over the VCF docs
        # (https://samtools.github.io/hts-specs/VCFv4.3.pdf) and copying how
        # LoFreq organizes their header.
        # ... (It took me an embarrassing amount of time to realize that VCF
        # had already moved on to 4.3.)

        vcf_file.write(
            "##fileformat=VCFv4.3\n"
            f"##fileDate={time.strftime('%Y%m%d')}\n"
            f'##source="strainFlye v{__version__}: {call_str}"\n'
            f"##reference={os.path.abspath(contigs)}\n"
            f"{contig_header}"
            f"{info_header}"
            f"{filter_header}"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
    fancylog("Wrote out diversity index and VCF file headers.", prefix="")
    fancylog(
        "(We'll convert the VCF file to an indexed BCF afterwards.)", prefix=""
    )

    # Finally, we can get to the meat (and the most computationally expensive
    # part) of this -- go through each position in each contig and call
    # mutations, as well as observe coverages / etc.
    fancylog(f"Running {call_str} and computing diversity indices...")
    any_muts_called_across_any_contigs = False
    for si, contig in enumerate(contig_name2len, 1):
        contig_len = contig_name2len[contig]
        if verbose:
            cli_utils.proglog(
                contig, si, num_fasta_contigs, fancylog, contig_len=contig_len
            )

        # Bail out early if this contig has no linear alignments. Speeds things
        # up a bit, and lets us be clear about this to the user.
        at_least_one_aln_to_this_contig = False
        for aln in bf.fetch(contig):
            at_least_one_aln_to_this_contig = True
            break
        if not at_least_one_aln_to_this_contig:
            # We lock this behind verbose because a fair amount of contigs in,
            # for example, SheepGut are completely uncovered -- this can happen
            # if contigs are shorter than read length (weird alignment
            # artifact)
            if verbose:
                fancylog(
                    (
                        f"No linear alignments to contig {contig} exist in "
                        "the BAM file. Not performing mutation calling within "
                        "this contig."
                    ),
                    prefix="",
                )

            # Write out a line in the div idx file.
            di_line = f"{contig}\t0\t{contig_len}"
            for di in range(len(di_list)):
                di_line += "\tNA"
            with open(output_diversity_indices, "a") as di_file:
                di_file.write(f"{di_line}\n")
            # We can just exit the loop at this point. We've already included
            # this contig in the BCF header, so downstream analyses will see
            # that this contig is included in the BCF but just doesn't have any
            # mutations called within it. And we have added this contig to the
            # diversity index file, so downstream analyses can see that this
            # contig has completely undefined diversity indices and an average
            # coverage of 0x.
            continue

        # keep a running sum of coverages, so we can get the average coverage
        # for each contig. ideally i guess we'd use a fancy algorithm for
        # computing a running average (e.g.
        # https://stackoverflow.com/a/1934266), but that's hard and i don't
        # think this will be a problem (esp since Python supports arbitrarily
        # large numbers).
        coverage_sum = 0

        num_any_muts = 0
        vcf_text = ""

        # For each DI threshold, keep track of how many positions are
        # sufficiently covered for it
        msc_pos_ct = [0] * len(min_suff_coverages)
        # ... and keep track of how many mutations we observe in these
        # positions
        msc_mut_ct = [0] * len(min_suff_coverages)

        # pad=True forces us to observe uncovered positions, also
        # Note that we treat positions here as 1-indexed; this matches how
        # VCF/BCF files handle positions
        for pos, rec in enumerate(
            pysamstats.stat_variation(
                bf,
                chrom=contig,
                fafile=contigs,
                pad=True,
                max_depth=config.MAX_DEPTH_PYSAM,
            ),
            1,
        ):
            rpos = rec["pos"] + 1
            if rpos != pos:
                raise WeirdError(
                    f"Found discontinuity in traversal: {pos:,}-th pos, but "
                    f"rec['pos'] + 1 is {rpos:,}."
                )

            cov, alt_freq, alt_nt, ref_freq, ref_nt = get_alt_pos_info(rec)

            # Call a mutation naively
            if using_p:
                is_mut = call_p_mutation(alt_freq, cov, min_p, min_alt_pos)
            else:
                is_mut = call_r_mutation(alt_freq, min_r)

            if is_mut:
                any_muts_called_across_any_contigs = True
                info = get_pos_info_str(alt_freq, cov)
                # 1. CHROM = contig (aka contig name)
                #
                # 2. POS = pos (1-indexed position on this contig)
                #
                # 3. ID = . (this can be a unique identifier for this variant,
                #            but... we don't have that sort of information)
                #
                # 4. REF = ref_nt (consensus nucleotide -- this can disagree
                #                  with the actual nucleotide at this position
                #                  on the contig in the case of "unreasonable"
                #                  positions)
                #
                # 5. ALT = alt_nt (second-most-common nucleotide at this
                #                  position, for which we are calling a
                #                  mutation; note that we currently just call
                #                  at most one mutation per position, although
                #                  this could in theory be generalized to work
                #                  with multiallelic mutations)
                #
                # 6. QUAL = . (In lieu of providing a single probability here,
                #              we use the FDR estimation stuff described in the
                #              paper)
                #
                # 7. FILTER = . (We only have one "filter" right now: given
                #               the minimum p or r value, is this position a
                #               p- or r-mutation for this minimum value? If so,
                #               this implies that this position is also a
                #               p- or r-mutation for all higher values of p
                #               or r. So no need to encode discrete "filter"
                #               info for higher values of p or r.)
                #
                # 8. INFO = info (Information about coverage, alt(pos), etc.:
                #                useful when plotting FDR curves and doing
                #                other downstream analyses)
                vcf_text += (
                    f"{contig}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t{info}\n"
                )
                num_any_muts += 1

            # Compute diversity index info, also
            coverage_sum += cov
            # You might say "can't we break from the loops below early if
            # we notice early on that a position isn't a p- or r-mutation for
            # some value of p or r, since we've already sorted the threshold
            # values?" but we still wanna keep track of the positions seen with
            # some minimum sufficient coverage, so we've gotta keep going
            # through. (Although it's still probs possible to speed this up.)
            if using_p:
                for pi, p in enumerate(div_index_p_list):
                    if cov >= min_suff_coverages[pi]:
                        msc_pos_ct[pi] += 1
                        if call_p_mutation(alt_freq, cov, p, min_alt_pos):
                            msc_mut_ct[pi] += 1
            else:
                for ri, r in enumerate(div_index_r_list):
                    if cov >= min_suff_coverages[ri]:
                        msc_pos_ct[ri] += 1
                        if call_r_mutation(alt_freq, r):
                            msc_mut_ct[ri] += 1

        # Now that we've examined all positions in this contig...
        if pos != contig_len:
            raise WeirdError(
                f"For contig {contig}, the final position = {pos:,}, but the "
                f"contig length = {contig_len:,}. Something went really wrong."
            )

        # Output VCF info, if we called any mutations
        if num_any_muts > 0:
            with open(output_vcf, "a") as vcf_file:
                vcf_file.write(vcf_text)

        # Output diversity index info
        avg_cov = coverage_sum / contig_len
        di_line = f"{contig}\t{avg_cov}\t{contig_len}"

        # For each threshold value, did we observe enough sufficiently-covered
        # positions in order to compute the diversity index for this threshold
        # for this contig? If so, compute and report this diversity index.
        num_defined_di = 0
        half_contig_len = contig_len / 2
        for di in range(len(di_list)):
            if msc_pos_ct[di] >= half_contig_len:
                # Note that we don't apply any formatting to this ratio -- we
                # just export it as is. If desired we could limit the precision
                # to 6 (?) digits or something, but I don't see a strong need
                # for that right now.
                di_line += f"\t{msc_mut_ct[di] / msc_pos_ct[di]}"
                num_defined_di += 1
            else:
                di_line += "\tNA"

        with open(output_diversity_indices, "a") as di_file:
            di_file.write(f"{di_line}\n")

        if verbose:
            fancylog(
                (
                    f"Called {num_any_muts:,} {param_name}-mutation(s) "
                    f"(using {min_str}) in contig {contig}."
                ),
                prefix="",
            )
            fancylog(
                (
                    f"{num_defined_di:,} / {len(di_list):,} diversity "
                    f"indices were defined for contig {contig}."
                ),
                prefix="",
            )
    if not any_muts_called_across_any_contigs:
        raise ParameterError(
            "Could not call any mutations across any of the "
            f"{num_fasta_contigs:,} contigs. Try adjusting the minimum "
            "threshold value, or other parameters?"
        )
    fancylog(
        f"Done running {call_str} and computing diversity indices.", prefix=""
    )
    output_bcf = os.path.join(output_dir, "naive-calls.bcf")
    bcf_utils.compress_vcf(output_vcf, output_bcf, fancylog)


def get_pos_info_str(alt_pos, cov):
    info_str = f"MDP={cov};AAD={alt_pos}"
    return info_str


def call_p_mutation(alt_pos, cov, p, min_alt_pos):
    """Calls a p-mutation at a position.

    Parameters
    ----------
    alt_pos: int >= 0

    cov: int >= 0

    p: int >= 1 in (0, 50,000]
        (Still scaled up by 100, so p = 1 corresponds to 1 / 100 = 0.01%.)

    min_alt_pos: int >= 1

    Returns
    -------
    bool
    """
    # This implicitly avoids the messed-up case where alt(pos) and cov are both
    # zero, since min_alt_pos must be at least 1.
    if alt_pos >= min_alt_pos:
        # We call a p-mutation if alt(pos) / coverage(pos) >= p / 10,000.
        # (The 10,000 is because p is scaled up by 100, then by 100 again,
        # relative to (0, 50]. It's now on the order of (0, 50,000].)
        #
        # Equivalently: we call a p-mutation if
        # 10,000 * alt(pos) >= p * coverage(pos).
        #
        # This way, we avoid division and mucking around with floats.
        # The numbers might be large ints, but they're not THAT large.
        lhs = 10000 * alt_pos
        rhs = p * cov
        if lhs >= rhs:
            return True
    return False


def call_r_mutation(alt_pos, r):
    """Calls a r-mutation at a position.

    Parameters
    ----------
    alt_pos: int >= 0

    r: int >= 1

    Returns
    -------
    bool
    """
    return alt_pos >= r
