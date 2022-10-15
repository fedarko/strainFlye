# We define multiple commands in the strainFlye "group" using Click.
# This lets the user do things like "strainFlye align", "strainFlye call",
# etc. See https://click.palletsprojects.com/en/8.0.x/commands/ for details.
import click
from strainflye import (
    cli_utils,
    align_utils,
    graph_utils,
    call_utils,
    fdr_utils,
    spot_utils,
    smooth_utils,
    link_utils,
    dynam_utils,
    config,
    __version__,
)
from strainflye import param_descriptions as desc


# By default, Click's help info (shown when running e.g. "strainFlye -h")
# lists strainFlye's commands in alphabetical order; this isn't ideal for
# pipelines like ours, where there is a reasonably natural ordering of the
# commands. We can get around this, and sort commands in the help info in
# the order in which we define them within this file, by adjusting how
# the strainFlye Click.Group works; see
# https://github.com/pallets/click/issues/513#issuecomment-504158316.
class ClickGroupWithOrderedCommands(click.Group):
    def list_commands(self, ctx):
        return self.commands.keys()


grp_params = {"cls": ClickGroupWithOrderedCommands}

cmd_params = {
    # github.com/marbl/MetagenomeScope/blob/master/metagenomescope/_cli.py:
    # make -h work as an alternative to --help
    "context_settings": {"help_option_names": ["-h", "--help"]},
    # https://click.palletsprojects.com/en/8.0.x/api/: if the user just
    # types "strainflye" or "strainflye align", show them the help for this
    # command rather than yelling at them for not providing required options
    "no_args_is_help": True,
}


@click.group(**grp_params, **cmd_params)
@click.version_option(__version__, "-v", "--version")
def strainflye():
    """Pipeline for the analysis of rare mutations in metagenomes.

    Please consult https://github.com/fedarko/strainFlye
    if you have any questions, comments, etc. about strainFlye.
    Thank you for using this tool!
    """
    pass


@strainflye.command(**cmd_params)
# We use arguments to allow for multiple read files -- this works well with
# globbing in unix, i.e. "strainflye align read*.fasta [other stuff here...]"
# see https://stackoverflow.com/a/34763795
@click.argument("reads", required=True, type=click.Path(exists=True), nargs=-1)
# I think a GZIP'd file is ok for minimap2 here, also, but I'm not 100% sure
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of contigs to which reads will be aligned.",
)
@click.option(
    "-g",
    "--graph",
    default=None,
    show_default="no graph",
    required=False,
    type=click.Path(exists=True),
    help=(
        "GFA 1-formatted file describing an assembly graph of the contigs. "
        'This is used in the "filter partially-mapped reads" step to make the '
        "filter less strict for reads mapped to adjacent contigs in the "
        "graph. This isn't required; if it isn't passed, adjacent contigs "
        "will not be considered in this filter."
    ),
)
@click.option(
    "-p",
    "--minimap2-params",
    required=False,
    default=config.DEFAULT_MM2_PARAMS,
    show_default=True,
    type=click.STRING,
    help=(
        "Additional parameters to pass to minimap2, besides the contig and "
        "read information. Depending on the size of your dataset and the "
        "amount of memory your system has, you may want to adjust the -I "
        "parameter; see minimap2's documentation for details. Please note "
        "that we do not perform any validation on this string before passing "
        "it to minimap2 (so if you are allowing users to run strainFlye "
        "through a web server, be careful about shell injection)."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which an output BAM file (final.bam) and BAM index file "
        "will be written. Some temporary files will also be written to this "
        "directory."
    ),
)
@click.option(
    "--rm-tmp-bam/--no-rm-tmp-bam",
    is_flag=True,
    default=True,
    show_default=True,
    help=(
        "Remove temporary (before filtering, and after just filtering OSAs) "
        "BAM files, to reduce space requirements."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display extra details during alignment filtering.",
)
# Regarding the \b marker, see https://stackoverflow.com/a/53302580 -- this is
# apparently needed to get the formatting to look the way I want (otherwise all
# of the four steps are smooshed into a paragraph)
def align(
    reads, contigs, graph, minimap2_params, output_dir, rm_tmp_bam, verbose
):
    """Align reads to contigs, and filter the resulting alignment.

    Files of reads should be in the FASTA or FASTQ formats; GZIP'd files
    are allowed.

    This command involves multiple steps, including:

    \b
      1) Align reads to contigs (using minimap2) to generate a SAM file
      2) Convert this SAM file to a sorted and indexed BAM file
      3) Filter overlapping supplementary alignments (OSAs) from this BAM file
      4) Filter partially-mapped reads from this BAM file

    Note that we only sort the alignment file once, although we do re-index it
    after the two filtering steps. This decision is motivated by
    https://www.biostars.org/p/131333/#131335.
    """

    # Convert collection of reads files into something more easy to "read"
    # (I'm here all night, folks)
    reads_info = ""
    for i, rf in enumerate(reads, 1):
        reads_info += f"\n{i}. {rf}"

    # Get a snazzy logging function we can use
    fancylog = cli_utils.fancystart(
        "align",
        (
            ("file(s) of reads", reads_info),
            ("contig file", contigs),
            ("graph file", graph),
        ),
        (("directory", output_dir),),
        extra_info=(
            f"minimap2 parameters: {minimap2_params}",
            f"Remove temporary BAM files?: {cli_utils.b2y(rm_tmp_bam)}",
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    align_utils.run(
        reads,
        contigs,
        graph,
        minimap2_params,
        output_dir,
        rm_tmp_bam,
        fancylog,
        verbose,
    )
    fancylog("Done.")


@click.group(name="call", **grp_params, **cmd_params)
def call():
    """[+] Call mutations in contigs na\u00efvely & compute diversity indices.

    Consider a position "pos" in a contig. Using the alignment, we can count
    how many reads have a (mis)match operation to "pos" with one of the four
    nucleotides (A, C, G, T; we ignore degenerate nucleotides in reads).
    We represent these four nucleotides' counts at pos as follows:

    \b
        N1 = # reads of the most-common aligned nucleotide at pos,
        N2 = # reads of the second-most-common aligned nucleotide at pos,
        N3 = # reads of the third-most-common aligned nucleotide at pos,
        N4 = # reads of the fourth-most-common aligned nucleotide at pos.

    (We break ties arbitrarily.)

    strainFlye supports two types of na\u00efve mutation calling based on
    these counts: p-mutations and r-mutations. These are described below.

    p-mutations (na\u00efve percentage-based mutation calling)
    -----------------------------------------------------

    This takes as input some percentage p in the range (0%, 50%].
    Define freq(pos) = N2 / (N1 + N2 + N3 + N4). This value, which is
    inherently constrained to the range [0%, 50%], is an estimate of the
    mutation frequency of this position. We classify pos as a p-mutation
    if freq(pos) \u2265 p, AND if N2 \u2265 the --min-alt-pos parameter
    (an integer representing a minimum number of reads that must support
    the alternate nucleotide).

    r-mutations (na\u00efve read-count-based mutation calling)
    -----------------------------------------------------

    This takes as input some integer r > 0. We classify pos as an
    r-mutation if N2 \u2265 r.

    \b
    Diversity indices
    -----------------

    Later on in the pipeline, if we perform FDR estimation on these called
    mutations using the "strainFlye fdr" module, we will need to select
    a "decoy contig" with relatively few mutations. A contig with low diversity
    indices may represent a promising decoy contig; so, for the sake of
    convenience, both commands output this information.
    """
    pass


# Nest command groups -- so we can, for example, put our "utility" commands
# under a single group, to un-clutter the top-level strainFlye CLI.
# Use of add_command() based on https://stackoverflow.com/a/61353240.
strainflye.add_command(call)


@call.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_CONTIGS_NAIVE_CALL,
)
@click.option(
    "-b",
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BAM,
)
@click.option(
    "--min-p",
    required=False,
    show_default=True,
    default=50,
    type=click.IntRange(min=0, max=5000, min_open=True),
    help=(
        "Minimum value of p for which to call p-mutations. This is scaled "
        "up by 100 (i.e. the default of 50 corresponds to 50 / 100 = 0.5%) "
        "in order to bypass floating-point precision issues."
    ),
)
@click.option(
    "--min-alt-pos",
    default=2,
    required=False,
    show_default=True,
    type=click.IntRange(min=1),
    help=(
        "In order for us to call a p-mutation at a position, this position's "
        "alternate nucleotide must be supported by at least this many reads."
    ),
)
@click.option(
    "--div-index-p-list",
    default="50,100,200,500,1000,2500,5000",
    required=False,
    show_default=True,
    type=click.STRING,
    help=(
        "List of values of p for which we'll compute diversity indices. "
        "These should all be separated by commas; and, as with --min-p, "
        "these are scaled up by 100. Please don't use commas as thousands "
        "separators."
    ),
)
@click.option(
    "-m",
    "--min-read-number",
    default=5,
    required=False,
    show_default=True,
    type=click.FloatRange(min=0, min_open=True),
    help=(
        "Parameter that impacts the minimum (mis)match coverage needed in "
        'order to consider "counting" a position / mutation towards the '
        "diversity index. Given a value of p (converted to the range "
        "(0, 0.5]), a position must have a coverage of at least "
        ' (--min-read-number / p) in order to be "sufficiently covered" and '
        "thus counted towards the diversity index."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=desc.OUTPUT_DIR_NAIVE_CALL,
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def p_mutation(
    contigs,
    bam,
    min_p,
    min_alt_pos,
    div_index_p_list,
    min_read_number,
    output_dir,
    verbose,
):
    """Call p-mutations and compute diversity indices.

    The primary parameter for this command is the lower bound of p, defined by
    --min-p. The BCF output will include "mutations" for all positions that
    pass this (likely very low) threshold; this BCF can be filtered
    using the utilities contained in the "strainFlye fdr" module.
    """
    fancylog = cli_utils.fancystart(
        "call p-mutation",
        (
            ("contig file", contigs),
            ("BAM file", bam),
        ),
        (("directory", output_dir),),
        extra_info=(
            f"Minimum p: {min_p:,}",
            (
                "Minimum number of supporting reads needed to call a "
                f"p-mutation: {min_alt_pos:,}"
            ),
            (
                "Values of p for which we'll compute diversity indices: "
                f"{div_index_p_list}"
            ),
            (
                "Minimum read number (used in diversity index computation): "
                f"{min_read_number:,}"
            ),
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    di_list = call_utils.parse_di_list(div_index_p_list, param="p")
    call_utils.run(
        contigs,
        bam,
        output_dir,
        fancylog,
        verbose,
        min_p=min_p,
        min_alt_pos=min_alt_pos,
        div_index_p_list=di_list,
        min_read_number=min_read_number,
    )
    fancylog("Done.")


@call.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_CONTIGS_NAIVE_CALL,
)
@click.option(
    "-b",
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BAM,
)
@click.option(
    "--min-r",
    required=False,
    default=5,
    show_default=True,
    type=click.IntRange(min=1),
    help="Minimum value of r for which to call r-mutations.",
)
@click.option(
    "--div-index-r-list",
    default="2,3,4,5,10,20,50,100,250,500",
    required=False,
    show_default=True,
    type=click.STRING,
    help=(
        "List of values of r for which we'll compute diversity indices. "
        "These should all be separated by commas. Please don't use commas "
        "as thousands separators."
    ),
)
@click.option(
    "-f",
    "--min-coverage-factor",
    default=2,
    required=False,
    show_default=True,
    type=click.FloatRange(min=1),
    help=(
        "Parameter that impacts the minimum (mis)match coverage needed in "
        'order to consider "counting" a position / mutation towards the '
        "diversity index. Given a value of r, a position must have a "
        "coverage of at least (--min-coverage-factor \u00d7 r) in order to be "
        '"sufficiently covered" and thus counted towards the diversity index.'
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=desc.OUTPUT_DIR_NAIVE_CALL,
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def r_mutation(
    contigs,
    bam,
    min_r,
    div_index_r_list,
    min_coverage_factor,
    output_dir,
    verbose,
):
    """Call r-mutations and compute diversity indices.

    The primary parameter for this command is the lower bound of r, defined by
    --min-r. The BCF output will include "mutations" for all positions that
    pass this (likely very low) threshold; this BCF can be filtered
    using the utilities contained in the "strainFlye fdr" module.

    (However, as discussed in the Supplemental Material of the strainFlye
    paper, we note that FDR estimation using r-mutations can be problematic
    due to coverage variations between contigs. So, if you plan to make use of
    the FDR estimation / fixing utilities in strainFlye, we suggest using
    p-mutations instead.)
    """
    fancylog = cli_utils.fancystart(
        "call r-mutation",
        (
            ("contig file", contigs),
            ("BAM file", bam),
        ),
        (("directory", output_dir),),
        extra_info=(
            f"Minimum r: {min_r:,}",
            (
                "Values of r for which we'll compute diversity indices: "
                f"{div_index_r_list}"
            ),
            (
                "Minimum coverage factor (used in diversity index "
                f"computation): {min_coverage_factor:,.2f}"
            ),
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    di_list = call_utils.parse_di_list(div_index_r_list, param="r")
    call_utils.run(
        contigs,
        bam,
        output_dir,
        fancylog,
        verbose,
        min_r=min_r,
        div_index_r_list=di_list,
        min_cov_factor=min_coverage_factor,
    )
    fancylog("Done.")


@click.group(name="fdr", **grp_params, **cmd_params)
def fdr():
    """[+] Estimate and fix FDRs for contigs' na\u00efve mutation calls."""


strainflye.add_command(fdr)


@fdr.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "FASTA file of contigs for which we will estimate mutation calls' "
        "FDRs. All contigs in this FASTA file should also be "
        "contained in the BAM and BCF files. It's ok if the BAM file "
        "contains contigs not in this FASTA file (we'll ignore them), "
        "but all contigs described in the BCF file's header must also be in "
        "this FASTA file."
    ),
)
@click.option(
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BAM,
)
@click.option(
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Indexed BCF file describing na\u00efvely called p- or r-mutations in "
        'the FASTA file\'s contigs, produced by "strainFlye call".'
    ),
)
@click.option(
    "-di",
    "--diversity-indices",
    default=None,
    show_default="nothing",
    required=False,
    type=click.Path(exists=True),
    help=(
        "TSV file describing the diversity indices of a set of contigs. "
        "Used to automatically select a decoy contig. "
        "This option is mutually exclusive with --decoy-contig."
    ),
)
@click.option(
    "-dc",
    "--decoy-contig",
    default=None,
    show_default="nothing",
    required=False,
    type=click.STRING,
    help=(
        "Name of a specific contig to use as the decoy contig for FDR "
        "estimation. This option is mutually exclusive with "
        "--diversity-indices."
    ),
)
@click.option(
    "-dctx",
    "--decoy-contexts",
    required=True,
    default=["CP2"],
    show_default=True,
    multiple=True,
    type=click.Choice(config.DCTX_ALL, case_sensitive=False),
    help=(
        '"Context-dependent" types of positions and/or mutations to which '
        "the computation of mutation rates in the decoy contig will be "
        'limited. The "Full" option will consider the entire decoy contig; '
        "all other options will limit the computation to certain types of "
        "positions and/or potential mutations. (CP2 will focus on positions "
        "in the second codon position of genes predicted in the decoy contig "
        "using Prodigal; Nonsyn will focus on potential nonsynonymous "
        "mutations in these genes; Nonsense will focus on potential nonsense "
        "mutations in these genes; and Tv will focus on potential "
        "transversion mutations.) You can specify this option multiple times "
        "to generate multiple sets of FDR estimates. (You can also specify "
        f'"{config.DCTX_EVERYTHING}" to use all available contexts.)'
    ),
)
@click.option(
    "-hp",
    "--high-p",
    required=False,
    default=500,
    show_default=True,
    type=click.IntRange(min=0, max=5000, min_open=True),
    help=(
        "(Only applies if the input BCF file describes p-mutations.) "
        "p-mutations with a mutation rate "
        "(freq(pos)) greater than or equal to this are considered "
        '"indisputable," and will not '
        "be included in FDR estimation. Like the values of p used as "
        'parameters to "strainFlye call p-mutation", this is in the range '
        "(0, 5000], such that h = N corresponds to (N / 100)%. Corresponds "
        'to the "highFrequency" threshold mentioned in the paper.'
    ),
)
@click.option(
    "-hr",
    "--high-r",
    required=False,
    default=100,
    show_default=True,
    type=click.IntRange(min=0, min_open=True),
    help=(
        "(Only applies if the input BCF file describes r-mutations.) "
        "r-mutations with an alternate nucleotide read "
        "coverage greater than or equal to this are considered "
        "indisputable, and will not be included in FDR estimation."
    ),
)
@click.option(
    "-dml",
    "--decoy-min-length",
    required=False,
    default=1000000,
    show_default=True,
    type=click.IntRange(min=0, min_open=True),
    help=(
        "Minimum length of a potential decoy contig. Only used if "
        "--diversity-indices is specified."
    ),
)
@click.option(
    "-dmac",
    "--decoy-min-average-coverage",
    required=False,
    default=500,
    show_default=True,
    type=click.FloatRange(min=0, min_open=True),
    help=(
        "Minimum average (mis)match coverage of a potential decoy contig. "
        "Only used if --diversity-indices is specified."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which TSV files describing the estimated FDRs and "
        "numbers of na\u00efvely called mutations per megabase will be "
        "written. In all TSV files, rows correspond to target contigs and "
        "columns correspond to p or r values. We will generate one estimated "
        "FDR TSV file for every time --decoy-context is specified, and just "
        "one number-of-na\u00efvely-called-mutations-per-megabase file."
    ),
)
def estimate(
    contigs,
    bam,
    bcf,
    diversity_indices,
    decoy_contig,
    decoy_contexts,
    high_p,
    high_r,
    decoy_min_length,
    decoy_min_average_coverage,
    output_dir,
):
    """Estimate the FDRs of contigs' mutation calls.

    We do this using the target-decoy approach (TDA). Given a set of C contigs,
    we select a "decoy contig" with relatively few called mutations. We then
    compute a mutation rate for this decoy contig, and use this mutation rate
    (along with the mutation rates of the other C - 1 "target" contigs) to
    estimate the FDRs of all of these target contigs' mutation calls.

    We can produce multiple FDR estimates for a single target contig's calls by
    varying the p or r threshold used (from the --min-p or --min-r threshold
    used to generate the input BCF file, up to the --high-p or --high-r
    threshold given here). Using this information (and information about the
    numbers of mutations called per megabase), we can plot an FDR curve for
    a given target contig's mutation calls.

    This command accepts an input BCF file of p- or r-mutations; however, in
    general we recommend using p-mutations (rather than r-mutations) for FDR
    estimation / fixing. Please see the Supplemental Material section of the
    strainFlye paper named "Identifying mutations based solely on read
    counts" for details.
    """
    fancylog = cli_utils.fancystart(
        "fdr estimate",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("BCF file", bcf),
            ("diversity indices file", diversity_indices),
        ),
        (("directory", output_dir),),
        extra_info=(
            f"Manually-set decoy contig: {decoy_contig}",
            (
                "Decoy contig context-dependent position / mutation type(s): "
                f"{decoy_contexts}"
            ),
            (
                "High p threshold (only used if the BCF describes "
                f"p-mutations): {high_p:,}"
            ),
            (
                "High r threshold (only used if the BCF describes "
                f"r-mutations): {high_r:,}"
            ),
            (
                "Minimum length of a potential decoy contig (only used if "
                f"diversity indices are specified): {decoy_min_length:,} bp"
            ),
            (
                "Minimum average coverage of a potential decoy contig (only "
                "used if diversity indices are specified): "
                f"{decoy_min_average_coverage:,}x"
            ),
        ),
    )
    fdr_utils.run_estimate(
        contigs,
        bam,
        bcf,
        diversity_indices,
        decoy_contig,
        decoy_contexts,
        high_p,
        high_r,
        decoy_min_length,
        decoy_min_average_coverage,
        output_dir,
        fancylog,
    )
    fancylog("Done.")


@fdr.command(**cmd_params)
@click.option(
    "-b",
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Indexed BCF file describing na\u00efvely called p- or r-mutations, "
        'produced by "strainFlye call".'
    ),
)
@click.option(
    "-fi",
    "--fdr-info",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "One of the estimated FDR TSV files produced by "
        '"strainFlye fdr estimate".'
    ),
)
@click.option(
    "--fdr",
    required=False,
    default=1,
    show_default=True,
    type=click.FloatRange(min=0, min_open=True),
    help=(
        "False discovery rate at which identified mutations will be fixed. "
        "This is interpreted as a percentage: the default of 1 corresponds "
        "to an FDR of 1%. We do not restrict this to an upper limit, since "
        "estimated FDRs can technically exceed 100%."
    ),
)
@click.option(
    "-o",
    "--output-bcf",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output indexed BCF file (describing the called "
        "mutations at the optimal p or r value for each contig, along with "
        'any "indisputable" mutations relative to the --high-p or --high-r '
        "parameter used for strainFlye fdr estimate) will be written."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display extra details while writing the filtered BCF.",
)
def fix(bcf, fdr_info, fdr, output_bcf, verbose):
    """Fix contigs' mutation calls' estimated FDRs to an upper limit.

    This takes as input the estimated FDRs from "strainFlye fdr estimate" (if
    you used multiple decoy contexts, then you will need to choose which set of
    FDR estimates to use here) to guide us on how to fix the FDR for each
    contig. Note that mutations that passed the "high" p or r threshold
    specified for "strainFlye fdr estimate", and thus were not used for FDR
    estimation, will all be included in the output BCF file from this command;
    these mutations are considered "indisputable."

    We include indisputable mutations from the decoy contig and from all
    target contigs in our output BCF file. We will only consider including
    non-indisputable mutations from the target contigs: the decision of
    which non-indisputable mutations will be included is based on the lowest
    p or r parameter for a target contig that yields an estimated FDR \u2264
    the fixed FDR given here.
    """
    fancylog = cli_utils.fancystart(
        "fdr fix",
        (
            ("BCF file", bcf),
            ("FDR estimate file", fdr_info),
        ),
        (("BCF file with mutation calls at the fixed FDR", output_bcf),),
        extra_info=(
            f"FDR at which to fix mutation calls: {fdr:.2f}%",
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    fdr_utils.run_fix(bcf, fdr_info, fdr, output_bcf, fancylog, verbose)
    fancylog("Done.")


@click.group(name="spot", **grp_params, **cmd_params)
def spot():
    """[+] Identify putative mutational hotspots or coldspots in contigs.

    Many methods exist for identifying these sorts of hotspots or coldspots;
    strainFlye's implementations of these methods are intended mostly as a
    quick proof-of-concept for replicating the results shown in our
    paper, and are not extremely "feature-rich" quite yet.
    """


strainflye.add_command(spot)


@spot.command(**cmd_params)
@click.option(
    "-b",
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BCF_DOWNSTREAM,
)
@click.option(
    "-f",
    "--features",
    required=True,
    type=click.Path(exists=True),
    help=(
        'Generic Feature Format version 3 (GFF3) file describing "features" '
        "(e.g. predicted protein-coding genes) in the contigs. We will ignore "
        "features located on contigs that are not described in the BCF file's "
        "header."
    ),
)
@click.option(
    "-mn",
    "--min-num-mutations",
    required=False,
    type=click.IntRange(min=0, min_open=True),
    default=None,
    show_default="no check",
    help=(
        'Label a feature as a "hotspot" if it contains at least this many '
        "mutations."
    ),
)
@click.option(
    "-mp",
    "--min-perc-mutations",
    required=False,
    type=click.FloatRange(min=0, max=100, min_open=True),
    default=None,
    show_default="no check",
    help=(
        'Label a feature as a "hotspot" if its percentage of mutations ((# '
        "mutations / feature length) \u00d7 100) is at least "
        "this value."
    ),
)
@click.option(
    "-o",
    "--output-hotspots",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output tab-separated values (TSV) file "
        "describing hotspot features across all contigs will be written."
    ),
)
def hot_features(
    bcf,
    features,
    min_num_mutations,
    min_perc_mutations,
    output_hotspots,
):
    """Identify hotspot features (for example, genes).

    By "feature", we refer to a single continuous region within a contig,
    as described in the file given for --features. These regions could describe
    anything: predicted protein-coding genes, introns or exons, intergenic
    regions of interest, etc. For now, we treat each feature independently
    (e.g. we don't lump together exons from the same "Parent" gene; each
    feature is considered separately as a potential "hotspot").

    You can configure whether or not we classify a feature as a hotspot by
    adjusting the --min-num-mutations and --min-perc-mutations parameters; at
    least one of these parameters must be specified. If both parameters are
    specified, then both checks (number of mutations in a feature, and
    percentage of mutations in a feature) will need to pass in order for us to
    label a feature as a hotspot.
    """
    fancylog = cli_utils.fancystart(
        "spot hot-features",
        (
            ("BCF file", bcf),
            ("feature file", features),
        ),
        (("file describing hotspot features", output_hotspots),),
        extra_info=(
            (
                "Minimum number of mutations needed to call a feature a "
                f"hotspot: {min_num_mutations:,}"
            ),
            (
                "Minimum % of mutations needed to call a feature a hotspot: "
                f"{min_perc_mutations:.2f}%"
            ),
        ),
    )
    spot_utils.run_hotspot_feature_detection(
        bcf,
        features,
        min_num_mutations,
        min_perc_mutations,
        output_hotspots,
        fancylog,
    )
    fancylog("Done.")


@spot.command(**cmd_params)
@click.option(
    "-b",
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BCF_DOWNSTREAM,
)
@click.option(
    "-l",
    "--min-length",
    required=True,
    type=click.IntRange(min=0, min_open=True),
    default=5000,
    show_default=True,
    help=(
        'Label a gap between mutations in a contig as a "coldspot" if the gap '
        "is at least this long."
    ),
)
@click.option(
    "--circular/--no-circular",
    is_flag=True,
    default=False,
    show_default=True,
    required=False,
    help=(
        "If --circular is specified, we'll assume that all contigs are "
        'circular: we\'ll consider the gap "looping around" from the '
        "rightmost mutation in each contig to the leftmost mutation in "
        "the contig as a potential coldspot. Otherwise, we will assume all "
        "contigs are linear."
    ),
)
@click.option(
    "--exact-pvals/--no-exact-pvals",
    is_flag=True,
    default=True,
    show_default=True,
    required=False,
    help=(
        "If --exact-pvals is specified, we'll use the method for "
        "computing exact longest-gap p-values given in equation (3.1) of "
        "(Naus 1982). We use Python's decimal library with default parameters "
        "to try to deal with large numbers, but we can't guarantee that this "
        "won't cause the program to fail in certain cases. "
        "If --exact-pvals is not specified, we'll use equation (3.3) of (Naus "
        "1982), which gives an approximation of the p-value."
    ),
)
@click.option(
    "-o",
    "--output-coldspots",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output tab-separated values (TSV) file "
        "describing coldspots will be written. The longest gap in each "
        "contig will also be given a p-value, indicating the probability "
        "of the longest gap being at least this long (given the length of "
        "and number of mutated positions in this contig, and assuming that "
        "mutations occur with a constant probability at any position in the "
        "contig). See (Naus 1982) for details."
    ),
)
def cold_gaps(bcf, min_length, circular, exact_pvals, output_coldspots):
    """Identify long coldspot "gaps" without any mutations.

    To clarify, we define a "gap" of length L on a contig as
    a collection of continuous positions [N, N + 1, ... N + L - 2, N + L - 1]
    in which no positions are mutations (based on the input BCF file).

    If the --circular flag is specified, then we can loop around the contig
    from right to left; otherwise, the left and right sides of the contig are
    hard boundaries. To give an example of this, consider a 9-nucleotide
    contig that has mutations at positions 4 and 6:

    \b
                               Mut.    Mut.
                    1   2   3   4   5   6   7   8   9

    If --circular is specified, then this contig has two gaps: one gap of
    length 1 (covering just position 5, between the two mutations), and another
    gap of length 6 (starting at position 7 and looping around to position 3:
    [7, 8, 9, 1, 2, 3]).

    If --circular is not specified, then this contig has three gaps: [1, 2, 3],
    [5], and [7, 8, 9].
    """
    fancylog = cli_utils.fancystart(
        "spot cold-gaps",
        (("BCF file", bcf),),
        (("file describing coldspot gaps", output_coldspots),),
        extra_info=(
            f"Minimum coldspot gap length: {min_length:,} bp",
            f"Check for circular coldspot gaps?: {cli_utils.b2y(circular)}",
            (
                "Compute exact longest-gap p-values?: "
                f"{cli_utils.b2y(exact_pvals)}"
            ),
        ),
    )
    spot_utils.run_coldspot_gap_detection(
        bcf, min_length, circular, exact_pvals, output_coldspots, fancylog
    )
    fancylog("Done.")


@click.group(name="smooth", **grp_params, **cmd_params)
def smooth():
    """[+] Create and assemble smoothed and virtual reads."""


strainflye.add_command(smooth)


@smooth.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "FASTA file of contigs for which we will create smoothed and virtual "
        f"reads. {desc.FASTA_SUBSET_BAM_BCF}"
    ),
)
@click.option(
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BAM,
)
@click.option(
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BCF_DOWNSTREAM,
)
@click.option(
    "-di",
    "--diversity-indices",
    default=None,
    show_default="nothing",
    required=False,
    type=click.Path(exists=True),
    help=(
        "TSV file describing the diversity indices of a set of contigs, "
        'produced by one of the "strainFlye call" commands. '
        "Only used if --virtual-reads is specified. Along with "
        "diversity indices, these files list contigs' average coverages. "
        "If --virtual-reads is specified, we will need to know contigs' "
        "average coverages. So, if a diversity indices file is provided here, "
        "then we can avoid re-computing average coverages. "
        "(Otherwise, if --virtual-reads is specified but no diversity indices "
        "file is provided, we'll need to compute average coverages; this will "
        "take some extra time.)"
    ),
)
@click.option(
    "--virtual-reads/--no-virtual-reads",
    is_flag=True,
    default=True,
    show_default=True,
    required=False,
    help=(
        "If --virtual-reads is specified, we'll create virtual reads "
        'covering "low-coverage" regions in contigs.'
    ),
)
@click.option(
    "-vrp",
    "--virtual-read-well-covered-perc",
    required=True,
    type=click.FloatRange(min=0, max=100),
    default=95,
    show_default=True,
    help=(
        "Only used if --virtual-reads is specified. In a contig with average "
        "coverage (based on the BAM file, and only considering match + "
        "mismatch counts) C, we will define a position in this contig (with "
        "coverage P, based only on the smoothed reads) as low-coverage "
        "if ((P / C) \u00d7 100) is less than this percentage."
    ),
)
@click.option(
    "-vrf",
    "--virtual-read-flank",
    required=True,
    type=click.IntRange(min=0),
    default=100,
    show_default=True,
    help=(
        "Only used if --virtual-reads is specified. When we add virtual "
        "reads spanning a single continuous low-coverage region, these reads "
        "will start and end this many positions before and after the region. "
        "(For example, the default of 100 means that the virtual reads "
        "created for a low-coverage region of 5,000 bp will all be "
        "5,200 bp.)"
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which smoothed reads (and virtual reads, if "
        "--virtual-reads is specified) will be written. Each contig's reads "
        "will be written to a gzipped FASTA file in this directory named "
        "[contig].fasta.gz."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def create(
    contigs,
    bam,
    bcf,
    diversity_indices,
    virtual_reads,
    virtual_read_well_covered_perc,
    virtual_read_flank,
    output_dir,
    verbose,
):
    """Create smoothed and virtual reads for each contig."""
    fancylog = cli_utils.fancystart(
        "smooth create",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("BCF file", bcf),
            ("diversity indices file", diversity_indices),
        ),
        (("directory", output_dir),),
        extra_info=(
            f"Create virtual reads?: {cli_utils.b2y(virtual_reads)}",
            (
                '"Well-covered" percentage used for virtual read creation: '
                f"{virtual_read_well_covered_perc:.2f}%"
            ),
            (
                'Create "flanks" of this length on either side of each '
                f"virtual read: {virtual_read_flank:,}"
            ),
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    smooth_utils.run_create(
        contigs,
        bam,
        bcf,
        diversity_indices,
        virtual_reads,
        virtual_read_well_covered_perc,
        virtual_read_flank,
        output_dir,
        verbose,
        fancylog,
    )
    fancylog("Done.")


@smooth.command(**cmd_params)
@click.option(
    "-r",
    "--reads-dir",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help=(
        'Directory produced by "strainFlye smooth create" containing smoothed '
        "(and optionally virtual) reads. We will use LJA to assemble each "
        "*.fasta.gz file in this directory (representing reads from different "
        "contigs) independently."
    ),
)
@click.option(
    "-p",
    "--lja-params",
    required=False,
    default=config.DEFAULT_LJA_PARAMS,
    show_default=True,
    type=click.STRING,
    help=(
        "Additional parameters to pass to LJA, besides the --reads and "
        "--output-dir parameters. To explain our defaults: the --simpleec "
        "flag is currently only available on the simple_ec branch of "
        'LJA. This flag tells LJA to perform "simple" error correction by '
        "removing all edges in the de Bruijn graph with k-mer coverage less "
        "than --Cov-threshold (10), as well as reads passing through these "
        'edges. This "simple" error correction is done instead of running '
        "the default LJA error correction module, mowerDBG. Please note that "
        "we do not perform any validation on this string before passing it to "
        "LJA (so if you are allowing users to run strainFlye through a web "
        "server, be careful about shell injection)."
    ),
)
@click.option(
    "-b",
    "--lja-bin",
    required=False,
    default=None,
    show_default="Look for LJA's bin in the $PATH environment variable",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, executable=True
    ),
    help=(
        'Location of LJA\'s "lja" binary file, which can be used to run LJA. '
        "This should be located in the bin/ directory constructed when you "
        "compiled LJA. If this is not provided, we will check to see if LJA "
        "is available in your $PATH."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which output LJA assemblies (one top-level directory "
        " per *.fasta.gz file in the input --reads-dir) will be written."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def assemble(reads_dir, lja_params, lja_bin, output_dir, verbose):
    """Assemble contigs' smoothed and virtual reads using LJA.

    Please note that this command relies on the "simple_ec" branch of LJA being
    installed on your system. See strainFlye's README (and/or LJA's manual) for
    details on installing LJA.
    """
    fancylog = cli_utils.fancystart(
        "smooth assemble",
        (("reads directory", reads_dir),),
        (("directory", output_dir),),
        extra_info=(
            f"LJA parameters: {lja_params}",
            f"LJA binary (if None, we'll look in $PATH): {lja_bin}",
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    smooth_utils.run_assemble(
        reads_dir, lja_params, lja_bin, output_dir, verbose, fancylog
    )
    fancylog("Done.")


@click.group(name="link", **grp_params, **cmd_params)
def link():
    """[+] Create link graphs showing co-occurring alleles.

    The "nt" command should be run first; this generates information about
    the counts and co-occurrences of nucleotides at mutated positions in a
    contig.

    The "graph" command takes as input the information produced by "nt" and
    creates link graphs (one link graph per contig) from it. There are many
    parameters that impact the graph creation, so we have split this into
    two steps in order to make it easy to rerun the "graph" step with different
    parameter settings if needed (the "nt" command will likely take longer to
    run).
    """


strainflye.add_command(link)


@link.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "FASTA file of contigs for which we will create link graphs. "
        f"{desc.FASTA_SUBSET_BAM_BCF}"
    ),
)
@click.option(
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BAM,
)
@click.option(
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BCF_DOWNSTREAM,
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which information about nucleotide frequencies at the "
        "mutated positions in each contig, as well as co-occurrence "
        "information between pairs of mutated positions, will be written. "
        'Each contig will be represented by two "pickle" files within this '
        "directory, one describing nucleotide frequencies at each position "
        "and one describing nucleotide pair frequencies at pairs of "
        "positions; each file's name will be prefixed with the corresponding "
        "contig name."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def nt(contigs, bam, bcf, output_dir, verbose):
    """Compute (co-)occurrence information for mutations' nucleotides."""
    fancylog = cli_utils.fancystart(
        "link nt",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("BCF file", bcf),
        ),
        (("directory", output_dir),),
        extra_info=(f"Verbose?: {cli_utils.b2y(verbose)}",),
    )
    link_utils.run_nt(contigs, bam, bcf, output_dir, verbose, fancylog)
    fancylog("Done.")


@link.command(**cmd_params)
@click.option(
    "-n",
    "--nt-dir",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    help=(
        'Directory produced by "strainFlye link nt" containing "pickle" files '
        "describing nucleotide (co-)occurrence information for certain "
        "contigs."
    ),
)
@click.option(
    "--min-nt-count",
    required=False,
    default=2,
    show_default=True,
    type=click.IntRange(min=1),
    help=(
        "We will only create a node representing a given nucleotide at a "
        "mutated position if this nucleotide is supported by at least this "
        "many reads at this position."
    ),
)
@click.option(
    "--min-span",
    required=False,
    default=501,
    show_default=True,
    type=click.IntRange(min=1),
    help=(
        "One of the prerequisites for creating an edge between two nodes "
        "(representing nucleotides Ni, Nj at mutated positions i and j) is "
        "that the total number of reads spanning both i and j, across all "
        "combinations of (non-degenerate) nucleotides at both i and j, is at "
        "least this parameter."
    ),
)
@click.option(
    "--low-link",
    required=False,
    default=0,
    show_default=True,
    type=click.FloatRange(min=0, max=1, max_open=True),
    help=(
        "The other prerequisite for creating an edge between two nodes "
        "(representing nucleotides Ni, Nj at mutated positions i and j) is "
        "that link(i, j, Ni, Nj) > this parameter. Note this check is not "
        "inclusive (i.e. the default of 0 means that, if reads(i, j, Ni, Nj) "
        "= 0, we will never create an edge between (i, Ni) and (j, Nj)."
    ),
)
@click.option(
    "-f",
    "--output-format",
    required=False,
    type=click.Choice(["dot", "nx"], case_sensitive=False),
    default="dot",
    show_default=True,
    help=(
        'Format in which to write out each contig\'s link graph. "dot" will '
        'write out graphs in Graphviz\' DOT language; "nx" will write out '
        'graphs to "pickle" files, in which each link graph is a NetworkX '
        "object."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which graphs will be written. We'll write one graph "
        "file per contig; each file's name will be prefixed with the "
        "corresponding contig name."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def graph(
    nt_dir,
    min_nt_count,
    min_span,
    low_link,
    output_format,
    output_dir,
    verbose,
):
    """Convert (co-)occurrence information into a link graph structure.

    We create one link graph per contig. These graphs are undirected. A node
    (i, Ni) in the graph represents an allele: in other words, the occurrence
    of nucleotide Ni (one of {A, C, G, T}) at position i (1-indexed).

    What do edges between nodes represent? To answer that, let's define
    reads(i, Ni) as the number of reads at which Ni is aligned to i. Let's also
    define reads(i, j, Ni, Nj) as the number of reads at which Ni is aligned to
    i, and Nj is aligned to j.

    Given these definitions, we define the link weight between two alleles (i,
    Ni) and (j, Nj) as

    \b
                                 reads(i, j, Ni, Nj)
     link(i, j, Ni, Nj) = -----------------------------------
                            max(reads(i, Ni), reads(j, Nj))

    We create an edge between nodes (i, Ni) and (j, Nj) if both of the
    following conditions hold:

    1. The number of reads spanning both i and j with any (non-degenerate)
    nucleotides aligned to either position is \u2265 the --min-span parameter.

    2. link(i, j, Ni, Nj) > the --low-link parameter.
    """
    fancylog = cli_utils.fancystart(
        "link graph",
        (("nucleotide information directory", nt_dir),),
        (("directory", output_dir),),
        extra_info=(
            (
                "To make a node (i, Ni): reads(i, Ni) must be \u2265 "
                f"{min_nt_count:,}"
            ),
            (
                "To make an edge btwn (i, Ni) & (j, Nj): spanCount(i, j) must "
                f"be \u2265 {min_span:,}"
            ),
            (
                "To make an edge btwn (i, Ni) & (j, Nj): link(i, j, Ni, Nj) "
                f"must be > {low_link:,.2f}"
            ),
            f"Output graph file format: {output_format}",
            f"Verbose?: {cli_utils.b2y(verbose)}",
        ),
    )
    link_utils.run_graph(
        nt_dir,
        min_nt_count,
        min_span,
        low_link,
        output_format,
        output_dir,
        verbose,
        fancylog,
    )
    fancylog("Done.")


@click.group(name="dynam", **grp_params, **cmd_params)
def dynam():
    """[+] Compute simple information about growth dynamics."""


strainflye.add_command(dynam)


@dynam.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help=(
        "FASTA file of contigs for which we will compute coverage and skew "
        f"information. {desc.FASTA_SUBSET_BAM}"
    ),
)
@click.option(
    "-b",
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help=desc.INPUT_BAM,
)
@click.option(
    "--bin-length",
    required=False,
    default=10000,
    show_default=True,
    type=click.IntRange(min=1),
    help=(
        "Bin length, used for both coverage and skew. If a contig's length is "
        "larger than this length but not evenly divisible by it, the "
        'rightmost bin will be a smaller one that contains all "overflowing" '
        "positions. If a contig's length is less than or equal to this bin "
        "length, we'll just create one bin for this contig."
    ),
)
@click.option(
    "-e",
    "--norm-coverage-epsilon",
    required=False,
    default=0.3,
    show_default=True,
    type=click.FloatRange(min=0, min_open=True),
    help=(
        "We clamp all bins' normalized coverages in a contig to within this "
        "distance of 1.0. For example, the default of 0.3 means we clamp to "
        "the range [0.7, 1.3]."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which TSV files describing contigs' coverage and skew "
        "will be written. We'll write one TSV file per contig; each file "
        "will be named [contig]_covskew.tsv."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help=desc.VERBOSE_BASIC,
)
def covskew(contigs, bam, bin_length, norm_coverage_epsilon, output_dir):
    """Compare normalized coverage and skew within contigs."""
    fancylog = cli_utils.fancystart(
        "dynam covskew",
        (
            ("contig file", contigs),
            ("BAM file", bam),
        ),
        (("directory", output_dir),),
        extra_info=(
            f"Bin length: {bin_length:,} bp",
            f"Normalized coverage epsilon: {norm_coverage_epsilon}",
        ),
    )
    dynam_utils.run_covskew(
        contigs,
        bam,
        bin_length,
        norm_coverage_epsilon,
        output_dir,
        verbose,
        fancylog,
    )
    fancylog("Done.")


# @strainflye.command(**cmd_params)
# def matrix():
#     """Computes mutation matrices of a MAG."""
#     print("MM")


@click.group(name="utils", **grp_params, **cmd_params)
def utils():
    """[+] Various utility commands provided with strainFlye."""
    pass


strainflye.add_command(utils)


@utils.command(**cmd_params)
@click.option(
    "-g",
    "--graph",
    required=True,
    type=click.Path(exists=True),
    help=(
        "GFA 1-formatted file describing an assembly graph. This assumes "
        "that the assembly graph contains a sequence for each segment "
        "(i.e. there are no segment lines where the sequence is omitted "
        "and replaced with a * symbol)."
    ),
)
@click.option(
    "-o",
    "--output-fasta",
    required=True,
    type=click.Path(),
    help="Output FASTA file containing the contigs in this input graph.",
)
def gfa_to_fasta(graph, output_fasta):
    """Convert a GFA 1 file to a FASTA file.

    Many of strainFlye's analyses accept a FASTA file as input (rather
    than a graph), and this is a small command that performs this conversion.
    """
    fancylog = cli_utils.fancystart(
        "utils gfa-to-fasta",
        (("GFA file", graph),),
        (("FASTA file", output_fasta),),
    )
    with open(output_fasta, "w") as output_fasta_file:
        num_seqs = graph_utils.gfa_to_fasta(graph, output_fasta_file)
    fancylog(f"Done.\nOutput FASTA file contains {num_seqs:,} sequences.")
