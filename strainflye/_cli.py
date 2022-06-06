# We define multiple commands in the strainFlye "group" using Click.
# This lets the user do things like "strainFlye align", "strainFlye call",
# etc. See https://click.palletsprojects.com/en/8.0.x/commands/ for details.
import click
from . import cli_utils, align_utils, graph_utils, call_utils, fdr_utils
from . import param_descriptions as desc


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
    "-o",
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "Directory to which an output BAM file and BAM index file will be "
        "written. Some temporary files will also be written to this directory."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display extra details for each contig during alignment filtering.",
)
# Regarding the \b marker, see https://stackoverflow.com/a/53302580 -- this is
# apparently needed to get the formatting to look the way I want (otherwise all
# of the four steps are smooshed into a paragraph)
def align(reads, contigs, graph, output_dir, verbose):
    """Aligns reads to contigs, then filters this alignment.

    Files of reads should be in the FASTA or FASTQ formats; GZIP'd files
    are allowed.

    This command involves multiple steps, including:

    \b
      1) Align reads to contigs (using minimap2) to generate a SAM file
      2) Convert this SAM file to a sorted and indexed BAM file
      3) Filter overlapping supplementary alignments within this BAM file
      4) Filter partially-mapped reads within this BAM file

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
        "strainFlye align",
        (
            ("file(s) of reads", reads_info),
            ("contig file", contigs),
            ("graph file", graph),
        ),
        (("directory", output_dir),),
    )
    align_utils.run(reads, contigs, graph, output_dir, fancylog, verbose)
    fancylog("Done.")


@click.group(name="call", **grp_params, **cmd_params)
def call():
    """[+] Na\u00efve mutation calling and diversity index computation.

    Consider a position "pos" in a contig. A given read with a (mis)match
    operation at "pos" must have one of four nucleotides (A, C, G, T) aligned
    to pos. We represent these nucleotides' counts at pos as follows:

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
    Define freq(pos) = N2 / (N1 + N2 + N3 + N4). This value, constrained to
    the range [0%, 50%], is an estimate of the mutation frequency of this
    position. We classify pos as a p-mutation if freq(pos) \u2265 p, AND if
    N2 \u2265 the --min-alt-pos parameter (an integer representing a minimum
    number of reads that must support the alternate nucleotide).

    r-mutations (na\u00efve read-count-based mutation calling)
    -----------------------------------------------------

    This takes as input some integer r > 0. We classify pos as an
    r-mutation if N2 \u2265 r.

    \b
    Diversity indices
    -----------------

    Later on in the pipeline, we'll need to select a decoy contig in order
    to perform FDR estimation for our called mutations. Contigs with low
    diversity indices may indicate promising decoy contigs; so, for the sake of
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
    type=click.IntRange(min=1),
    help=(
        "Parameter that impacts the minimum (mis)match coverage needed in "
        'order to consider "counting" a position / mutation towards the '
        "diversity index. Given a value of p (converted to the range "
        "[0, 0.5)), a position must have a coverage of at least "
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
    help=desc.VERBOSE_CALL,
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
    """Calls p-mutations and computes diversity indices.

    The primary parameter for this command is the lower bound of p, defined by
    --min-p. The BCF output will include "mutations" for all positions that
    pass this (likely very low) threshold, but this BCF should be adjusted
    using the utilities contained in the "strainFlye fdr" module.
    """
    fancylog = cli_utils.fancystart(
        "strainFlye call p-mutation",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("minimum p", min_p),
            ("--min-alt-pos", min_alt_pos),
            ("--div-index-p-list", div_index_p_list),
            ("minimum read number", min_read_number),
        ),
        (("directory", output_dir),),
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
    default="5,10,20,50,100,250,500",
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
    help=desc.VERBOSE_CALL,
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
    """Calls r-mutations and computes diversity indices.

    The primary parameter for this command is the lower bound of r, defined by
    --min-r. The BCF output will include "mutations" for all positions that
    pass this (likely very low) threshold, but this BCF should be adjusted
    using the utilities contained in the "strainFlye fdr" module.
    """
    fancylog = cli_utils.fancystart(
        "strainFlye call r-mutation",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("minimum r", min_r),
            ("--div-index-r-list", div_index_r_list),
            ("minimum coverage factor", min_coverage_factor),
        ),
        (("directory", output_dir),),
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
    """[+] FDR estimation and fixing for contigs' mutation calls."""


strainflye.add_command(fdr)


@fdr.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of contigs.",
)
@click.option(
    "-b",
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Indexed BCF file describing na\u00efvely called p- or r-mutations in "
        "the FASTA file's contigs."
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
    "--decoy-context",
    required=False,
    default="CP2",
    show_default=True,
    type=click.Choice(
        ["Full", "CP2", "Nonsyn", "Nonsense", "CP2Nonsyn", "CP2Nonsense"],
        case_sensitive=False,
    ),
    help=(
        '"Context-dependent" types of mutations to which the '
        "computation of mutation rates in the decoy contig will be limited. "
        'The "Full" option will consider the entire decoy contig, and all '
        "other options will limit the computation to certain types of "
        "positions and/or potential mutations."
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
        "Minimum average coverage of a potential decoy contig. Only used if "
        "--diversity-indices is specified."
    ),
)
@click.option(
    "-of",
    "--output-fdr-info",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output tab-separated values (TSV) file "
        "describing estimated FDRs (the ratio of the decoy contig mutation "
        "rate to the target contig mutation rate, multiplied by 100) will be "
        "written. Rows correspond to target contigs, and columns correspond "
        "to p or r values."
    ),
)
@click.option(
    "-on",
    "--output-num-info",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output tab-separated values (TSV) file "
        "describing the number of naively called mutations per megabase will "
        "be written. Has the same dimensions as the output FDR info file."
    ),
)
def estimate(
    contigs,
    bcf,
    diversity_indices,
    decoy_contig,
    decoy_context,
    high_p,
    high_r,
    decoy_min_length,
    decoy_min_average_coverage,
    output_fdr_info,
    output_num_info,
):
    """Estimates contigs' mutation calls' FDRs.

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
    """
    fancylog = cli_utils.fancystart(
        "strainFlye fdr estimate",
        (
            ("contig file", contigs),
            ("BCF file", bcf),
            ("diversity indices file", diversity_indices),
            ("manually-set decoy contig", decoy_contig),
            ("decoy contig context-dependent mutation type", decoy_context),
            (
                (
                    "high p threshold (only used if the BCF describes "
                    "p-mutations)"
                ),
                high_p,
            ),
            (
                (
                    "high r threshold (only used if the BCF describes "
                    "r-mutations)"
                ),
                high_r,
            ),
            (
                (
                    "min length of a potential decoy contig (only used if "
                    "diversity indices are specified)"
                ),
                decoy_min_length,
            ),
            (
                (
                    "min average coverage of a potential decoy contig (only "
                    "used if diversity indices are specified)"
                ),
                decoy_min_average_coverage,
            ),
        ),
        (
            ("FDR estimate file", output_fdr_info),
            ("number of mutations per megabase file", output_num_info),
        ),
    )
    fdr_utils.run_estimate(
        contigs,
        bcf,
        diversity_indices,
        decoy_contig,
        decoy_context,
        high_p,
        high_r,
        decoy_min_length,
        decoy_min_average_coverage,
        output_fdr_info,
        output_num_info,
        fancylog,
    )
    fancylog("Done.")


@fdr.command(**cmd_params)
@click.option(
    "-v",
    "--bcf",
    required=True,
    type=click.Path(exists=True),
    help="Indexed BCF file describing na\u00efvely called p- or r-mutations.",
)
@click.option(
    "-fi",
    "--fdr-info",
    required=True,
    type=click.Path(dir_okay=False),
    help='Estimated FDR TSV file produced by "strainFlye fdr estimate".',
)
@click.option(
    "--fdr",
    required=False,
    default=100,
    show_default=True,
    type=click.IntRange(min=0, min_open=True),
    help=(
        "False discovery rate at which identified mutations will be fixed. "
        "This is interpreted as a scaled-up percentage, such that f = N "
        "corresponds to (N / 100)% (i.e. the default of f = 100 corresponds "
        "to an FDR of 1%). No upper limit is imposed, since estimated FDRs "
        "can technically exceed (10000 / 100)% = 100% in uncommon cases."
    ),
)
@click.option(
    "-o",
    "--output-bcf",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output BCF file (describing the called "
        "mutations at the optimal p or r value for each contig) will be "
        "written. These mutations will be a subset of those described in the "
        "input BCF file."
    ),
)
def fix(bcf, fdr_info, fdr, output_bcf):
    """Fixes contigs' mutation calls' FDRs to an upper limit."""
    fancylog = cli_utils.fancystart(
        "strainFlye fdr fix",
        (
            ("BCF file", bcf),
            ("FDR estimate file", fdr_info),
            ("FDR to fix mutation calls at", fdr),
        ),
        (("BCF file with fixed FDR", output_bcf),),
    )
    fdr_utils.run_fix(bcf, fdr_info, fdr, output_bcf, fancylog)
    fancylog("Done.")


# @strainflye.command(**cmd_params)
# def spots():
#     """Identifies hot- and/or cold-spots in MAGs."""
#     print("H")
#
#
# @strainflye.command(**cmd_params)
# def covskew():
#     """Visualizes MAG coverage and GC skew."""
#     # input: FASTA of contigs, BAM file mapping reads to contigs
#     # output: cov / skew plots; PTR estimates, if requested?
#     print("SMOOTH")
#
#
# @strainflye.command(**cmd_params)
# def matrix():
#     """Computes mutation matrices of a MAG."""
#     print("MM")
#
#
# @strainflye.command(**cmd_params)
# def link_graph():
#     """Constructs the link graph structure for a MAG."""
#     print("LG")
#
#
# @strainflye.command(**cmd_params)
# def smooth():
#     """Generates smoothed haplotypes."""
#     # input: contigs, reads, bcf of mutations
#     # output: contigs / graph / etc. assembled by LJA
#     print("SMOOTH")


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
    """Converts a GFA 1 file to a FASTA file.

    Many of strainFlye's analyses accept a FASTA file as input (rather
    than a graph), and this is a small command that performs this conversion.
    """
    fancylog = cli_utils.fancystart(
        "strainFlye utils gfa-to-fasta",
        (("GFA file", graph),),
        (("FASTA file", output_fasta),),
    )
    with open(output_fasta, "w") as output_fasta_file:
        num_seqs = graph_utils.gfa_to_fasta(graph, output_fasta_file)
    fancylog(f"Done.\nOutput FASTA file contains {num_seqs:,} sequences.")
