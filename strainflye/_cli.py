# We define multiple commands in the strainFlye "group" using Click.
# This lets the user do things like "strainFlye align", "strainFlye call",
# etc. See https://click.palletsprojects.com/en/8.0.x/commands/ for details.
import os
import subprocess
import click
from . import cli_utils, align_utils, graph_utils, call_utils
from .errors import ParameterError


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
    """Pipeline for the analysis of rare mutations in metagenomes."""
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
    required=False,
    type=click.Path(exists=True),
    show_default="no graph",
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
        "written. Some temporary files may also be written to this directory."
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

    # Make the output dir if it doesn't already exist
    os.makedirs(output_dir, exist_ok=True)
    first_output_bam = os.path.join(output_dir, "sorted-unfiltered.bam")

    # There isn't really a need to store the SAM file from minimap2, or the
    # unsorted BAM file from "samtools view". So we use piping.
    # SAMtools stuff based on:
    # https://davetang.org/wiki/tiki-index.php?page=SAMTools#Converting_SAM_directly_to_a_sorted_BAM_file # noqa
    # Python stuff based on https://stackoverflow.com/a/4846923 and
    # https://stackoverflow.com/a/9655939

    # reads will be a tuple
    if type(reads) != tuple:
        # if we just have a single string (e.g. click changes something in the
        # future about how variadic arguments work) then we can fix this, but
        # for now let's be defensive. (If you encounter this error, go yell at
        # marcus)
        raise ParameterError("Collection of input reads should be a tuple.")

    # There's probably a way to print stuff after each individual command in
    # the chain finishes, but I don't think that sorta granularity is super
    # necessary right now tbh

    threesteps = "minimap2 --> samtools view --> samtools sort"
    fancylog(f"Running {threesteps}...")

    # NOTE: the -ax asm20 preset is what we use in the paper, but later
    # versions of minimap2 have added in "-ax map-hifi" which is probs a better
    # option in most cases. Shouldn't make too much of a difference; for
    # simplicity's sake we just stick with asm20 here, but we could definitely
    # change this (or add the option to configure it) if desired
    minimap2_run = subprocess.Popen(
        [
            "minimap2",
            "-ax",
            "asm20",
            "--secondary=no",
            "--MD",
            contigs,
            *reads,
        ],
        stdout=subprocess.PIPE,
    )
    sam_to_bam_run = subprocess.Popen(
        ["samtools", "view", "-b", "-"],
        stdin=minimap2_run.stdout,
        stdout=subprocess.PIPE,
    )
    minimap2_run.stdout.close()
    bam_to_sorted_bam_run = subprocess.Popen(
        ["samtools", "sort", "-", "-o", first_output_bam],
        stdin=sam_to_bam_run.stdout,
    )
    sam_to_bam_run.stdout.close()
    bam_to_sorted_bam_run.communicate()

    fancylog(f"Done running {threesteps}.", prefix="")

    align_utils.index_bam(first_output_bam, "sorted BAM", fancylog)

    osa_filter_bam = os.path.join(output_dir, "sorted-osa-filtered.bam")
    align_utils.filter_osa_reads(
        first_output_bam, osa_filter_bam, fancylog, verbose
    )
    align_utils.index_bam(osa_filter_bam, "OSA-filtered BAM", fancylog)

    # Now that we've finished the OSA filter step, we can remove its input BAM
    # -- the first one we generated -- to save space. These files are big
    # enough (e.g. upwards of 70 GB for the SheepGut dataset) that keeping all
    # three around at the same time might exceed the space limit on a user's
    # system.
    os.remove(first_output_bam)
    os.remove(first_output_bam + ".bai")

    pm_filter_bam = os.path.join(output_dir, "final.bam")
    align_utils.filter_pm_reads(
        graph, osa_filter_bam, pm_filter_bam, fancylog, verbose
    )
    align_utils.index_bam(pm_filter_bam, "final BAM", fancylog)

    # Similarly, we can remove the OSA-filtered (but not PM-filtered) BAM now.
    # The PM-filtered BAM represents the "final" BAM produced by the alignment
    # step, and is the one that should be used in downstream analyses.
    os.remove(osa_filter_bam)
    os.remove(osa_filter_bam + ".bai")

    fancylog("Done.")


@click.group(name="call", **grp_params, **cmd_params)
def call():
    """Na\u00efve mutation calling and diversity index computation.

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


# Use of add_command() based on https://stackoverflow.com/a/61353240.
strainflye.add_command(call)


@call.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of contigs in which to na\u00efvely call mutations.",
)
@click.option(
    "-b",
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help="BAM file representing an alignment of reads to contigs.",
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
        "Parameter that impacts the minimum coverage needed in order "
        'to consider "counting" a position / mutation towards the diversity '
        "index. Larger values increase the minimum coverage."
    ),
)
@click.option(
    "-ov",
    "--output-vcf",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output VCF file (describing the called "
        "mutations) will be written."
    ),
)
@click.option(
    "-od",
    "--output-diversity-indices",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output tab-separated values (TSV) file "
        "describing diversity indices for the values of p given in "
        "--div-index-p-list will be written."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display extra details for each contig while running.",
)
def p_mutation(
    contigs,
    bam,
    min_p,
    min_alt_pos,
    div_index_p_list,
    min_read_number,
    output_vcf,
    output_diversity_indices,
    verbose,
):
    """Calls p-mutations and computes diversity indices.

    The primary parameter for this command is the lower bound of p, defined by
    --min-p. The VCF output will include "mutations" for all positions that
    pass this (likely very low) threshold, but this VCF should be adjusted
    using the utilities contained in the "strainFlye fdr" module.
    """
    fancylog = cli_utils.fancystart(
        "strainFlye call",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("minimum p", min_p),
            ("--min-alt-pos", min_alt_pos),
            ("--div-index-p-list", div_index_p_list),
            ("minimum read number", min_read_number),
        ),
        (
            ("VCF file", output_vcf),
            ("Diversity indices file", output_diversity_indices),
        ),
    )
    di_list = call_utils.parse_di_list(div_index_p_list, param="p")
    call_utils.run(
        contigs,
        bam,
        output_vcf,
        fancylog,
        verbose,
        min_p=min_p,
        min_alt_pos=min_alt_pos,
        div_index_p_list=di_list,
    )
    fancylog("Done with p-mutation calling and diversity index computation.")


@call.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of contigs in which to na\u00efvely call mutations.",
)
@click.option(
    "-b",
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help="BAM file representing an alignment of reads to contigs.",
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
    "-ov",
    "--output-vcf",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output VCF file (describing the called "
        "mutations) will be written."
    ),
)
@click.option(
    "-od",
    "--output-diversity-indices",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output tab-separated values (TSV) file "
        "describing diversity indices for the values of r given in "
        "--div-index-r-list will be written."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display extra details for each contig while running.",
)
def r_mutation(
    contigs,
    bam,
    min_r,
    div_index_r_list,
    output_vcf,
    output_diversity_indices,
    verbose,
):
    """Calls r-mutations and computes diversity indices.

    The primary parameter for this command is the lower bound of r, defined by
    --min-r. The VCF output will include "mutations" for all positions that
    pass this (likely very low) threshold, but this VCF should be adjusted
    using the utilities contained in the "strainFlye fdr" module.
    """
    fancylog = cli_utils.fancystart(
        "strainFlye call r-mutation",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("minimum r", min_r),
            ("--div-index-r-list", div_index_r_list),
        ),
        (
            ("VCF file", output_vcf),
            ("Diversity indices file", output_diversity_indices),
        ),
    )
    di_list = call_utils.parse_di_list(div_index_r_list, param="r")
    call_utils.run(
        contigs,
        bam,
        output_vcf,
        fancylog,
        verbose,
        min_r=min_r,
        div_index_r_list=di_list,
    )
    fancylog("Done with r-mutation calling.")


# @strainflye.command(**cmd_params)
# @click.option(
#     "-c",
#     "--contigs",
#     required=True,
#     type=click.Path(exists=True),
#     help="FASTA file of contigs to compute diversity indices for.",
# )
# @click.option(
#     "-v",
#     "--vcf",
#     required=True,
#     type=click.Path(exists=True),
#     help="VCF file describing called mutations in the contigs.",
# )
# @click.option(
#     "-mc",
#     "--min-cov",
#     required=False,
#     default=10,
#     type=click.INT,
#     help=(
#         "Minimum coverage: positions with coverage less than this will not "
#         "be included in the diversity index calculation. If you don't want "
#         "to impose a minimum coverage, you can set this to 0."
#     ),
# )
# def diversity(contigs, vcf):
#     """Computes the diversity index for MAGs.
#
#     For an arbitrary MAG, let's define C as the number of (well-covered)
#     positions at which a mutation is called, and let's define T as the total
#     number of (well-covered) positions in this MAG. The diversity index is
#     then defined as the percentage C / T.
#     """
#     print("DI")
#
#
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
#     # input: contigs, reads, vcf of mutations
#     # output: contigs / graph / etc. assembled by LJA
#     print("SMOOTH")


@click.group(name="utils", **grp_params, **cmd_params)
def utils():
    """Various utility commands provided with strainFlye."""
    pass


# Nest command groups -- so we can, for example, put our "utility" commands
# under a single group, to un-clutter the top-level strainFlye CLI.
# Use of add_command() based on https://stackoverflow.com/a/61353240.
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
    # graph_utils.gfa_to_fasta() assum
    with open(output_fasta, "w") as output_fasta_file:
        num_seqs = graph_utils.gfa_to_fasta(graph, output_fasta_file)
    fancylog(f"Done.\nOutput FASTA file contains {num_seqs:,} sequences.")
