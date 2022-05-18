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


@click.group(name="utils", **grp_params, **cmd_params)
def utils():
    """Various utility commands provided with strainFlye."""
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

    pm_filter_bam = os.path.join(output_dir, "final.bam")
    align_utils.filter_pm_reads(
        graph, osa_filter_bam, pm_filter_bam, fancylog, verbose
    )
    align_utils.index_bam(pm_filter_bam, "final BAM", fancylog)

    fancylog("Done.")


@strainflye.command(**cmd_params)
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
    "-p",
    required=False,
    default=None,
    type=click.FloatRange(min=0, max=50, min_open=True),
    help=(
        "Main parameter used in p-mutation (percentage-based) calling. "
        "If this is specified, r cannot be specified."
    ),
)
@click.option(
    "-r",
    required=False,
    default=None,
    type=click.IntRange(min=0, min_open=True),
    help=(
        "Main parameter used in r-mutation (read-count-based) calling. "
        "If this is specified, p cannot be specified."
    ),
)
@click.option(
    "--min-alt-pos",
    default=2,
    required=False,
    show_default=True,
    type=click.IntRange(min=0),
    help=(
        "Additional parameter for p-mutation calling: the alternate "
        "nucleotide for a p-mutation must be supported by at least this many "
        "reads. If you want to completely disable this check, you can set "
        "this to zero (but this is not recommended)."
    ),
)
@click.option(
    "-o",
    "--output-vcf",
    required=True,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which an output VCF file (describing the called "
        "mutations) will be written."
    ),
)
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Display extra details for each contig.",
)
def call(contigs, bam, p, r, min_alt_pos, output_vcf, verbose):
    """Performs na\u00efve mutation calling.

    Consider a position "pos" in a contig. A given read with a (mis)match
    operation at "pos" must have one of four nucleotides (A, C, G, T) aligned
    to pos. We represent these nucleotides' counts at pos as follows:

    \b
        N1 = # reads of the most-common aligned nucleotide at pos,
        N2 = # reads of the second-most-common aligned nucleotide at pos,
        N3 = # reads of the third-most-common aligned nucleotide at pos,
        N4 = # reads of the fourth-most-common aligned nucleotide at pos.

    (We break ties arbitrarily.)

    This command supports two types of na\u00efve mutation calling based on
    these counts: p-mutations and r-mutations. These are described below.
    You can specify either p or r to trigger p- or r-mutation calling,
    respectively; however, you can't specify both at once here.

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
    """
    fancylog = cli_utils.fancystart(
        "strainFlye call",
        (
            ("contig file", contigs),
            ("BAM file", bam),
            ("p-mutation parameter p", p),
            ("r-mutation parameter r", r),
            ("--min-alt-pos", min_alt_pos),
        ),
        (("VCF file", output_vcf),),
    )
    call_str = call_utils.run(
        contigs, bam, output_vcf, min_alt_pos, fancylog, verbose, p=p, r=r
    )
    fancylog(f"Done with {call_str}.")


@strainflye.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of contigs for which diversity indices will be computed.",
)
@click.option(
    "-v",
    "--vcf",
    required=True,
    type=click.Path(exists=True),
    help="VCF file describing called mutations in the contigs.",
)
@click.option(
    "-mc",
    "--min-cov",
    required=False,
    default=10,
    type=click.INT,
    help=(
        "Minimum coverage: positions with coverage less than this will not be "
        "included in the diversity index calculation. If you don't want to "
        "impose a minimum coverage, you can set this to 0."
    ),
)
def diversity(contigs, vcf):
    """Computes the diversity index for MAGs.

    For an arbitrary MAG, let's define C as the number of (well-covered)
    positions at which a mutation is called, and let's define T as the total
    number of (well-covered) positions in this MAG. The diversity index is then
    defined as the percentage C / T.
    """
    print("DI")


@strainflye.command(**cmd_params)
def spots():
    """Identifies hot- and/or cold-spots in MAGs."""
    print("H")


@strainflye.command(**cmd_params)
def covskew():
    """Visualizes MAG coverage and GC skew."""
    # input: FASTA of contigs, BAM file mapping reads to contigs
    # output: cov / skew plots; PTR estimates, if requested?
    print("SMOOTH")


@strainflye.command(**cmd_params)
def matrix():
    """Computes mutation matrices of a MAG."""
    print("MM")


@strainflye.command(**cmd_params)
def link_graph():
    """Constructs the link graph structure for a MAG."""
    print("LG")


@strainflye.command(**cmd_params)
def smooth():
    """Generates smoothed haplotypes."""
    # input: contigs, reads, vcf of mutations
    # output: contigs / graph / etc. assembled by LJA
    print("SMOOTH")


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
