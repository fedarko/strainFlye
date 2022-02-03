# We define multiple commands in the strainFlye "group" using Click.
# This lets the user do things like "strainFlye align", "strainFlye call",
# etc. See https://click.palletsprojects.com/en/8.0.x/commands/ for details.
import click


cmd_params = {
    # github.com/marbl/MetagenomeScope/blob/master/metagenomescope/_cli.py:
    # make -h work as an alternative to --help
    "context_settings": {"help_option_names": ["-h", "--help"]},
    # https://click.palletsprojects.com/en/8.0.x/api/: if the user just
    # types "strainflye" or "strainflye align", show them the help for this
    # command rather than yelling at them for not providing required options
    "no_args_is_help": True,
}


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


@click.group(cls=ClickGroupWithOrderedCommands, **cmd_params)
def strainflye():
    """Pipeline for the analysis of rare mutations in metagenomes."""
    pass


@strainflye.command(**cmd_params)
# Ideally we'd let the user specify multiple reads files without having to
# repeatedly say "-r", but this doesn't seem possible with options in Click;
# see https://stackoverflow.com/q/48391777. Oh well! This state of affairs is
# at least consistent with LJA + jumboDBG as of writing.
@click.option(
    "-r",
    "--reads",
    required=True,
    type=click.Path(exists=True),
    multiple=True,
    help=(
        "FASTA or FASTQ file(s). GZIP'd files are allowed. "
        "You can use this option multiple times if you'd like to specify "
        "multiple files of reads at once."
    ),
)
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
    required=True,
    type=click.Path(exists=True),
    help=(
        "GFA-formatted file describing an assembly graph of the contigs. "
        'This is used in the "filter partially-mapped reads" step to make the '
        "filter less strict for reads mapped to adjacent contigs in the graph."
    ),
)
@click.option(
    "-o",
    "--output-bam",
    required=True,
    type=click.Path(dir_okay=False),
    help="Filepath to which an output BAM file will be written.",
)
# Regarding the \b marker, see https://stackoverflow.com/a/53302580 -- this is
# apparently needed to get the formatting to look the way I want (otherwise all
# of the four steps are smooshed into a paragraph)
def align(reads, contigs, output_bam):
    """Aligns reads to contigs, and filters this alignment.

    This involves multiple steps, including:

    \b
      1) Align reads to contigs (using minimap2) to generate a SAM file
      2) Convert this SAM file to a sorted and indexed BAM file
      3) Filter overlapping supplementary alignments within this BAM file
      4) Filter partially-mapped reads within this BAM file
    """
    print("Aligning")
    # ... TODO actually do this!
    # I think we can merge the SAM -> BAM step together with the
    # sorting-and-indexing step, and begin this step immediately from
    # minimap2's output using piping; not sure if we'll need to re-(sort and
    # index) after performing filtering
    #
    # See https://github.com/lh3/minimap2/issues/350,
    # http://www.htslib.org/doc/samtools-sort.html
    print(f"reads are: {reads}")


@strainflye.command(**cmd_params)
@click.option(
    "-c",
    "--contigs",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of contigs in which to call mutations.",
)
@click.option(
    "-b",
    "--bam",
    required=True,
    type=click.Path(exists=True),
    help="BAM file representing an alignment of reads to contigs.",
)
@click.option(
    "-f",
    "--fdr",
    required=False,
    type=click.FloatRange(min=0),
    default=5,
    help=(
        "False Discovery Rate (FDR) to fix mutation calls at. This "
        "corresponds to a percentage, so this should usually be a value "
        "within the range [0, 100]. In theory the target-decoy approach "
        "can estimate FDRs greater than 100%, which is why we have left "
        "the upper bound here open."
    ),
)
@click.option(
    "-dc",
    "--decoy-contig",
    required=False,
    type=click.STRING,
    help=(
        "A contig name in the --contigs FASTA file. "
        "You can use this option to force strainFlye to use a specific "
        'decoy contig (e.g. "edge_6104").'
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
    "-ofc",
    "--output-fdr-curve",
    required=False,
    type=click.Path(dir_okay=False),
    help="Filepath to which an output FDR curve PNG image will be written.",
)
@click.option(
    "-ofd",
    "--output-fdr-data",
    required=False,
    type=click.Path(dir_okay=False),
    help=(
        "Filepath to which a TSV file describing the FDR curve data will be "
        "written. This can be useful if you'd prefer to handle FDR curve "
        "plotting yourself."
    ),
)
def call_naive(
    contigs,
    bam,
    fdr,
    decoy_contig,
    output_vcf,
    output_fdr_curve,
    output_fdr_data,
):
    """Performs naive mutation calling with controlled FDR.

    The FDR (false discovery rate) is estimated based on the target-decoy
    approach.
    """
    # This should contain stuff from bam-to-pileup.py. Need to:
    #
    # save mutation info about MAGs to a space-efficient format (ideally
    # better than the "pileup" format I set up for the jupyter notebooks),
    #
    # then quickly
    # (greedily?) find a decoy genome [high coverage, lowest #muts at a given
    # threshold],
    #
    # then, given the specified FDR -- use binary search or
    # something to control mutations with this FDR.
    #
    # TODO?: Either just select best decoy filter combination (e.g. CP2) or
    # set it manually in advance. likely do that and expose as parameter here
    # (with transversions, CP2, nonsyn, nonsense), maybe with
    # multiple=True so that these can be added together
    print("Calling")


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
def matrix():
    """Computes mutation matrices of a MAG."""
    print("MM")


@strainflye.command(**cmd_params)
def link_graph():
    """Constructs the link graph structure for a MAG."""
    print("LG")


@strainflye.command(**cmd_params)
def smooth_jumbo():
    """Generates smoothed haplotypes using jumboDBG."""
    print("SJ")
