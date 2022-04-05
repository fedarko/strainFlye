# We define multiple commands in the strainFlye "group" using Click.
# This lets the user do things like "strainFlye align", "strainFlye call",
# etc. See https://click.palletsprojects.com/en/8.0.x/commands/ for details.
import os
import time
import subprocess
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
    "--output-dir",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False),
    help="Directory to which an output BAM file and index will be written.",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Show extra details while running.",
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
    """
    t0 = time.time()

    def fancylog(msg):
        t1 = time.time()
        print(f"--------\n{t1 - t0:.2f} sec: {msg}")

    fancylog("Starting strainFlye align...")
    if verbose:
        # Print list of reads files
        reads_info = ""
        for i, rf in enumerate(reads, 1):
            reads_info += f"\n{i}. {rf}"
        fancylog(f"Input file(s) of reads:{reads_info}")
        fancylog(f"Input contig file: {contigs}")
        fancylog(f"Output directory: {output_dir}")

    # Make the output dir if it doesn't already exist
    os.makedirs(output_dir, exist_ok=True)
    output_bam = os.path.join(output_dir, "sorted-unfiltered.bam")

    # There isn't really a need to store the SAM file from minimap2, or the
    # unsorted BAM file from "samtools view". So we use piping.
    # SAMtools stuff based on:
    # https://davetang.org/wiki/tiki-index.php?page=SAMTools#Converting_SAM_directly_to_a_sorted_BAM_file # noqa
    # Python stuff based on https://stackoverflow.com/a/4846923 and
    # https://stackoverflow.com/a/9655939

    # reads will be a tuple
    if type(reads) == tuple:
        reads_str = " ".join(reads)
    else:
        # if we just have a single string (e.g. click changes something in the
        # future about how variadic arguments work) then we can fix this, but
        # for now let's be defensive. (If you encounter this error, go yell at
        # marcus)
        raise ValueError("The reads should be a tuple, right?")

    if verbose:
        fancylog(
            f"FYI, we're telling minimap2 that the reads are:\n{reads_str}"
        )

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
            reads_str,
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
        ["samtools", "sort", "-", "-o", output_bam],
        stdin=sam_to_bam_run.stdout,
    )
    sam_to_bam_run.stdout.close()
    bam_to_sorted_bam_run.communicate()

    if verbose:
        fancylog("Done running {threesteps}.")
        fancylog("Indexing this BAM...")

    subprocess.run(["samtools", "index", output_bam])

    if verbose:
        fancylog("Indexed the BAM.")

    # TODO! Invoke the filters, re-indexing after each
    fancylog("Filtering overlapping supplementary alignments...")
    fancylog("Filtering partially-mapped reads...")

    fancylog("Done running strainFlye align.")


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
def smooth():
    """Generates smoothed haplotypes."""
    # input: contigs, reads, vcf of mutations
    # output: contigs / graph / etc. assembled by LJA
    print("SMOOTH")
