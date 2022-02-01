# We define multiple commands in the strainFlye "group" using Click.
# This lets the user do things like "strainFlye align", "strainFlye call",
# etc. See https://click.palletsprojects.com/en/8.0.x/commands/ for details.
import click


@click.group()
def strainflye():
    """Pipeline for the analysis of rare mutations in metagenomes."""
    pass


@strainflye.command()
def align():
    """Performs alignment."""
    print("Aligning")


@strainflye.command()
def call_naive():
    """Performs naive mutation identification using NaiveFreq."""
    print("Calling")


@strainflye.command()
def est_fdr():
    """Estimates the FDR of identified mutations.

    This is done using the target-decoy approach."""
    print("Estimating")


@strainflye.command()
def div_idx():
    """Computes the diversity index of MAGs."""
    print("DI")


@strainflye.command()
def spots():
    """Identifies hot- and/or cold-spots in MAGs."""
    print("H")


@strainflye.command()
def mut_matrix():
    """Computes mutation matrices of a MAG."""
    print("MM")


@strainflye.command()
def link_graph():
    """Constructs the link graph structure for a MAG."""
    print("LG")


@strainflye.command()
def smooth_jumbo():
    """Generates smoothed haplotypes using jumboDBG."""
    print("SJ")
