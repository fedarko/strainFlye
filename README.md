# strainFlye

[![strainFlye CI](https://github.com/fedarko/strainFlye/actions/workflows/main.yml/badge.svg)](https://github.com/fedarko/strainFlye/actions/workflows/main.yml) [![Code Coverage](https://codecov.io/gh/fedarko/strainFlye/branch/main/graph/badge.svg)](https://codecov.io/gh/fedarko/strainFlye)

strainFlye is a pipeline for calling, analyzing, and phasing rare mutations
in metagenome-assembled genomes produced from HiFi sequencing data.

![strainFlye pipeline diagram](https://github.com/fedarko/strainFlye/raw/main/docs/strainflye-pipeline.png)

## Note about this source code

This repository is under active development; we're working on porting our code
from a set of ad hoc analyses (https://github.com/fedarko/sheepgut) to an
easier-to-use pipeline. In the meantime, if you have any questions,
please feel free to open an issue.

## Installation

Long story short, this is a normal Python package. However, it dependends on a
few external non-Python tools (e.g. minimap2, Prodigal, SAMtools) which should
be installed via Conda. The steps below should work for most systems, but
please let me know if you encounter any problems.

```bash
wget https://raw.githubusercontent.com/fedarko/strainFlye/main/environment.yml
conda env create -f environment.yml
conda activate strainflye
pip install git+https://github.com/fedarko/strainFlye.git
```

## Command-line usage

If you installed strainFlye using the steps shown above (creating a conda
environment and installing strainFlye and its dependencies into this
environment), then you will need to activate this environment in order to use
strainFlye (e.g. `conda activate strainflye`).

<!-- STARTDOCS -->
```
Usage: strainFlye [OPTIONS] COMMAND [ARGS]...

  Pipeline for the analysis of rare mutations in metagenomes.

Options:
  -h, --help  Show this message and exit.

Commands:
  align       Aligns reads to contigs, then filters this alignment.
  call-naive  Performs naive mutation calling with controlled FDR.
  diversity   Computes the diversity index for MAGs.
  spots       Identifies hot- and/or cold-spots in MAGs.
  matrix      Computes mutation matrices of a MAG.
  link-graph  Constructs the link graph structure for a MAG.
  smooth      Generates smoothed haplotypes.
```
