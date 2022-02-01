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

```bash
wget https://raw.githubusercontent.com/fedarko/strainFlye/main/setup_requirements.txt
pip install -r setup_requirements.txt
pip install git+https://github.com/fedarko/strainFlye.git
```

## Command-line usage

<!-- STARTDOCS -->
```
Usage: strainFlye [OPTIONS] COMMAND [ARGS]...

  Pipeline for the analysis of rare mutations in metagenomes.

Options:
  -h, --help  Show this message and exit.

Commands:
  align         Aligns reads to contigs, and processes this alignment.
  call-naive    Performs naive mutation identification using NaiveFreq.
  est-fdr       Estimates the FDR of identified mutations.
  div-idx       Computes the diversity index of MAGs.
  spots         Identifies hot- and/or cold-spots in MAGs.
  mut-matrix    Computes mutation matrices of a MAG.
  link-graph    Constructs the link graph structure for a MAG.
  smooth-jumbo  Generates smoothed haplotypes using jumboDBG.
```
