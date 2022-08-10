# strainFlye

<a href="https://github.com/fedarko/strainFlye/actions/workflows/main.yml"><img src="https://github.com/fedarko/strainFlye/actions/workflows/main.yml/badge.svg" alt="strainFlye CI" /></a>
<a href="https://codecov.io/gh/fedarko/strainFlye"><img src="https://codecov.io/gh/fedarko/strainFlye/branch/main/graph/badge.svg" alt="Code Coverage" /></a>

strainFlye is a pipeline for calling, analyzing, and phasing rare mutations
in metagenome-assembled genomes produced from HiFi sequencing data.

<img src="https://github.com/fedarko/strainFlye/raw/main/docs/strainflye-pipeline.png" alt="strainFlye pipeline diagram" />

## Note about this source code

This repository is under active development; we're working on porting our code
from a set of ad hoc analyses (https://github.com/fedarko/sheepgut) to an
easier-to-use pipeline. In the meantime, if you have any questions,
please feel free to open an issue.

## Installation

Long story short, strainFlye is an ordinary Python package. However, it
dependends on a few external non-Python tools (e.g. minimap2, Prodigal,
SAMtools) which it is probably easiest to install via
[conda](https://conda.io).

Assuming that you already have conda installed, the following installation
instructions should work for most Linux or macOS systems. (However, please feel
free to open an issue if you encounter any problems.)

```bash
# Download the YAML file describing the conda packages we'll install
# (if you don't have wget, you could also just download this file manually)
wget https://raw.githubusercontent.com/fedarko/strainFlye/main/environment.yml

# Create a new conda environment based on this YAML file
# (by default, it'll be named "strainflye")
conda env create -f environment.yml

# Activate this conda environment
conda activate strainflye

# Install the actual strainFlye package
pip install git+https://github.com/fedarko/strainFlye.git
```

## Documentation

### Tutorial

We provide a (work in progress) tutorial using the SheepGut dataset
**[here](https://nbviewer.org/github/fedarko/strainFlye/blob/main/docs/SheepGutExample.ipynb)**.

### Command-line usage

If you installed strainFlye using the steps shown above (creating a conda
environment and installing strainFlye and its dependencies into this
environment), then you will need to activate this environment in order to use
strainFlye (e.g. `conda activate strainflye`).

<!-- STARTDOCS -->
```
Usage: strainFlye [OPTIONS] COMMAND [ARGS]...

  Pipeline for the analysis of rare mutations in metagenomes.

  Please consult https://github.com/fedarko/strainFlye if you have any
  questions, comments, etc. about strainFlye. Thank you for using this tool!

Options:
  -h, --help  Show this message and exit.

Commands:
  align  Align reads to contigs, and filter the resulting alignment.
  call   [+] Call mutations in contigs naïvely & compute diversity indices.
  fdr    [+] Estimate and fix FDRs for contigs' naïve mutation calls.
  spot   [+] Identify hotspots or coldspots in contigs.
  utils  [+] Various utility commands provided with strainFlye.
```

## Acknowledgements

### Test data
`sample1.gfa` (located in `strainflye/tests/inputs/`)
was downloaded from the [gfalint](https://github.com/sjackman/gfalint)
repository. The other GFA files in this folder beginning with `sample1` are
also based on this GFA file.
