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

Long story short, strainFlye is an ordinary Python package. However, it
dependends on a few external non-Python tools (e.g. minimap2, Prodigal,
SAMtools) which it is probably easiest to install via
[conda](https://conda.io).

We offer two sets of instructions for installing strainFlye, depending on
if you want to just use the tool (Scenario 1) or if you would like to also
develop its code (Scenario 2). For most users, Scenario 1 is probably the
easiest option. (You can always install strainFlye from source later, if you'd
like!)

Both sets of installation steps below should work for most Linux or macOS
systems, but please let me know if you encounter any problems.

### Scenario 1: I just want to use strainFlye on my data!

```bash
# Download the YAML file describing the conda packages we'll install
wget https://raw.githubusercontent.com/fedarko/strainFlye/main/environment.yml

# Create a new conda environment based on this YAML file
# (by default, it'll be named "strainflye")
conda env create -f environment.yml

# Activate this conda environment
conda activate strainflye

# Install the actual strainFlye package
pip install git+https://github.com/fedarko/strainFlye.git
```

### Scenario 2: I'm interested in developing strainFlye's source code

```bash
# Clone strainFlye's source code
# (It's probably best to fork the repository first, and then clone your fork)
git clone https://github.com/fedarko/strainFlye.git

# Create and activate a conda environment, like in Scenario 1
cd strainFlye
conda env create -f environment.yml
conda activate strainflye

# Install strainFlye from the cloned source in "editable mode"
pip install -e .[dev]
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
