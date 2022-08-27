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
SAMtools); almost all of these can be installed using
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

### Optional: install LJA in order to run `strainFlye smooth assemble`

strainFlye's `smooth` module includes two commands. The first,
`strainFlye smooth create`, creates smoothed and virtual reads for each contig;
the second, `strainFlye smooth assemble`, assembles these reads using
[LJA](https://github.com/AntonBankevich/LJA). LJA is not installed using the
conda installation instructions above, so—in order to run the
`strainFlye smooth assemble` command—you will need to
install the [LJA](https://github.com/AntonBankevich/LJA) software (in
particular, the
[`simple_ec` branch](https://github.com/AntonBankevich/LJA/tree/simple_ec) of
LJA's code).

Please see [LJA's manual](https://github.com/AntonBankevich/LJA/blob/main/docs/lja_manual.md)
for the most up-to-date installation instructions. Assuming that you have all
of LJA's requirements installed, something like the following should work:

```bash
git clone https://github.com/AntonBankevich/LJA.git
cd LJA
git checkout simple_ec
cmake .
make
```

... but this is subject to change as LJA is updated.

## Documentation

### Tutorial

We provide a (work in progress) tutorial using the SheepGut dataset
**[in this Jupyter Notebook](https://nbviewer.org/github/fedarko/strainFlye/blob/main/docs/SheepGutExample.ipynb)**.

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
  align   Align reads to contigs, and filter the resulting alignment.
  call    [+] Call mutations in contigs naïvely & compute diversity indices.
  fdr     [+] Estimate and fix FDRs for contigs' naïve mutation calls.
  spot    [+] Identify putative mutational hotspots or coldspots in contigs.
  smooth  [+] Create and assemble smoothed and virtual reads.
  link    [+] Create link graphs in order to show mutations' co-occurrences.
  utils   [+] Various utility commands provided with strainFlye.
```

### Development documentation

If you're interested in making changes to strainFlye's code, please see
[`CONTRIBUTING.md`](https://github.com/fedarko/strainFlye/blob/main/CONTRIBUTING.md)
for some tips on getting started.

## Acknowledgements

### Test data
`sample1.gfa` (located in `strainflye/tests/inputs/`)
was downloaded from the [gfalint](https://github.com/sjackman/gfalint)
repository. The other GFA files in this folder beginning with `sample1` are
also based on this GFA file.
