# strainFlye

<a href="https://github.com/fedarko/strainFlye/actions/workflows/main.yml"><img src="https://github.com/fedarko/strainFlye/actions/workflows/main.yml/badge.svg" alt="strainFlye CI" /></a>
<a href="https://codecov.io/gh/fedarko/strainFlye"><img src="https://codecov.io/gh/fedarko/strainFlye/branch/main/graph/badge.svg" alt="Code Coverage" /></a>

strainFlye is a pipeline for calling, analyzing, and phasing rare mutations
in metagenome-assembled genomes produced from HiFi sequencing data. The main
inputs to strainFlye are 1) reads and 2) contigs.
However, most steps in the pipeline can be "jumped to" if you already have
other files (e.g. alignment, mutation calls) ready.

<img src="https://github.com/fedarko/strainFlye/raw/main/docs/strainflye-pipeline.png" alt="strainFlye pipeline diagram" />

## Installation

Long story short, strainFlye is an ordinary Python package. However, it
dependends on a few external non-Python tools (e.g. minimap2, Prodigal,
SAMtools); most of these can be installed using [conda](https://conda.io).

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

This should be set up as a conda package soon.

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

We provide a tutorial using the SheepGut dataset **[in this Jupyter Notebook](https://nbviewer.org/github/fedarko/strainFlye/blob/main/docs/SheepGutExample.ipynb)**.

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
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  align   Align reads to contigs, and filter the resulting alignment.
  call    [+] Call mutations in contigs naïvely & compute diversity indices.
  fdr     [+] Estimate and fix FDRs for contigs' naïve mutation calls.
  spot    [+] Identify putative mutational hotspots or coldspots.
  smooth  [+] Create and assemble smoothed and virtual reads.
  link    [+] Create link graphs showing co-occurring alleles.
  matrix  [+] Create codon and amino acid mutation matrices.
  dynam   [+] Compute simple information about growth dynamics.
  utils   [+] Miscellaneous utility commands provided with strainFlye.
```

### Development documentation

If you're interested in making changes to strainFlye's code, please see
[`CONTRIBUTING.md`](https://github.com/fedarko/strainFlye/blob/main/CONTRIBUTING.md)
for some tips on getting started.

## Miscellaneous notes about minor details you probably don't need to care about

### Shell injection (only relevant if this is hosted on a web server)

Some of strainFlye's commands use Python's
[`subprocess` module](https://docs.python.org/3/library/subprocess.html) to run
non-Python software: minimap2, samtools, bcftools, Prodigal, LJA, etc.
Most of the time, we do this using `subprocess.run()` with `shell=False`:
long story short, this helps prevent the problem of
[shell injection](https://en.wikipedia.org/wiki/Code_injection#Shell_injection).

However, as of writing, there are two places where strainFlye uses
`subprocess.run()` with `shell=True`: in `strainFlye align` (when running
minimap2 and samtools), and in `strainFlye smooth assemble` (when running LJA).
This is for convenience's sake, since we allow the user to pass in extra
parameters to these commands (the `--minimap2-params` option for `strainFlye
align`, and the `--lja-params` option for `strainFlye smooth assemble`).

Our use of `shell=True` in these two cases means that it's possible to make
these commands do unexpected things (see
[Python's documentation here](https://docs.python.org/3/library/subprocess.html#security-considerations) for details).
*This should not be a problem if you are running strainFlye directly.* However,
if you decide to host strainFlye on a server somewhere (and you allow users to
upload files, specify parameters, etc.) then you should be careful about
preventing shell injection in these cases. Feel free to open an issue if you
have any questions about this.

## Acknowledgements

### Test data
`sample1.gfa` (located in `strainflye/tests/inputs/`)
was downloaded from the [gfalint](https://github.com/sjackman/gfalint)
repository. The other GFA files in this folder beginning with `sample1` are
also based on this GFA file.
