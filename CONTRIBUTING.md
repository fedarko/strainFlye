# strainFlye development documentation

Thanks for your interest in improving strainFlye's code!

## Setting up a development environment

Assuming that you already have conda installed, the following installation
instructions should work for most Linux or macOS systems. (However, please feel
free to open an issue if you encounter any problems.)

First, you should [fork strainFlye](https://docs.github.com/en/get-started/quickstart/fork-a-repo);
this will make contributing your changes back up to the main branch easier.

After this, the following instructions should work (assuming that conda is
installed):

```bash
# Clone your fork of strainFlye's source code
git clone https://github.com/your-github-username-goes-here/strainFlye.git

# Create and activate a conda environment based on the "environment.yml" file
# (by default, it'll be named "strainflye")
cd strainFlye
conda env create -f environment.yml
conda activate strainflye

# Install strainFlye from the cloned source in "editable mode"
pip install -e .[dev]
```

### Sidenote: using mamba instead of conda

You may prefer to create the environment using
[mamba](https://mamba.readthedocs.io) instead of conda, since mamba can be much
faster. As of writing, probably the easiest way to do this is to create an
empty environment and then
[update it based on `environment.yml`](https://github.com/mamba-org/mamba/issues/633#issuecomment-812272143).

### Sidenote: installing LJA in order to use `strainFlye smooth assemble`

If you want to work on `strainFlye smooth assemble` then you may also want to
install [LJA](https://github.com/AntonBankevich/LJA), as is discussed in the
strainFlye README. (That said: as of writing strainFlye's tests do not rely on
LJA being installed, so this step is optional.)

## Development commands

strainFlye's
[`Makefile`](https://github.com/fedarko/strainFlye/blob/main/Makefile),
in the root of its repository, provides a few useful commands for development.
These are described here. All of these commands assume that strainFlye has been
installed as described above, that `Make` is installed, and that you are
running these commands from the root of strainFlye's repository.

### Running tests

strainFlye uses `pytest` and `pytest-cov` to run its automated tests.
You can run these tests using:

```bash
make test
```

Thanks to our use of `pytest-cov`, this command will create a few formats
of coverage report after running:

1. it'll display a summary on the terminal,
2. it'll create an XML file in the root of strainFlye's repository that is
   useful for passing coverage information to CodeCov,
3. and it'll create a directory in the root of strainFlye's repository
   called `htmlcov`.

The last of these (`htmlcov`) is probably the most useful for local
development.

### Linting, style-checking, and auto-formatting

strainFlye uses `flake8` to lint and style-check its source code, and `black`
to auto-format its source code.

To run linting and style-checking, you can use the following command:

```bash
make stylecheck
```

To auto-format the source code:

```bash
make style
```

Note that, in some rare cases, `flake8` and `black` can disagree. If you run
into trouble trying to get the style-checking to pass (or with other
development steps), feel free to open an issue in this repository.
