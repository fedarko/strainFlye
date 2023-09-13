#!/usr/bin/env python
# This file adapted from https://github.com/biocore/qurro/blob/master/setup.py

from setuptools import find_packages, setup
from strainflye import __version__

classifier_str = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classifier_str.split("\n") if s]

description = (
    "Pipeline for analyzing rare mutations in metagenomes assembled using "
    "long and accurate reads"
)

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="strainFlye",
    version=__version__,
    license="BSD",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Marcus Fedarko",
    author_email="mfedarko@ucsd.edu",
    maintainer="Marcus Fedarko",
    maintainer_email="mfedarko@ucsd.edu",
    url="https://github.com/fedarko/strainFlye",
    packages=find_packages(),
    install_requires=[
        "scikit-bio >= 0.5.8",
        "networkx",
        # Also not a hard requirement.
        "pandas >= 1.0",
        # We use the min_open and max_open parameters of click.IntRange and
        # click.FloatRange, which need click >= 8.0. We explicitly avoid
        # needing to rely on click >= 8.1.0 (e.g. using the newly-added
        # "executable" parameter of click.Path), because click 8.1.0 removes
        # support for python 3.6.
        "click >= 8.0",
    ],
    setup_requires=[
        "cython",
        "numpy",
        "scipy",
        "pysam",
        "pysamstats",
    ],
    # Based on how Altair splits up its requirements:
    # https://github.com/altair-viz/altair/blob/master/setup.py
    extras_require={
        "dev": [
            "pytest >= 4.2",
            "pytest-cov >= 2.0",
            "flake8",
            # black 22.10 stopped supporting python 3.6; ideally we'd adjust
            # our CI to use later versions of black on the 3.7 build, but as a
            # hacky solution pinning black is fine
            "black < 22.10",
            # This isn't a hard version requirement -- I'm not sure what the
            # absolute minimum version is (maybe lower, probably not higher).
            # This should be good enough. (Needed for tutorial notebooks.)
            "matplotlib >= 3.0",
        ]
    },
    classifiers=classifiers,
    entry_points={
        "console_scripts": [
            "strainFlye=strainflye._cli:strainflye",
            "strainflye=strainflye._cli:strainflye",
        ],
    },
    python_requires=">=3.6,<3.8",
)
