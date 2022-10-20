#!/usr/bin/env python
# This file adapted from https://github.com/biocore/qurro/blob/master/setup.py

from setuptools import find_packages, setup
from strainflye import __version__

classifier_str = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: MIT License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classifier_str.split("\n") if s]

description = (
    "Pipeline for analyzing rare mutations in metagenome-assembled genomes"
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
        "scikit-bio",
        "networkx",
        # This isn't a hard version requirement -- I'm not sure what the
        # absolute minimum version is (maybe lower, probably not higher).
        # This should be good enough.
        "matplotlib >= 3.0",
        # Also not a hard requirement.
        "pandas >= 1.0",
        # We use the min_open and max_open parameters of click.IntRange and
        # click.FloatRange, which need click >= 8.0; and we use the executable
        # parameter of click.Path, which needs click >= 8.1. If desired we
        # could avoid using these newer click features and re-implement the
        # checks ourselves (to avoid reliance on newer click versions), but I
        # don't think that is worth the extra effort right now.
        "click >= 8.1",
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
            "black",
        ]
    },
    classifiers=classifiers,
    entry_points={
        "console_scripts": [
            "strainFlye=strainflye._cli:strainflye",
            "strainflye=strainflye._cli:strainflye",
        ],
    },
    python_requires=">=3.6",
)
