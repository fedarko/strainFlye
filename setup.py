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

setup_reqs = []
with open("pip_requirements.txt", "r") as f:
    for line in f:
        if not line.startswith("#"):
            setup_reqs.append(line.strip())

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
        # Needed by recent versions of black, apparently:
        # https://pythonissues.com/issues/2903160 (click 7.0 results in the
        # call to black on GitHub Actions failing with this error)
        "click >= 8.0",
    ],
    setup_requires=setup_reqs,
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
