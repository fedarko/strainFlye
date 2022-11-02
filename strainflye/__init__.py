"""
strainFlye is a pipeline for analyzing rare mutations in metagenomes.

(It's mainly intended to be used as a command-line tool, but I guess it could
be useful to import it as a Python module if you'd like to use some of its
utility functions.)
"""

import os


__version__ = "0.2.0"

# Kind of a hack -- make it so that all python files in this directory, except
# those starting with _, can be imported when we say "from strainflye import *"
# See https://stackoverflow.com/a/1057534
neighbors = os.listdir(os.path.dirname(__file__))
__all__ = [
    m[:-3] for m in neighbors if m.endswith(".py") and not m.startswith("_")
]
