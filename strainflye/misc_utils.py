import os


def make_output_dir(output_dir):
    """Creates an output directory, if it doesn't already exist."""
    os.makedirs(output_dir, exist_ok=True)
