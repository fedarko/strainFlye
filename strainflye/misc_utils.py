import os
from strainflye.errors import ParameterError


def make_output_dir(output_dir):
    """Creates an output directory, if it doesn't already exist."""
    os.makedirs(output_dir, exist_ok=True)


def verify_contigs_subset(child, parent, child_desc, parent_desc, exact=False):
    """Verifies that one set of contig names is a subset of another set.

    Parameters
    ----------
    child: set
        Set of contig names.

    parent: set
        Set of contig names.

    child_desc: str
        Human-readable description of the child set of contig names.

    parent_desc: str
        Human-readable description of the parent set of contig names.

    exact: bool
        If True, this makes the check stricter -- it now ensures that the two
        sets are identical (so parent must also be a subset of child).

    Returns
    -------
    None

    Raises
    ------
    ParameterError
        - If the child set is not a subset of the parent set.
        - If exact is True, and the parent set and child set are not identical.
        The resulting error message will include child_desc and parent_desc.
    """
    if not child.issubset(parent):
        raise ParameterError(
            f"All contigs in {child_desc} must also be contained in "
            f"{parent_desc}."
        )
    # Equivalently, we could check "exact and not parent.issubset(child)", but
    # I think this way of writing it is clearer
    if exact and parent != child:
        raise ParameterError(
            f"All contigs in {parent_desc} must also be contained in "
            f"{child_desc}."
        )
