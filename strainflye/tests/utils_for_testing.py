import os
import tempfile
import contextlib
from strainflye.bcf_utils import compress_vcf


def mock_log(text, prefix="PREFIX\n"):
    print(f"{prefix}MockLog: {text}")


def write_vcf(text, delete=True):
    """Writes out text to a tempfile that can be used as a test VCF file.

    Parameters
    ----------
    text: str
        Text to be written out to the VCF file. We don't add or remove
        anything, so this string should cover the entirety of what you want to
        include in the file.

    delete: bool
        Will be passed to tempfile.NamedTemporaryFile() as the "delete"
        parameter there. Long story short, if you want this file to be
        automatically deleted when it's closed, keep this as True; if you plan
        to delete it yourself, then set this as False.

    Returns
    -------
    fh: "file-like object" (in practice, type tempfile._TemporaryFileWrapper)
        Represents this tempfile. For testing, you should probably use this as
        a context manager (see
        https://docs.python.org/3/library/tempfile.html#examples) to make it
        easy to automatically close (and thus delete) this file after a test.
    """
    # ... just for reference, using StringIO didn't seem to work with pysam's
    # BCF reader, hence the use of tempfiles here
    fh = tempfile.NamedTemporaryFile(suffix=".vcf", delete=delete)
    with open(fh.name, "w") as f:
        f.write(text)
    return fh


@contextlib.contextmanager
def write_indexed_bcf(vcf_text):
    """Writes out a tempfile that can be used as a test BCF file.

    Parameters
    ----------
    text: str
        Text to be written out to a VCF file. This VCF file will then be
        converted to BCF and indexed. (The "middle-man" VCF file will
        automatically be deleted.)

    Returns
    -------
    fh: "file-like object" (in practice, type tempfile._TemporaryFileWrapper)
        Represents this tempfile. Usable as a context manager; see write_vcf()
        docs for more details.

        I guess this technically "yields" this fh to the user??? But apparently
        these docstrings should still say "returns"
        (https://stackoverflow.com/a/39962779). But, like, idk. Look I still
        think yield statements are magic. They already took Santa Claus and the
        Easter Bunny from me, please just let me have this innocence of not
        knowing how yields work.

    References
    -----
    See
    https://docs.python.org/3/library/contextlib.html#contextlib.contextmanager
    for some details on context managers. Luckily, the basic context manager
    functionality worked pretty well out of the box for letting us easily clean
    up *.bcf.csi and *.vcf files.
    """
    # When we're testing stuff like mutation rate computation (for which a VCF
    # / BCF index needs to be present), we actually need to write out a BCF and
    # index it -- we can't just use a VCF file instead. Notably: the ##contig
    # line needs to be there, otherwise we can't convert to BCF.
    #
    # (We set delete to False because compress_vcf() will delete the VCF file
    # -- if we leave delete as True then we get an error when the tempfile
    # module later tries to automatically delete the already-deleted VCF.
    # However, we leave delete as True for the BCF so that it will be deleted
    # at the end of the test for which this function is called. Update, never
    # mind, we don't actually close these files, so they aren't getting deleted
    # anyway -- need to fix this -- see
    # https://github.com/fedarko/strainFlye/issues/38.)
    try:
        vcf_fh = write_vcf(vcf_text, delete=False)
        bcf_fh = tempfile.NamedTemporaryFile(suffix=".bcf")
        compress_vcf(vcf_fh.name, bcf_fh.name, mock_log)
        yield bcf_fh
    finally:
        # Clean up tempfiles that may or may not exist. (Vulnerable to race
        # conditions of files being deleted then added back or some nonsense
        # like that, but that shouldn't matter unless someone is trying to
        # break these tests intentionally.)

        # The *.vcf file should have been deleted during compress_vcf(), but if
        # an error popped up than this file may still be around.
        if os.path.exists(vcf_fh.name):
            os.remove(vcf_fh.name)

        # The *.bcf file (represented by bcf_fh.name) should already be
        # removed (or soon to be removed), since we yield bcf_fh to the caller
        # (and they should be using it as its own context manager, as
        # documented in the tempfile API.)
        #
        # However, we gotta take care of the *.bcf.csi file created during
        # indexing. (Ordinarily, this file should exist, but it may not exist
        # if for example we ran into an error before indexing. So we check,
        # just to be safe.)
        idx_name = bcf_fh.name + ".csi"
        if os.path.exists(idx_name):
            os.remove(idx_name)
