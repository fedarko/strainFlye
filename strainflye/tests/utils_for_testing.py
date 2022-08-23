import tempfile
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


def write_indexed_bcf(vcf_text, delete_vcf=False):
    """Writes out a tempfile that can be used as a test BCF file.

    Parameters
    ----------
    text: str
        Text to be written out to a VCF file. This VCF file will then be
        converted to BCF and indexed. (The "middle-man" VCF file will
        automatically be deleted.)

    delete_vcf: bool
        Will be passed to write_vcf() as the delete parameter there. This
        should usually be False (because when we call compress_vcf() on the VCF
        file, it'll be deleted); however, in the silly case where we expect
        compress_vcf() to fail before closing the VCF file, it makes sense to
        keep this as True. (Sorry, this is inelegant, but I'm in a rush.)

    Returns
    -------
    fh: "file-like object" (in practice, type tempfile._TemporaryFileWrapper)
        Represents this tempfile. Usable as a context manager; see write_vcf()
        docs for more details.

    Notes
    -----
    We need to figure out how to delete the *.bcf.csi files created from
    indexing... Maybe figure out how to bind their deletion to the closing of
    the BCF file?
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
    fh = write_vcf(vcf_text, delete=delete_vcf)
    bcf_fh = tempfile.NamedTemporaryFile(suffix=".bcf")
    compress_vcf(fh.name, bcf_fh.name, mock_log)
    return bcf_fh
