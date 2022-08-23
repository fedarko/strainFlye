import tempfile
from strainflye.bcf_utils import compress_vcf


def mock_log(text, prefix="PREFIX\n"):
    print(f"{prefix}MockLog: {text}")


def write_vcf_tempfile(text, delete=True):
    # ... just for reference, using StringIO didn't seem to work with pysam's
    # BCF reader, hence the use of tempfiles here
    fh = tempfile.NamedTemporaryFile(suffix=".vcf", delete=delete)
    with open(fh.name, "w") as f:
        f.write(text)
    return fh


def write_indexed_bcf(vcf_text):
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
    fh = write_vcf_tempfile(vcf_text, delete=False)
    bcf_name = tempfile.NamedTemporaryFile(suffix=".bcf").name
    compress_vcf(fh.name, bcf_name, mock_log)
    return bcf_name
