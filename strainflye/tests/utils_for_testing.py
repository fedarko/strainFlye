import tempfile
import subprocess


def mock_log(text, prefix="PREFIX\n"):
    print(f"{prefix}MockLog: {text}")


def write_vcf_tempfile(text):
    fh = tempfile.NamedTemporaryFile(suffix=".vcf")
    with open(fh.name, "w") as f:
        f.write(text)
    return fh


def write_indexed_bcf(vcf_text):
    # When we're testing stuff like mutation rate computation (for which a VCF
    # / BCF index needs to be present), we actually need to write out a BCF and
    # index it -- we can't just use a VCF file instead. Notably: the ##contig
    # line needs to be there, otherwise we can't convert to BCF.
    fh = write_vcf_tempfile(vcf_text)
    bcf_name = fh.name[:-4] + ".bcf"
    subprocess.run(["bcftools", "view", "-O", "b", fh.name, "-o", bcf_name])
    subprocess.run(["bcftools", "index", bcf_name])
    return bcf_name
