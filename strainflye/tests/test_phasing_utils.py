import os
import pytest
from io import StringIO
import strainflye.phasing_utils as pu
from strainflye.errors import ParameterError
from .utils_for_testing import mock_log, write_indexed_bcf


TEST_DIR = os.path.join(
    "strainflye",
    "tests",
    "inputs",
    "small",
)
FASTA = os.path.join(TEST_DIR, "contigs.fasta")
BAM = os.path.join(TEST_DIR, "alignment.bam")
BCF = os.path.join(TEST_DIR, "call-r-min3-di12345", "naive-calls.bcf")


def test_load_triplet_good(capsys):
    n2l, bam_obj, bcf_obj = pu.load_triplet(FASTA, BAM, BCF, mock_log)
    # The actual fasta / bam / bcf loading done in this file isn't super
    # important to test, since these are either already tested elsewhere (the
    # fasta and bcf loading) or are completely done by an external library (bam
    # loading).
    #
    # So we just cursorily check that these loaded files look ok. the main
    # thing we wanna check is that the logged output is correct.
    assert n2l == {"c1": 23, "c2": 12, "c3": 16}
    exp_contigs = set(["c1", "c2", "c3"])
    assert set(bam_obj.references) == exp_contigs
    assert set(bcf_obj.header.contigs) == exp_contigs

    captured = capsys.readouterr()
    assert captured.out == (
        "PREFIX\nMockLog: Loading and checking FASTA, BAM, and BCF files...\n"
        "MockLog: The FASTA file describes 3 contig(s).\n"
        "MockLog: All FASTA contig(s) are included in the BAM file (this BAM "
        "file has 3 reference(s)).\n"
        "MockLog: All FASTA contig(s) are included in the BCF file (the "
        "header of this BCF file describes 3 contig(s)).\n"
        "MockLog: The lengths of all contig(s) in the FASTA file match the "
        "corresponding lengths in the BAM and BCF files.\n"
        "MockLog: So far, these files seem good.\n"
    )


def test_load_triplet_fasta_not_in_bcf_all():
    # This is actually a valid VCF file (we use it in the parse_sf_bcf()
    # tests); the problem with it is that it doesn't describe c1, c2, or c3.
    with write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=1000>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    ) as fh:
        with pytest.raises(ParameterError) as ei:
            pu.load_triplet(FASTA, BAM, fh.name, mock_log)
        assert str(ei.value) == (
            "All contigs in the FASTA file must also be contained in the BCF "
            "file."
        )


def test_load_triplet_fasta_not_in_bcf_partial():
    with write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=c2,length=12>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "c2\t1\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "c2\t4\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "c2\t5\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "c2\t9\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    ) as fh:
        with pytest.raises(ParameterError) as ei:
            pu.load_triplet(FASTA, BAM, fh.name, mock_log)
        assert str(ei.value) == (
            "All contigs in the FASTA file must also be contained in the BCF "
            "file."
        )


def test_load_triplet_fasta_not_in_bam_all():
    # The BAM check should get performed before the BCF check, so we can test
    # the case where FASTA contigs aren't in the BAM by just using a different
    # FASTA file
    with pytest.raises(ParameterError) as ei:
        pu.load_triplet(StringIO(">missing_contig\nACGT"), BAM, BCF, mock_log)
    assert str(ei.value) == (
        "All contigs in the FASTA file must also be contained in the BAM "
        "file."
    )
