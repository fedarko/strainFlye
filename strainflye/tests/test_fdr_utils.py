import os
import tempfile
import subprocess
import pytest
import pandas as pd
import strainflye.fdr_utils as fu
from io import StringIO
from strainflye.errors import ParameterError, SequencingDataError
from strainflye.config import DI_PREF
from .utils_for_testing import mock_log

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


def test_parse_bcf_good():
    # Header + first four mutations called on SheepGut using p = 0.15%
    # ... just for reference, using StringIO didn't seem to work with pysam's
    # BCF reader, hence the use of tempfiles here
    # NOTE FROM MARCUS: pysam accepts either VCF or BCF files, and -- since we
    # just recently switched from working with VCF to BCF, mainly -- these
    # tests still output VCF files.
    fh = write_vcf_tempfile(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(fh.name)

    # Let's get the easy checks out of the way first...
    assert thresh_type == "p"
    assert thresh_min == 15
    # This part of things -- parsing the BCF -- pysam handles, so we don't go
    # too in depth. Just verify it gets the positions, alleles, and info right.
    correct_pos = [43, 255, 356, 387]
    correct_mdp = [380, 385, 403, 395]
    correct_aad = [2, 2, 2, 2]
    correct_ref = ["A", "T", "A", "T"]
    correct_alt = ["T", "C", "T", "C"]
    for ri, rec in enumerate(bcf_obj.fetch()):
        assert rec.pos == correct_pos[ri]
        assert rec.info.get("MDP") == correct_mdp[ri]
        assert rec.ref == correct_ref[ri]

        # this part we'll have to change if we start calling multiple mutations
        # at a single position
        assert len(rec.alts) == 1
        assert rec.alts[0] == correct_alt[ri]

        # AAD is treated as a tuple containing one element by pysam, not as a
        # single value. This is because I've labelled the AAD INFO field with a
        # Number of "A", meaning that there is one of these per alternate
        # allele. Currently we don't call more than one mutation at a position,
        # so this doesn't do anything, but -- in case we modify this in the
        # future -- this leaves the door open for us to add this in.
        rec_aad = rec.info.get("AAD")
        assert len(rec_aad) == 1
        assert rec_aad[0] == correct_aad[ri]


def test_parse_bcf_multi_sf_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fh = write_vcf_tempfile(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_20, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    with pytest.raises(ParameterError) as ei:
        fu.parse_bcf(fh.name)
    exp_patt = (
        f"BCF file {fh.name} has multiple strainFlye threshold filter headers."
    )
    assert str(ei.value) == exp_patt


def test_parse_bcf_no_sf_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fh = write_vcf_tempfile(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    with pytest.raises(ParameterError) as ei:
        fu.parse_bcf(fh.name)
    exp_patt = (
        f"BCF file {fh.name} doesn't seem to be from strainFlye: no threshold "
        "filter headers."
    )
    assert str(ei.value) == exp_patt


def test_parse_bcf_missing_info_fields():
    # Header + first four mutations called on SheepGut using p = 0.15%
    # kind of a gross way of doing this -- should ideally set this up as a list
    # from the start, theeen join it up, rather than going from string to list
    # to string again. but whatevs
    text = (
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    split_text = text.splitlines()

    # remove just the MDP line
    fh = write_vcf_tempfile("\n".join(split_text[:4] + split_text[5:]))
    with pytest.raises(ParameterError) as ei:
        fu.parse_bcf(fh.name)
    exp_patt = f"BCF file {fh.name} needs to have MDP and AAD info fields."
    assert str(ei.value) == exp_patt

    # remove just the AAD line
    fh = write_vcf_tempfile("\n".join(split_text[:5] + split_text[6:]))
    with pytest.raises(ParameterError) as ei:
        fu.parse_bcf(fh.name)
    exp_patt = f"BCF file {fh.name} needs to have MDP and AAD info fields."
    assert str(ei.value) == exp_patt

    # remove both the MDP and AAD lines
    fh = write_vcf_tempfile("\n".join(split_text[:4] + split_text[6:]))
    with pytest.raises(ParameterError) as ei:
        fu.parse_bcf(fh.name)
    exp_patt = f"BCF file {fh.name} needs to have MDP and AAD info fields."
    assert str(ei.value) == exp_patt

    # finally, sanity check that putting back in all the lines works (doesn't
    # throw an error) -- we check this dataset in detail above, so no need to
    # do that again here
    fh = write_vcf_tempfile("\n".join(split_text))
    fu.parse_bcf(fh.name)


def test_check_decoy_selection():
    # "Good" cases
    assert fu.check_decoy_selection(None, "contig_name") == "DC"
    assert fu.check_decoy_selection("totally_real_file.tsv", None) == "DI"

    # terrible no good very bad cases
    with pytest.raises(ParameterError) as ei:
        fu.check_decoy_selection("file.tsv", "contig")
    assert str(ei.value) == (
        "Both the diversity indices file and a decoy contig are specified. "
        "These options are mutually exclusive."
    )

    with pytest.raises(ParameterError) as ei:
        fu.check_decoy_selection(None, None)
    assert str(ei.value) == (
        "Either the diversity indices file or a decoy contig must be "
        "specified."
    )


def test_autoselect_decoy_not_enough_contigs():
    # fortunately, pandas.read_csv() is cool with StringIO, so we don't have to
    # write out tempfiles
    tsv = StringIO(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}\t{DI_PREF}\n"
        "edge_1\t35.2\t100\t0.5\t0.2\n"
    )
    with pytest.raises(ParameterError) as ei:
        fu.autoselect_decoy(tsv, 1e6, 1e2, mock_log)
    assert str(ei.value) == "Diversity indices file describes < 2 contigs."


def test_autoselect_decoy_missing_required_cols():
    def run_check(tsv_text):
        with pytest.raises(ParameterError) as ei:
            fu.autoselect_decoy(StringIO(tsv_text), 1e6, 1e2, mock_log)
        assert str(ei.value) == (
            "Diversity indices file must include the Length and "
            "AverageCoverage columns."
        )

    # Check 1: omit AverageCoverage
    run_check(
        f"Contig\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35\t0.5\t0.2\n"
        "edge_2\t36\t0.5\t0.2\n"
    )

    # Check 2: omit Length
    run_check(
        f"Contig\tAverageCoverage\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35.2\t0.5\t0.2\n"
        "edge_2\t36\t0.5\t0.2\n"
    )

    # Check 3: omit both
    run_check(
        f"Contig\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t0.5\t0.2\n"
        "edge_2\t0.5\t0.2\n"
    )


def test_autoselect_decoy_no_passing():
    def run_check(tsv_text):
        with pytest.raises(SequencingDataError) as ei:
            fu.autoselect_decoy(StringIO(tsv_text), int(1e6), 100.5, mock_log)
        assert str(ei.value) == (
            "No contigs pass the min length \u2265 1,000,000 and min average "
            "cov \u2265 100.5x checks."
        )

    # Check 1: neither contig passes either the length or coverage check
    run_check(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35\t1000\t0.5\t0.2\n"
        "edge_2\t36\t1000\t0.6\t0.1\n"
    )
    # Check 2: One passes the coverage check, the other passes the length
    # check, neither passes both
    run_check(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t3500\t100\t0.5\t0.2\n"
        "edge_2\t36\t100000000000\t0.6\t0.1\n"
    )


def test_autoselect_decoy_only_one_passing(capsys):
    tsv_text = (
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35\t1000\t0.5\t0.2\n"
        "edge_2\t202.53\t10000000\t0.6\t0.1\n"
    )
    dc = fu.autoselect_decoy(StringIO(tsv_text), int(1e7), 202.53, mock_log)
    assert dc == "edge_2"
    captured = capsys.readouterr()
    assert captured.out == (
        "MockLog: Warning: Only one contig passes the min length \u2265 "
        "10,000,000 and min average cov \u2265 202.53x checks. Selecting it.\n"
    )


def test_autoselect_decoy_all_passing_undefined_di():
    def run_check(tsv_text):
        with pytest.raises(SequencingDataError) as ei:
            fu.autoselect_decoy(StringIO(tsv_text), int(1e7), 202.53, mock_log)
        assert str(ei.value) == (
            "No diversity index column has at least two contigs that (1) pass "
            "the min length \u2265 10,000,000 and min average cov \u2265 "
            "202.53x checks and (2) have defined and distinct diversity "
            "indices in this column."
        )

    # Check 1: One NA in both columns
    run_check(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35000\t100000000000\tNA\t0.2\n"
        "edge_2\t202.53\t10000000\t0.6\tNA\n"
    )
    # Check 2: All NAs
    run_check(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35000\t100000000000\tNA\tNA\n"
        "edge_2\t202.53\t10000000\tNA\tNA\n"
    )


def test_autoselect_decoy_all_passing_identical_di():
    # The error message received is the same as the "all passing undefined"
    # case tested above, but the cause is slightly different.
    # I could actually see this happening in practice (not likely, but if there
    # are a lot of small contigs / the users adjust the threshold params /
    # etc), so we should test it.
    def run_check(tsv_text):
        with pytest.raises(SequencingDataError) as ei:
            fu.autoselect_decoy(StringIO(tsv_text), int(1e7), 202.53, mock_log)
        assert str(ei.value) == (
            "No diversity index column has at least two contigs that (1) pass "
            "the min length \u2265 10,000,000 and min average cov \u2265 "
            "202.53x checks and (2) have defined and distinct diversity "
            "indices in this column."
        )

    run_check(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35000\t100000000000\t0\t3\n"
        "edge_2\t202.53\t10000000\t0\t3\n"
        "edge_3\t202.53\t10000000\t0\t3\n"
    )
    run_check(
        f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
        "edge_1\t35000\t100000000000\tNA\t3\n"
        "edge_2\t202.53\t10000000\tNA\t3\n"
        "edge_3\t202.53\t10000000\tNA\t3\n"
    )


def test_autoselect_decoy_good():
    # Total scores:
    #
    # A: 1 + 0 = 1 (the 1 is because its div idx is undefined in col 1)
    #
    # B: 0 + 1 = 1
    #
    # C:  ((0.01 / 0.6) = 0.01666)
    #   + (((0.2 - 0.14) / (0.31 - 0.14)) = 0.3529) = 0.3696
    #
    # D:  ((0.1 / 0.6) = 0.1666)
    #   + (((0.15 - 0.14) / (0.31 - 0.14)) = 0.3529) = 0.2255
    #
    # E: 1 + 0 = 1
    #
    # F: No score (length doesn't pass)
    #
    # The "D" contig should be selected. It doesn't have the lowest diversity
    # index in either column, but it's close enough to the bottom.
    assert (
        fu.autoselect_decoy(
            StringIO(
                f"Contig\tAverageCoverage\tLength\t{DI_PREF}1\t{DI_PREF}2\n"
                "A\t35000\t100000000000\tNA\t0.14\n"
                "B\t202.53\t10000000\t0.0\t0.31\n"
                "C\t202.53\t10000000\t0.01\t0.2\n"
                "D\t202.53\t10000000\t0.1\t0.15\n"
                "E\t202.53\t10000000\t0.6\tNA\n"
                "F\t202.53\t5\t0\t0\n"
            ),
            1000,
            100,
            mock_log,
        )
        == "D"
    )


def test_normalize_series_identical():
    assert (
        fu.normalize_series(
            pd.Series([3, 3, 3], index=["a", "b", "c"], name="DI")
        )
        is None
    )

    # in practice the decoy autoselection should never pass a series with < 2
    # elements but you never know, may as well be safe
    assert fu.normalize_series(pd.Series([3], index=["a"], name="DI")) is None


def test_normalize_series_good():
    s_params = {"index": ["a", "b", "c", "d", "e"], "name": "DI"}
    assert (
        (fu.normalize_series(pd.Series([2, 5, 9, 2, 1], **s_params)))
        == pd.Series([1 / 8, 1 / 2, 1, 1 / 8, 0], **s_params)
    ).all()

    s_params = {"index": ["a", "b"], "name": "DI"}
    assert (
        (fu.normalize_series(pd.Series([2, 5], **s_params)))
        == pd.Series([0, 1], **s_params)
    ).all()


def test_compute_full_decoy_contig_mut_rates_p_simple():
    # This example is designed to be simple -- the MDPs are all exactly 10,000
    # to make it simple to see what AADs will trigger a p-mutation.
    bcf_name = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=500>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=10000;AAD=1\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=10000;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=10000;AAD=16\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=10000;AAD=19\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(bcf_name)
    mut_rates = fu.compute_full_decoy_contig_mut_rates(
        bcf_obj, thresh_type, range(15, 21), "edge_1", 500
    )
    denominator = 3 * 500
    try:
        # There are two p-mutations at p = 0.15% and p = 0.16%;
        # one p-mutation at p = 0.17%, p = 0.18%, and p = 0.19%;
        # and no p-mutations at p = 0.2%.
        assert mut_rates == [
            2 / denominator,
            2 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            0,
        ]
    finally:
        # Python should automatically clear the .vcf tempfile we created
        # earlier from the system, regardless of how this test goes, but
        # it won't do this for the .bcf and .bcf.csi files we generated using
        # bcftools view and bcftools index. So we clear these ourselves, to
        # avoid littering the cluster with this nonsense!
        os.remove(bcf_name)
        os.remove(bcf_name + ".csi")


def test_compute_full_decoy_contig_mut_rates_r_simple():
    bcf_name = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220608\n"
        '##source="strainFlye v0.0.1: r-mutation calling (--min-r = 5)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=500>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minr_5, Description="min r threshold">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=10000;AAD=1\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=10000;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=10000;AAD=5\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=10000;AAD=10\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(bcf_name)
    mut_rates = fu.compute_full_decoy_contig_mut_rates(
        bcf_obj, thresh_type, range(5, 13), "edge_1", 500
    )
    denominator = 3 * 500
    try:
        # Two r-mutations at r = 5, then only one r-mutation for 6 <= r <= 10,
        # then zero r-mutations from then up.
        assert mut_rates == [
            2 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            0,
            0,
        ]
    finally:
        os.remove(bcf_name)
        os.remove(bcf_name + ".csi")


def test_compute_number_of_mutations_in_full_contig_p_with_indisputable():
    bcf_name = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=999>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t100\t.\tA\tT\t.\t.\tMDP=100000;AAD=5\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=100000;AAD=151\n"
        "edge_1\t300\t.\tA\tT\t.\t.\tMDP=100000;AAD=4899\n"
        "edge_1\t301\t.\tA\tT\t.\t.\tMDP=100000;AAD=4900\n"
        "edge_1\t302\t.\tA\tT\t.\t.\tMDP=100000;AAD=4995\n"
        "edge_1\t303\t.\tA\tT\t.\t.\tMDP=100000;AAD=4999\n"
        "edge_1\t304\t.\tA\tT\t.\t.\tMDP=100000;AAD=5001\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=100000;AAD=5000\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(bcf_name)
    num_muts = fu.compute_number_of_mutations_in_full_contig(
        bcf_obj, thresh_type, range(15, 500), "edge_1"
    )
    try:
        assert len(num_muts) == 485
        # Ignore the two indisputable mutations (AAD = 5000 and 5001, which
        # come out to a frequency of 5,000 / 100,000 >= 5%). Also, ignore the
        # mutation with AAD = 5, which is below p = 0.15% (in practice we
        # should never see these sorts of mutations since the naive caller
        # wouldn't include them, but this checks that they're implicitly
        # ignored).
        assert num_muts[0] == 5

        # Once we get up to p = 16 (aka 0.16%), the AAD = 151 mutation (0.151%
        # frequency) no longer is counted. We stay at 4 mutations up until we
        # get to i = 475, aka p = 490, at which point we'll drop the AAD = 4899
        # (4.899% frequency) mutation.
        for i in range(1, 475):
            assert num_muts[i] == 4
        assert num_muts[475] == 3
        # When we get up to i = 476 (p = 491), we need to also drop the AAD =
        # 4900 (4.9% frequency) mutation. All that's left now are the AAD =
        # 4995 and 4999 mutations.
        for i in range(476, 485):
            assert num_muts[i] == 2
    finally:
        # Python should automatically clear the .vcf tempfile we created
        # earlier from the system, regardless of how this test goes, but
        # it won't do this for the .bcf and .bcf.csi files we generated using
        # bcftools view and bcftools index. So we clear these ourselves, to
        # avoid littering the cluster with this nonsense!
        os.remove(bcf_name)
        os.remove(bcf_name + ".csi")


def test_compute_number_of_mutations_in_full_contig_thresh_val_errors():
    bcf_name = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220608\n"
        '##source="strainFlye v0.0.1: r-mutation calling (--min-r = 5)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=500>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minr_5, Description="min r threshold">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=10000;AAD=1\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=10000;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=10000;AAD=5\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=10000;AAD=10\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(bcf_name)
    try:
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_full_contig(
                bcf_obj, thresh_type, range(5, 13, 2), "edge_1"
            )
        assert str(ei.value) == "thresh_vals must use a step size of 1."
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_full_contig(
                bcf_obj, thresh_type, range(15, 5, -1), "edge_1"
            )
        assert str(ei.value) == "thresh_vals must use a step size of 1."
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_full_contig(
                bcf_obj, thresh_type, range(5, 5), "edge_1"
            )
        assert str(ei.value) == "thresh_vals must have a positive length."
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_full_contig(
                bcf_obj, thresh_type, range(-9, -1), "edge_1"
            )
        assert str(ei.value) == "thresh_vals' start and stop must be positive."
    finally:
        os.remove(bcf_name)
        os.remove(bcf_name + ".csi")


def test_compute_target_contig_fdr_curve_info():
    bcf_name = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220608\n"
        '##source="strainFlye v0.0.1: r-mutation calling (--min-r = 5)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=100>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minr_5, Description="min r threshold">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=10000;AAD=1\n"
        "edge_1\t50\t.\tT\tC\t.\t.\tMDP=10000;AAD=2\n"
        "edge_1\t51\t.\tA\tT\t.\t.\tMDP=10000;AAD=5\n"
        "edge_1\t90\t.\tT\tC\t.\t.\tMDP=10000;AAD=10\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(bcf_name)
    try:
        fdr_line, num_line = fu.compute_target_contig_fdr_curve_info(
            bcf_obj,
            thresh_type,
            range(5, 12),
            "edge_1",
            100,
            [0.001, 0.002, 0.003, 0.004, 0.1, 0.006, 0.007],
        )
        assert (
            fdr_line == "edge_1\t15.0\t60.0\t90.0\t120.0\t3000.0\t180.0\tNA\n"
        )
        assert num_line == (
            "edge_1\t20000.0\t10000.0\t10000.0\t10000.0\t10000.0\t10000.0\t0\n"
        )
    finally:
        os.remove(bcf_name)
        os.remove(bcf_name + ".csi")


def test_compute_decoy_contig_mut_rates_full_r_simple():
    # this is the same as test_compute_full_decoy_contig_mut_rates_r_simple(),
    # but it calls compute_decoy_contig_mutation_rates() (which, in turn, calls
    # compute_full_decoy_contig_mut_rates(), which calls
    # compute_number_of_mutations_in_full_contig()...)
    bcf_name = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220608\n"
        '##source="strainFlye v0.0.1: r-mutation calling (--min-r = 5)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=500>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minr_5, Description="min r threshold">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=10000;AAD=1\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=10000;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=10000;AAD=5\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=10000;AAD=10\n"
    )
    bcf_obj, thresh_type, thresh_min = fu.parse_bcf(bcf_name)
    # where we're going, we don't need other nucleotides
    contigs_file = StringIO(f">edge_1\n{'A' * 500}")
    mut_rates = fu.compute_decoy_contig_mut_rates(
        contigs_file, bcf_obj, thresh_type, range(5, 13), "edge_1", "Full"
    )
    denominator = 3 * 500
    try:
        # Two r-mutations at r = 5, then only one r-mutation for 6 <= r <= 10,
        # then zero r-mutations from then up.
        assert mut_rates == [
            2 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            1 / denominator,
            0,
            0,
        ]
    finally:
        os.remove(bcf_name)
        os.remove(bcf_name + ".csi")


def test_load_and_sanity_check_fdr_file_basic_errors():
    ep = "Input FDR TSV file seems malformed"
    # Very real and legitimate TSV file
    tsv = StringIO(f"OIJsdoijsjdlasidjlasdj\taiodsfjoaidsfj,a123")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert str(ei.value) == f'{ep}: no "Contig" header?'

    tsv = StringIO(f"Contig\tp15\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert str(ei.value) == f'{ep}: no contigs described?'

    tsv = StringIO(f"Contig\nedge_1\nedge_2\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert str(ei.value) == f'{ep}: no threshold values described?'


def test_load_and_sanity_check_fdr_file_p_vs_r_mismatch():
    ep = "Input FDR TSV file seems malformed"

    # Check p vs. r mismatch
    tsv_str = f"Contig\tp15\tp16\nedge_1\t3\t5\n"
    tsv = StringIO(tsv_str)
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "r")
    assert str(ei.value) == f'{ep}: columns should start with r.'

    # ... And check that using p instead of r as the thresh_type works.
    # (We need to make a new StringIO, because I guess the old one already got
    # used up?)
    tsv = StringIO(tsv_str)
    obs_df = fu.load_and_sanity_check_fdr_file(tsv, "p")
    exp_df = pd.DataFrame({"p15": [3], "p16": [5]}, index=pd.Index(["edge_1"], name="Contig"))
    pd.testing.assert_frame_equal(obs_df, exp_df)
