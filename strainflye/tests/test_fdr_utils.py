import os
import pytest
import pandas as pd
import numpy as np
import strainflye.fdr_utils as fu
import strainflye.bcf_utils as bu
from io import StringIO
from strainflye.errors import ParameterError, SequencingDataError
from strainflye.config import DI_PREF
from .utils_for_testing import mock_log, write_indexed_bcf


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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(bcf_name)
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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(bcf_name)
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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(bcf_name)
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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(bcf_name)
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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(bcf_name)
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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(bcf_name)
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
    tsv = StringIO("OIJsdoijsjdlasidjlasdj\taiodsfjoaidsfj,a123")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert str(ei.value) == f'{ep}: no "Contig" header?'

    tsv = StringIO("Contig\tp15\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert str(ei.value) == f"{ep}: no contigs described?"

    tsv = StringIO("Contig\nedge_1\nedge_2\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert str(ei.value) == f"{ep}: no threshold values described?"


def test_load_and_sanity_check_fdr_file_p_vs_r_mismatch():
    ep = "Input FDR TSV file seems malformed"

    # Check p vs. r mismatch
    tsv_str = "Contig\tp15\tp16\nedge_1\t3\t5\n"
    tsv = StringIO(tsv_str)
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "r")
    assert str(ei.value) == f"{ep}: columns should start with r."

    # ... And check that using p instead of r as the thresh_type works.
    # (We need to make a new StringIO, because I guess the old one already got
    # used up?)
    tsv = StringIO(tsv_str)
    obs_df = fu.load_and_sanity_check_fdr_file(tsv, "p")
    exp_df = pd.DataFrame(
        {"p15": [3], "p16": [5]}, index=pd.Index(["edge_1"], name="Contig")
    )
    pd.testing.assert_frame_equal(obs_df, exp_df)


def test_load_and_sanity_check_fdr_file_weird_col_errors():
    ep = "Input FDR TSV file seems malformed"

    # Check that p / r vals should increase
    # p
    tsv = StringIO("Contig\tp15\tp14\nedge_1\t3\t5\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert (
        str(ei.value)
        == f"{ep}: values of p should increase from left to right."
    )

    # r
    tsv = StringIO("Contig\tr15\tr14\nedge_1\t3\t5\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "r")
    assert (
        str(ei.value)
        == f"{ep}: values of r should increase from left to right."
    )

    # Check that p / r vals should increase ins teps of 1
    # p
    tsv = StringIO("Contig\tp15\tp17\nedge_1\t3\t5\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "p")
    assert (
        str(ei.value)
        == f"{ep}: values of p should only increase in steps of 1."
    )

    # r
    tsv = StringIO("Contig\tr50\tr52\nedge_1\t3\t5\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "r")
    assert (
        str(ei.value)
        == f"{ep}: values of r should only increase in steps of 1."
    )


def test_load_and_sanity_check_fdr_file_invalid_fdrs():
    ep = "Input FDR TSV file seems malformed"
    tsv = StringIO("Contig\tr50\tr51\nedge_1\t3\t-1\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "r")
    assert (
        str(ei.value) == f"{ep}: Column r51 contains negative FDR estimates?"
    )

    ep = "Input FDR TSV file seems malformed"
    tsv = StringIO("Contig\tr50\tr51\nedge_1\tHotdog\t1\n")
    with pytest.raises(ParameterError) as ei:
        fu.load_and_sanity_check_fdr_file(tsv, "r")
    assert str(ei.value) == f"{ep}: Column r50 doesn't seem to be numeric?"


def test_load_and_sanity_check_fdr_file_normal():
    tsv = StringIO(
        "Contig\tr50\tr51\tr52\n"
        "edge_1\t0.5\t0.3\t0.2\n"
        "edge_2\t0.2\t0.1\t0.4\n"
        "edge_3\t1.5\tNA\tNA\n"
    )
    obs_df = fu.load_and_sanity_check_fdr_file(tsv, "r")
    exp_df = pd.DataFrame(
        {
            "r50": [0.5, 0.2, 1.5],
            "r51": [0.3, 0.1, np.nan],
            "r52": [0.2, 0.4, np.nan],
        },
        index=pd.Index(["edge_1", "edge_2", "edge_3"], name="Contig"),
    )
    pd.testing.assert_frame_equal(obs_df, exp_df)


def test_get_optimal_threshold_values():
    contig_idx = pd.Index(["edge_A", "edge_B", "edge_C"], name="Contig")
    # Please don't try to actually interpret this data, I made it up for these
    # tests
    fdr_df = pd.DataFrame(
        {
            "r5": [0.6, 0.61, 3],
            "r6": [1, 0.7, 0.8],
            "r7": [3, np.nan, 1],
            "r8": [4, 0.2, np.nan],
            "r9": [2, 0.1, np.nan],
            "r10": [1, 0.05, np.nan],
        },
        index=contig_idx,
    )

    obs_otv_10 = fu.get_optimal_threshold_values(fdr_df, 0.6)
    exp_otv_10 = pd.Series([5, 8, np.nan], index=contig_idx)
    pd.testing.assert_series_equal(obs_otv_10, exp_otv_10)

    obs_otv_10 = fu.get_optimal_threshold_values(fdr_df, 2)
    exp_otv_10 = pd.Series([5, 5, 6], index=contig_idx)
    pd.testing.assert_series_equal(obs_otv_10, exp_otv_10)

    # "Case 3": all of the FDRs are <= 10%
    obs_otv_10 = fu.get_optimal_threshold_values(fdr_df, 10)
    exp_otv_10 = pd.Series([5, 5, 5], index=contig_idx)
    pd.testing.assert_series_equal(obs_otv_10, exp_otv_10)

    # "Case 4": all of the FDRs are > 0.001%
    obs_otv_10 = fu.get_optimal_threshold_values(fdr_df, 0.001)
    exp_otv_10 = pd.Series([np.nan, np.nan, np.nan], index=contig_idx)
    pd.testing.assert_series_equal(obs_otv_10, exp_otv_10)


def test_log_optimal_threshold_value_stats_good(capsys):
    otv = pd.Series(
        [np.nan, 2, 5, 10, 3, 2, 9, 1, np.nan, 19],
        index=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],
    )
    fu.log_optimal_threshold_value_stats(otv, "p", 1, 500, 1.0, mock_log)
    captured = capsys.readouterr()
    assert captured.out == (
        "MockLog: For 8 / 10 contigs, there exist values of p (at least, "
        "considering the range from p = 1 to p = 500) that yield estimated "
        "FDRs \u2264 1.0%.\nMockLog: These values range from p = 1 (H) to p "
        "= 19 (J).\nMockLog: The mean of these values is p = 6.38.\n"
    )


def test_log_optimal_threshold_value_stats_same_thresh_extrema(capsys):
    otv = pd.Series(
        [np.nan, 5, 5, np.nan, 5, 5], index=["A", "B", "C", "D", "E", "F"]
    )
    fu.log_optimal_threshold_value_stats(otv, "r", 5, 5, 12.5, mock_log)
    captured = capsys.readouterr()
    # Yeah, this log message looks weird, but it's correct. I don't think it's
    # worth the effort to make a new case for this pathological scenario
    assert captured.out == (
        "MockLog: For 4 / 6 contigs, there exist values of r (at least, "
        "considering the range from r = 5 to r = 5) that yield estimated "
        "FDRs \u2264 12.5%.\nMockLog: These values range from r = 5 (B) to r "
        "= 5 (B).\nMockLog: The mean of these values is r = 5.00.\n"
    )


def test_log_optimal_threshold_value_stats_allnan(capsys):
    otv = pd.Series([np.nan, np.nan, np.nan], index=["A", "B", "C"])
    fu.log_optimal_threshold_value_stats(otv, "r", 1, 500, 1.0, mock_log)
    captured = capsys.readouterr()
    assert captured.out == (
        "MockLog: Warning: No values of r resulted in estimated FDRs \u2264 "
        "the fixed FDR, for all 3 contigs.\n"
    )
