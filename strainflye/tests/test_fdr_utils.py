import os
import tempfile
import collections
import pytest
import skbio
import pysam
import pandas as pd
import numpy as np
import strainflye.fdr_utils as fu
import strainflye.bcf_utils as bu
from io import StringIO
from strainflye.errors import ParameterError, SequencingDataError, WeirdError
from strainflye.config import DI_PREF
from .utils_for_testing import mock_log, write_indexed_bcf


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BAM = os.path.join(IN_DIR, "alignment.bam")
BCF = os.path.join(IN_DIR, "call-r-min3-di12345", "naive-calls.bcf")
DI = os.path.join(IN_DIR, "call-r-min3-di12345", "diversity-indices.tsv")


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
    assert str(ei.value) == "Diversity indices file describes < 2 contig(s)."


def test_autoselect_decoy_missing_required_cols():
    def run_check(tsv_text):
        with pytest.raises(ParameterError) as ei:
            fu.autoselect_decoy(StringIO(tsv_text), 1e6, 1e2, mock_log)
        assert str(ei.value) == (
            'Diversity indices file must include the "Length" and '
            '"AverageCoverage" columns.'
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
    with write_indexed_bcf(
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
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        mut_rates = fu.compute_any_mutation_decoy_contig_mut_rates(
            bcf_obj, thresh_type, range(15, 21), "edge_1"
        )
        denominator = 3 * 500
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


def test_compute_full_decoy_contig_mut_rates_r_simple():
    with write_indexed_bcf(
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
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        mut_rates = fu.compute_any_mutation_decoy_contig_mut_rates(
            bcf_obj, thresh_type, range(5, 13), "edge_1"
        )
        denominator = 3 * 500
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


def test_compute_number_of_mutations_in_contig_p_with_indisputable():
    with write_indexed_bcf(
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
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        num_muts = fu.compute_number_of_mutations_in_contig(
            bcf_obj, thresh_type, range(15, 500), "edge_1"
        )
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


def test_compute_number_of_mutations_in_contig_thresh_val_errors():
    with write_indexed_bcf(
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
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_contig(
                bcf_obj, thresh_type, range(5, 13, 2), "edge_1"
            )
        assert str(ei.value) == "thresh_vals must use a step size of 1."
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_contig(
                bcf_obj, thresh_type, range(15, 5, -1), "edge_1"
            )
        assert str(ei.value) == "thresh_vals must use a step size of 1."
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_contig(
                bcf_obj, thresh_type, range(5, 5), "edge_1"
            )
        assert str(ei.value) == "thresh_vals must have a positive length."
        with pytest.raises(ParameterError) as ei:
            fu.compute_number_of_mutations_in_contig(
                bcf_obj, thresh_type, range(-9, -1), "edge_1"
            )
        assert str(ei.value) == "thresh_vals' start and stop must be positive."


def test_compute_target_contig_fdr_curve_info():
    with write_indexed_bcf(
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
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        ctx2fdr_lines, num_line = fu.compute_target_contig_fdr_curve_info(
            bcf_obj,
            thresh_type,
            range(5, 12),
            "edge_1",
            100,
            # mostly the same -- just checking that multiple fdr lines get
            # output
            {
                "butt": [0.001, 0.002, 0.003, 0.004, 0.1, 0.006, 0.007],
                "lol": [0.001, 0.002, 0.003, 0.004, 0.1, 0.5, 7],
            },
        )
        assert ctx2fdr_lines == {
            "butt": "edge_1\t15.0\t60.0\t90.0\t120.0\t3000.0\t180.0\tNA\n",
            "lol": "edge_1\t15.0\t60.0\t90.0\t120.0\t3000.0\t15000.0\tNA\n",
        }
        assert num_line == (
            "edge_1\t20000.0\t10000.0\t10000.0\t10000.0\t10000.0\t10000.0\t0\n"
        )


def test_compute_decoy_contig_mut_rates_full_r_simple():
    # this is the same as test_compute_full_decoy_contig_mut_rates_r_simple(),
    # but it calls compute_decoy_contig_mutation_rates() (which, in turn, calls
    # compute_full_decoy_contig_mut_rates(), which calls
    # compute_number_of_mutations_in_contig()...)
    with write_indexed_bcf(
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
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        # where we're going, we don't need other nucleotides
        contigs_file = StringIO(f">edge_1\n{'A' * 500}")
        ctx2mr = fu.compute_decoy_contig_mut_rates(
            contigs_file,
            pysam.AlignmentFile(BAM),
            bcf_obj,
            thresh_type,
            range(5, 13),
            "edge_1",
            ["Full"],
        )
        denominator = 3 * 500
        # Two r-mutations at r = 5, then only one r-mutation for 6 <= r <= 10,
        # then zero r-mutations from then up.
        assert ctx2mr == {
            "Full": [
                2 / denominator,
                1 / denominator,
                1 / denominator,
                1 / denominator,
                1 / denominator,
                1 / denominator,
                0,
                0,
            ]
        }


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


def test_run_estimate_both_dividx_and_decoy_specified():
    # Integration test version of a unit test from above
    with tempfile.TemporaryDirectory() as td:

        with pytest.raises(ParameterError) as ei:
            fu.run_estimate(
                FASTA,
                BAM,
                BCF,
                DI,
                "c1",
                # decoy contexts
                ["Full"],
                # high p
                500,
                # high r
                100,
                # decoy min length
                10,
                # decoy avg coverage
                5,
                # output dir
                td,
                mock_log,
            )

        assert str(ei.value) == (
            "Both the diversity indices file and a decoy contig are "
            "specified. These options are mutually exclusive."
        )


def check_fdr_and_num_dfs_r3_high10(td):
    out_fdr_fp = os.path.join(td, "fdr-Full.tsv")
    out_num_fp = os.path.join(td, "num-mutations-per-mb.tsv")

    assert os.path.exists(out_fdr_fp)
    assert os.path.exists(out_num_fp)

    # Test that the FDR estimates look good
    obs_fdr_df = pd.read_csv(out_fdr_fp, sep="\t", index_col=0)
    # Since we let strainFlye auto-select the decoy from this test data,
    # the contig "c2" will be selected. this contig has zero mutations, so
    # the FDRs for any of the target contigs will be zero (if this target
    # contig has at least 1 mutation at this r-threshold) or undefined /
    # NaN (if this target contig has no mutations at this r-threshold).
    #
    # Note that in practice we prooobably won't see decoy contigs with
    # exactly zero mutations due to the default length and cov thresholds
    # (knock on wood). this is just a pathological case that is useful for
    # testing.
    #
    # For reference: c1 has three mutated positions (one at r=3, two at
    # r=5), and c3 has two mutated positions (one at r=3, one at r=6).
    # So the NaNs start to show up in the FDR estimates after the max of
    # these r-values for each of these target contigs.
    exp_fdr_df = pd.DataFrame(
        {
            "r3": [0.0, 0.0],
            "r4": [0.0, 0.0],
            "r5": [0.0, 0.0],
            "r6": [np.nan, 0.0],
            "r7": [np.nan, np.nan],
            "r8": [np.nan, np.nan],
            "r9": [np.nan, np.nan],
        },
        index=pd.Index(["c1", "c3"], name="Contig"),
    )
    pd.testing.assert_frame_equal(obs_fdr_df, exp_fdr_df)

    # Test that the (# mutations / Mb) values look good
    # These are really high because these are short contigs, so they have
    # relatively high mutation rates.
    obs_num_df = pd.read_csv(out_num_fp, sep="\t", index_col=0)
    c1_num_coeff = 1e6 / 23
    c3_num_coeff = 1e6 / 16
    exp_num_df = pd.DataFrame(
        {
            "r3": [c1_num_coeff * 3, c3_num_coeff * 2],
            "r4": [c1_num_coeff * 2, c3_num_coeff * 1],
            "r5": [c1_num_coeff * 2, c3_num_coeff * 1],
            "r6": [0, c3_num_coeff * 1],
            "r7": [0, 0],
            "r8": [0, 0],
            "r9": [0, 0],
        },
        index=pd.Index(["c1", "c3"], name="Contig"),
    )
    pd.testing.assert_frame_equal(obs_num_df, exp_num_df)


def test_run_estimate_with_autoselect_and_full_decoy_good(capsys):
    with tempfile.TemporaryDirectory() as td:
        fu.run_estimate(
            FASTA,
            BAM,
            BCF,
            DI,
            None,
            ["Full"],
            # high p
            500,
            # high r (note that this is not inclusive; so, our FDR estimates
            # will only go up to r = 9)
            10,
            # decoy min length
            10,
            # decoy avg coverage
            5,
            # output dir
            td,
            mock_log,
        )
        check_fdr_and_num_dfs_r3_high10(td)

        # Finally, verify that the logged output looks good
        assert capsys.readouterr().out == (
            "PREFIX\nMockLog: Loading and checking FASTA, BAM, and BCF "
            "files...\n"
            "MockLog: The FASTA file describes 3 contig(s).\n"
            "MockLog: The FASTA file's contig(s) and BAM file's reference(s) "
            "match.\n"
            "MockLog: The FASTA file's contig(s) and BCF file's contig(s) "
            "match.\n"
            "MockLog: Also, the input BCF file describes r-mutations (minimum "
            "r = 3).\n"
            "MockLog: The lengths of all contig(s) in the FASTA file match "
            "the corresponding lengths in the BAM and BCF files.\n"
            "MockLog: So far, these files seem good.\n"
            "PREFIX\nMockLog: Selecting a decoy contig based on the diversity "
            "indices...\n"
            "MockLog: Selected c2 as the decoy contig.\n"
            "MockLog: Verified that this decoy contig is present in the "
            "FASTA, BAM, and BCF files.\n"
            "MockLog: (Sorry for doubting you.)\n"
            "PREFIX\nMockLog: Determining range of value(s) of r to "
            "consider...\n"
            "MockLog: We'll consider 7 value(s) of r: from 3 to 9.\n"
            "MockLog: r-mutations for r \u2265 10 will be considered "
            '"indisputable."\n'
            'MockLog: These "indisputable" mutations won\'t be included in '
            "the FDR estimation results.\n"
            "PREFIX\nMockLog: Computing mutation rates for c2 at these "
            "threshold values, for each of the 1 decoy context(s)...\n"
            "MockLog: Done.\n"
            "PREFIX\nMockLog: Computing mutation rates and FDR estimates for "
            "the 2 target contig(s)...\n"
            "MockLog: Done.\n"
        )


def test_run_estimate_tiny_chunk_size_good():
    # same as the above test, just with tiny chunk sizes to verify output
    # remains the same
    for cs in range(1, 15):
        with tempfile.TemporaryDirectory() as td:
            fu.run_estimate(
                FASTA,
                BAM,
                BCF,
                DI,
                None,
                ["Full"],
                500,
                10,
                10,
                5,
                td,
                mock_log,
                chunk_size=cs,
            )
            check_fdr_and_num_dfs_r3_high10(td)


def test_run_estimate_with_preselected_full_decoy_good(capsys):
    # same pre-selected decoy contig (c2) as is given by the above tests
    with tempfile.TemporaryDirectory() as td:
        fu.run_estimate(
            FASTA,
            BAM,
            BCF,
            None,
            "c2",
            ["Full"],
            # high p
            500,
            # high r
            10,
            # decoy min length
            10,
            # decoy avg coverage
            5,
            td,
            mock_log,
        )
        check_fdr_and_num_dfs_r3_high10(td)

        # logged output is a bit different now due to no longer using
        # decoy autoselection
        assert capsys.readouterr().out == (
            "PREFIX\nMockLog: Loading and checking FASTA, BAM, and BCF "
            "files...\n"
            "MockLog: The FASTA file describes 3 contig(s).\n"
            "MockLog: The FASTA file's contig(s) and BAM file's reference(s) "
            "match.\n"
            "MockLog: The FASTA file's contig(s) and BCF file's contig(s) "
            "match.\n"
            "MockLog: Also, the input BCF file describes r-mutations (minimum "
            "r = 3).\n"
            "MockLog: The lengths of all contig(s) in the FASTA file match "
            "the corresponding lengths in the BAM and BCF files.\n"
            "MockLog: So far, these files seem good.\n"
            "PREFIX\nMockLog: The specified decoy contig is c2.\n"
            "MockLog: Verified that this decoy contig is present in the "
            "FASTA, BAM, and BCF files.\n"
            "MockLog: (Sorry for doubting you.)\n"
            "PREFIX\nMockLog: Determining range of value(s) of r to "
            "consider...\n"
            "MockLog: We'll consider 7 value(s) of r: from 3 to 9.\n"
            "MockLog: r-mutations for r \u2265 10 will be considered "
            '"indisputable."\n'
            'MockLog: These "indisputable" mutations won\'t be included in '
            "the FDR estimation results.\n"
            "PREFIX\nMockLog: Computing mutation rates for c2 at these "
            "threshold values, for each of the 1 decoy context(s)...\n"
            "MockLog: Done.\n"
            "PREFIX\nMockLog: Computing mutation rates and FDR estimates for "
            "the 2 target contig(s)...\n"
            "MockLog: Done.\n"
        )


def test_run_estimate_small_high_thresholds():
    with tempfile.TemporaryDirectory() as td:
        with pytest.raises(ParameterError) as ei:
            fu.run_estimate(
                FASTA,
                BAM,
                BCF,
                DI,
                None,
                ["Full"],
                # high p
                500,
                # high r
                2,
                10,
                5,
                td,
                mock_log,
            )
        assert str(ei.value) == (
            "--high-r = 2 must be larger than the minimum r used in the BCF "
            "(3)."
        )
        assert len(os.listdir(td)) == 0

        with pytest.raises(ParameterError) as ei:
            fu.run_estimate(
                FASTA,
                BAM,
                os.path.join(IN_DIR, "call-p-min100", "naive-calls.bcf"),
                os.path.join(IN_DIR, "call-p-min100", "diversity-indices.tsv"),
                None,
                ["Full"],
                # high p
                100,
                # high r
                10,
                10,
                5,
                td,
                mock_log,
            )
        assert str(ei.value) == (
            "--high-p = 100 must be larger than the minimum p used in the BCF "
            "(100)."
        )
        assert len(os.listdir(td)) == 0


def test_run_estimate_selected_decoy_not_in_fasta():
    with tempfile.TemporaryDirectory() as td:

        # Case 1: the auto-selected decoy contig isn't in the FASTA
        # (note that we only check this one contig from the div indices; we
        # don't make sure the other contigs in the div indices file are in the
        # FASTA or BCF.)
        with pytest.raises(ParameterError) as ei:
            fu.run_estimate(
                FASTA,
                BAM,
                BCF,
                StringIO(
                    f"Contig\tAverageCoverage\tLength\t{DI_PREF}\t{DI_PREF}\n"
                    "edge_1\t35.2\t100\t0.5\t0.2\n"
                    "edge_2\t35.2\t100\t0.8\t0.2\n"
                ),
                None,
                ["Full"],
                # high p
                500,
                # high r
                2,
                10,
                5,
                td,
                mock_log,
            )
        assert str(ei.value) == (
            f"Selected decoy contig edge_1 is not present in {FASTA}."
        )
        assert len(os.listdir(td)) == 0

        # Case 2: the pre-selected decoy contig isn't in the FASTA
        with pytest.raises(ParameterError) as ei:
            fu.run_estimate(
                FASTA,
                BAM,
                BCF,
                None,
                "c4",
                ["Full"],
                # high p
                500,
                # high r
                2,
                10,
                5,
                td,
                mock_log,
            )
        assert str(ei.value) == (
            f"Selected decoy contig c4 is not present in {FASTA}."
        )
        assert len(os.listdir(td)) == 0


def test_parse_sco_good():
    df = fu.parse_sco(
        os.path.join("strainflye", "tests", "inputs", "edge_6104.sco")
    )
    assert len(df.index) == 1297
    # no need to go crazy testing this i think, let's just check a few genes
    pd.testing.assert_series_equal(
        df.loc[1],
        pd.Series(
            {"LeftEnd": 266, "RightEnd": 712, "Length": 447, "Strand": "-"},
            name=1,
        ),
    )
    pd.testing.assert_series_equal(
        df.loc[1217],
        pd.Series(
            {
                "LeftEnd": 1208927,
                "RightEnd": 1210075,
                "Length": 1149,
                "Strand": "+",
            },
            name=1217,
        ),
    )


def test_parse_sco_weird_line_prefix(tmp_path):
    # Yo apparently you can just have pytest take care of making a temporary
    # directory for you??? see https://stackoverflow.com/a/4205449 and
    # https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html
    #
    # i've been doing it manually like a chump LOLLLLLLLL
    fp = os.path.join(tmp_path, "test.sco")
    with open(fp, "w") as fh:
        fh.write(">1_266_712_-\n" "2_713_715_+\n")
    with pytest.raises(WeirdError) as ei:
        fu.parse_sco(fp)
    assert str(ei.value) == (
        'Unrecognized line prefix in SCO: line = "2_713_715_+"'
    )


def test_parse_sco_empty_lines_ok(tmp_path):
    fp = os.path.join(tmp_path, "test.sco")
    with open(fp, "w") as fh:
        fh.write(">1_266_712_-\n" "    \n" ">2_713_715_+\n" "  \t \n")
    obs_df = fu.parse_sco(fp)
    exp_df = pd.DataFrame(
        {
            "LeftEnd": [266, 713],
            "RightEnd": [712, 715],
            "Length": [447, 3],
            "Strand": ["-", "+"],
        },
        index=pd.Index([1, 2]),
    )
    pd.testing.assert_frame_equal(obs_df, exp_df)


def test_parse_sco_non3div_length(tmp_path):
    fp = os.path.join(tmp_path, "test.sco")
    with open(fp, "w") as fh:
        fh.write(">1_266_711_-\n")
    with pytest.raises(WeirdError) as ei:
        fu.parse_sco(fp)
    assert str(ei.value) == (
        "Gene 1 in SCO is 446 bp long; lengths must be divisible by 3."
    )


def test_parse_sco_left_not_lt_right(tmp_path):
    fp = os.path.join(tmp_path, "test.sco")
    with open(fp, "w") as fh:
        fh.write(">1_711_266_-\n")
    with pytest.raises(WeirdError) as ei:
        fu.parse_sco(fp)
    assert str(ei.value) == (
        "Gene 1 in SCO has left end of 711, which is not < the right end of "
        "266."
    )


def test_parse_sco_no_genes(tmp_path):
    fp = os.path.join(tmp_path, "test.sco")
    with open(fp, "w") as fh:
        fh.write("\n")
    with pytest.raises(WeirdError) as ei:
        fu.parse_sco(fp)
    assert str(ei.value) == "No genes described in SCO."


def test_parse_sco_bad_strand(tmp_path):
    fp = os.path.join(tmp_path, "test.sco")
    with open(fp, "w") as fh:
        fh.write(">1_266_711_$\n")
    with pytest.raises(WeirdError) as ei:
        fu.parse_sco(fp)
    assert str(ei.value) == "Unrecognized strand in SCO: $"


def test_get_next_cp():
    gene_namedtuple_type = collections.namedtuple("Gene", ["Index", "Strand"])

    test_gene = gene_namedtuple_type(100, "+")

    assert fu.get_next_cp(1, test_gene) == 2
    assert fu.get_next_cp(2, test_gene) == 3
    assert fu.get_next_cp(3, test_gene) == 1

    assert fu.get_next_cp(3, test_gene, ascending=False) == 2
    assert fu.get_next_cp(2, test_gene, ascending=False) == 1
    assert fu.get_next_cp(1, test_gene, ascending=False) == 3

    for bad_cp in (-3, -2, -1, -0.5, 0, 0.5, 1.5, 100, "very real number"):
        for strand in ("+", "-"):
            test_gene = gene_namedtuple_type(1234, strand)
            if strand == "+":
                asc = True
            else:
                asc = False
            with pytest.raises(WeirdError) as ei:
                fu.get_next_cp(bad_cp, test_gene, ascending=asc)
            assert str(ei.value) == (
                "Codon position got out of whack: gene 1,234, strand "
                f"{strand}, CP {bad_cp}"
            )


def test_complain_about_cps():
    with pytest.raises(WeirdError) as ei:
        fu.complain_about_cps(1234, "+", 4)
    assert str(ei.value) == (
        "Codon position got out of whack: gene 1,234, strand +, CP 4"
    )


def test_get_single_gene_and_cp2_positions_overlapping_genes():
    # looks like
    #
    #   *     *     *
    # ----- ----- -----
    # 1 2 3 4 5 6 7 8 9
    #           6 7 8 9 A B
    #           ----- -----
    #             *     *
    #
    # ... with codons marked with a -----, CP2 positions marked with a *,
    # and 10 and 11 replaced with A and B for the sake of keeping widths
    # consistent in my insane ascii art that no one will read.
    #
    # ANYWAY 8 and 7 should get ignored since they're in multiple genes
    df = pd.DataFrame(
        {
            "LeftEnd": [1, 6],
            "RightEnd": [9, 11],
            "Length": [9, 6],
            "Strand": ["+", "-"],
        },
        index=pd.Index([1234, 56789]),
    )
    assert fu.get_single_gene_and_cp2_positions(df) == (
        {1, 2, 3, 4, 5, 10, 11},
        {2, 5, 10},
    )


def test_get_single_gene_and_cp2_positions_one_gene():
    #   *     *     *
    # ----- ----- -----
    # 1 2 3 4 5 6 7 8 9
    df = pd.DataFrame(
        {
            "LeftEnd": [1],
            "RightEnd": [9],
            "Length": [9],
            "Strand": ["-"],
        },
        index=pd.Index([1]),
    )
    assert fu.get_single_gene_and_cp2_positions(df) == (
        {1, 2, 3, 4, 5, 6, 7, 8, 9},
        {2, 5, 8},
    )


def test_get_single_gene_and_cp2_positions_no_single_gene_pos():
    #   *     *     *
    # ----- ----- -----
    # 1 2 3 4 5 6 7 8 9
    # 1 2 3 4 5 6 7 8 9
    # ----- ----- -----
    #   *     *     *
    df = pd.DataFrame(
        {
            "LeftEnd": [1, 1],
            "RightEnd": [9, 9],
            "Length": [9, 9],
            "Strand": ["-", "+"],
        },
        index=pd.Index([1, 2]),
    )
    with pytest.raises(ParameterError) as ei:
        fu.get_single_gene_and_cp2_positions(df)

    assert str(ei.value) == (
        "No single-gene positions exist (given the predicted genes)."
    )


def test_get_single_gene_and_cp2_positions_no_single_gene_cp2_pos():
    # (No single-gene positions implies no single-gene CP2 positions, but
    # no single-gene CP2 positions does not imply no single-gene positions)
    #   *     *     *
    # ----- ----- -----
    # 1 2 3 4 5 6 7 8 9
    #   2 3 4 5 6 7 8 9 A
    #   ----- ----- -----
    #     *     *     *
    df = pd.DataFrame(
        {
            "LeftEnd": [1, 2],
            "RightEnd": [9, 10],
            "Length": [9, 9],
            "Strand": ["-", "+"],
        },
        index=pd.Index([1, 2]),
    )
    with pytest.raises(ParameterError) as ei:
        fu.get_single_gene_and_cp2_positions(df)

    assert str(ei.value) == (
        "No single-gene CP2 positions exist (given the predicted genes)."
    )


def test_compute_cp2_decoy_contig_mut_rates():
    # bit of a silly example -- edge_1 has "length" 500, as given in the BCF
    # text below, but we only define genes on the first few bp in this edge to
    # make testing easier.
    with write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220905\n"
        '##source="strainFlye v0.0.1: r-mutation calling (--min-r = 5)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=500>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minr_5, Description="min r threshold">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t2\t.\tT\tC\t.\t.\tMDP=10000;AAD=8\n"
        "edge_1\t7\t.\tA\tT\t.\t.\tMDP=10000;AAD=5\n"
        "edge_1\t8\t.\tA\tT\t.\t.\tMDP=10000;AAD=50\n"
        "edge_1\t10\t.\tT\tC\t.\t.\tMDP=10000;AAD=10\n"
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        genes_df = pd.DataFrame(
            {
                "LeftEnd": [1, 6],
                "RightEnd": [9, 11],
                "Length": [9, 6],
                "Strand": ["-", "+"],
            },
            index=pd.Index([1, 2]),
        )
        sg_pos, sg_cp2_pos = fu.get_single_gene_and_cp2_positions(genes_df)
        mut_rates = fu.compute_any_mutation_decoy_contig_mut_rates(
            bcf_obj,
            thresh_type,
            range(5, 13),
            "edge_1",
            pos_to_consider=sg_cp2_pos,
        )
        # The "denominator" is 3 * length, which in this case is the number of
        # CP2 single-gene positions considered (here, there are three such
        # positions: 2, 5, and 10)
        assert mut_rates == [
            # Two CP2 single-gene r-mutations for 5 <= r <= 8 (2 and 10)
            # r = 5
            2 / 9,
            # r = 6
            2 / 9,
            # r = 7
            2 / 9,
            # r = 8
            2 / 9,
            # Now, just one valid mutation left for r = 9 and r = 10
            1 / 9,
            1 / 9,
            # "Nothing beside remains" -- Ozymandias, from the hit blockbuster
            # movie Ozymandias 2: It's Ozymandin' Time
            0,
            0,
        ]


def test_compute_tv_decoy_contig_mut_rates():
    with write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220905\n"
        '##source="strainFlye v0.0.1: r-mutation calling (--min-r = 5)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=500>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minr_5, Description="min r threshold">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t1\t.\tA\tC\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t2\t.\tA\tG\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t3\t.\tA\tT\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t4\t.\tC\tA\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t5\t.\tC\tG\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t6\t.\tC\tT\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t7\t.\tG\tA\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t8\t.\tG\tC\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t9\t.\tG\tT\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t10\t.\tT\tA\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t11\t.\tT\tC\t.\t.\tMDP=10000;AAD=10\n"
        "edge_1\t12\t.\tT\tG\t.\t.\tMDP=10000;AAD=10\n"
    ) as fh:
        bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fh.name)
        # Try for all positions
        mut_rates = fu.compute_specific_mutation_decoy_contig_mut_rates(
            bcf_obj,
            thresh_type,
            range(5, 13),
            "edge_1",
            set(range(1, 501)),
            ["Tv"],
            skbio.DNA("A" * 500),
            # genes DF (not needed)
            None,
            # codon2cp2mts (not needed)
            None,
        )
        assert mut_rates == ([8 / 1000] * 6) + [0, 0]

        # Now, let's say certain positions are "unreasonable" -- this should
        # impact pos_to_consider. We'll label positions 1, 6, and 7 as
        # unreasonable: 1 is a transversion, and 6 and 7 are transitions.
        mut_rates = fu.compute_specific_mutation_decoy_contig_mut_rates(
            bcf_obj,
            thresh_type,
            range(5, 13),
            "edge_1",
            set(range(1, 501)) - {1, 6, 7},
            ["Tv"],
            skbio.DNA("A" * 500),
            # genes DF (not needed)
            None,
            # codon2cp2mts (not needed)
            None,
        )
        assert mut_rates == ([7 / 994] * 6) + [0, 0]


def test_get_mutation_types_good():
    assert fu.get_mutation_types(["Tv"]) == (True, False, False)
    assert fu.get_mutation_types(["Nonsyn"]) == (False, True, False)
    assert fu.get_mutation_types(["Nonsense"]) == (False, False, True)
    assert fu.get_mutation_types(["Tv", "Nonsyn"]) == (True, True, False)
    assert fu.get_mutation_types(["Tv", "Nonsense"]) == (True, False, True)


def test_get_mutation_types_bad():
    with pytest.raises(ParameterError) as ei:
        fu.get_mutation_types([])
    assert str(ei.value) == "mutation_types must be nonempty."

    with pytest.raises(ParameterError) as ei:
        fu.get_mutation_types(["Tv", "hamburgor"])
    assert str(ei.value) == "Unrecognized mutation type: hamburgor"

    # not case sensitive because these values should be generated internally
    # so no need to do that
    with pytest.raises(ParameterError) as ei:
        fu.get_mutation_types(["Nonsense", "tv"])
    assert str(ei.value) == "Unrecognized mutation type: tv"

    with pytest.raises(ParameterError) as ei:
        fu.get_mutation_types(["Tv", "Nonsyn", "Nonsense"])
    assert str(ei.value) == (
        "No need to specify both Nonsyn and Nonsense: just say Nonsense."
    )

    with pytest.raises(ParameterError) as ei:
        fu.get_mutation_types(["Nonsyn", "Nonsense"])
    assert str(ei.value) == (
        "No need to specify both Nonsyn and Nonsense: just say Nonsense."
    )


def test_get_mutation_types_for_cp():
    assert fu.get_mutation_types_for_cp("AAA", 3, "G") == (True, True)
    assert fu.get_mutation_types_for_cp("AAA", 3, "C") == (False, True)
    assert fu.get_mutation_types_for_cp("AAA", 1, "T") == (False, False)

    assert fu.get_mutation_types_for_cp("TAA", 3, "C") == (False, None)
    assert fu.get_mutation_types_for_cp("TAA", 3, "G") == (True, None)

    with pytest.raises(WeirdError) as ei:
        fu.get_mutation_types_for_cp("AAA", 1, "A")
    assert str(ei.value) == "In codon AAA, trying to mutate CP 1 into itself?"


def test_is_transversion_good():
    assert fu.is_transversion("A", "C")
    assert not fu.is_transversion("A", "G")
    assert fu.is_transversion("A", "T")

    assert fu.is_transversion("C", "A")
    assert fu.is_transversion("C", "G")
    assert not fu.is_transversion("C", "T")

    assert not fu.is_transversion("G", "A")
    assert fu.is_transversion("G", "C")
    assert fu.is_transversion("G", "T")

    assert fu.is_transversion("T", "A")
    assert not fu.is_transversion("T", "C")
    assert fu.is_transversion("T", "G")


def test_is_transversion_same():
    for nt in "ACGT":
        with pytest.raises(WeirdError) as ei:
            fu.is_transversion(nt, nt)

        assert str(ei.value) == (
            f"is_transversion() called with nt1 == nt2 == {nt}"
        )


def test_is_transversion_non_nucleotide():
    for nt in "ACGT":
        with pytest.raises(WeirdError) as ei:
            fu.is_transversion(nt, skbio.DNA("C"))

        assert str(ei.value) == (
            "is_transversion() parameters are not both str nucleotides. "
            f"nt1 == {nt}, nt2 == C. Check types?"
        )

        with pytest.raises(WeirdError) as ei:
            fu.is_transversion(nt, "N")

        assert str(ei.value) == (
            "is_transversion() parameters are not both str nucleotides. "
            f"nt1 == {nt}, nt2 == N. Check types?"
        )

        with pytest.raises(WeirdError) as ei:
            fu.is_transversion(nt, "AC")

        assert str(ei.value) == (
            "is_transversion() parameters are not both str nucleotides. "
            f"nt1 == {nt}, nt2 == AC. Check types?"
        )


def test_CodonPositionMutationCounts_good():
    mc = fu.CodonPositionMutationCounts("ACG", 3)

    assert mc.si == 3
    assert mc.si_tv == 2

    assert mc.ni == 0
    assert mc.ni_tv == 0

    assert mc.nnsi == 3
    assert mc.nnsi_tv == 2

    assert mc.nsi == 0
    assert mc.nsi_tv == 0

    mc = fu.CodonPositionMutationCounts("TGA", 2)

    assert mc.si == 1
    assert mc.si_tv == 0

    assert mc.ni == 2
    assert mc.ni_tv == 2

    assert mc.nnsi is None
    assert mc.nnsi_tv is None

    assert mc.nsi is None
    assert mc.nsi_tv is None


def test_CodonPositionMutationCounts_bad():
    with pytest.raises(WeirdError) as ei:
        mc = fu.CodonPositionMutationCounts("ACS", 3)
    assert str(ei.value) == "Codon should only contain {A, C, G, T}"

    with pytest.raises(WeirdError) as ei:
        mc = fu.CodonPositionMutationCounts("ACGT", 3)
    assert str(ei.value) == "Codon must be exactly 3 nt long"

    for bad_cp in (-100, -2, -1, 0, 0.5, 1.5, 2.5, 4, 100):
        with pytest.raises(WeirdError) as ei:
            mc = fu.CodonPositionMutationCounts("ACG", 4)
        assert str(ei.value) == "CP must be one of 1, 2, or 3"
