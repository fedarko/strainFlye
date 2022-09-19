import io
import os
import tempfile
import pytest
import numpy as np
import pandas as pd
import strainflye.spot_utils as su
from pytest import approx
from strainflye.errors import ParameterError, WeirdError
from .utils_for_testing import mock_log


TEST_BCF_PATH = os.path.join(
    "strainflye",
    "tests",
    "inputs",
    "small",
    "call-r-min3-di12345",
    "naive-calls.bcf",
)


def test_hotspot_no_mins_specified():
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=2
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                None,
                None,
                out_fh.name,
                mock_log,
            )

        assert str(ei.value) == (
            "At least one of (--min-num-mutations, --min-perc-mutations) must "
            "be specified."
        )


def test_hotspot_one_feature_multi_rows():
    # TLDR, we don't support this type of "discontinuous" feature yet
    # so if we see it, complain
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c1	marcus	gene	12	20	.	+	0	ID=split_feature
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "The feature ID split_feature is used in multiple GFF3 rows for "
            "contig c1. Features of a contig must have unique IDs; this "
            'command does not support "discontinuous features" at the moment.'
        )


def test_hotspot_feature_with_no_id():
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	sussy=baka
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert (
            "A feature in the GFF3 file on contig c1 exists without a defined "
            "ID: "
        ) in str(ei.value)


def test_hotspot_all_contigs_in_gff3_but_not_bcf():
    # Here is our "table" of what's ok and what isn't here, where each row
    # represents the "status" of a contig:
    # (I put too much effort into formatting this)
    #
    # In BCF? | In GFF3? | Result
    # --------|----------|----------------------------------------------
    #   Yes   |   Yes    | OK -- this contig could have features defined
    # --------|----------|----------------------------------------------
    #   Yes   |   No     | OK -- this contig just doesn't have any features
    # --------|----------|----------------------------------------------
    #         |          | OK -- just skip this contig and move on.
    #   No    |   Yes    | (Previously, we raised an error, but that was
    #         |          | clunky and annoying.)
    # --------|----------|----------------------------------------------
    #   No    |   No     | Literally doesn't matter LOL
    # --------|----------|----------------------------------------------
    #
    # We raise an error if no features are described at all in the GFF3 file,
    # or if the features described in the GFF3 file are not located on any of
    # the BCF contigs. (So, if every feature in the GFF3 file falls into "row
    # 3" in the table above, then we'd throw an error.) We test this case here.

    gff_text = """##gff-version 3
c20	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c9	marcus	gene	6	7	100	-	.	ID=sus_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "None of the feature(s) described in the GFF3 file are located on "
            "contigs that are described in the BCF file."
        )


def test_hotspot_some_contigs_in_gff3_but_not_bcf(capsys):
    # Unlike the above test, just having a few GFF3 contigs not be in the BCF
    # is ok
    gff_text = """##gff-version 3
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c20	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c9	marcus	gene	6	7	100	-	.	ID=sus_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            1,
            None,
            out_fh.name,
            mock_log,
        )
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1"],
                "FeatureID": ["split_feature"],
                "FeatureStart_1IndexedInclusive": [1],
                "FeatureEnd_1IndexedInclusive": [5],
                "NumberMutatedPositions": [1],
                "PercentMutatedPositions": ["20.00%"],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Identified 1 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    assert capsys.readouterr().out == exp_out


def test_hotspot_empty_gff3():
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(""),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "The GFF3 file doesn't seem to describe any features."
        )


def test_hotspot_feature_starts_far_after_contig():
    gff_text = """##gff-version 3
c1	marcus	cds	50000	100000	.	+	0	ID=sus_feature
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "Feature sus_feature on contig c1 has a (1-indexed) start "
            "coordinate of 50,000, which is greater than the contig's "
            "length of 23."
        )


def test_hotspot_feature_starts_close_after_contig():
    gff_text = """##gff-version 3
c1	marcus	cds	5	10	.	+	0	ID=sus_feature
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	13	14	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "Feature worlds_shittiest_gene on contig c2 has a (1-indexed) "
            "start coordinate of 13, which is greater than the contig's "
            "length of 12."
        )


def test_hotspot_feature_ends_after_contig():
    gff_text = """##gff-version 3
c1	marcus	cds	5	10	.	+	0	ID=sus_feature
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	2	13	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "Feature worlds_shittiest_gene on contig c2 has a (1-indexed) "
            "end coordinate of 13, which is greater than the contig's "
            "length of 12. We do not support 'circular' features yet."
        )


def test_hotspot_feature_starts_before_contig():
    gff_text = """##gff-version 3
c1	marcus	cds	0	10	.	+	0	ID=impostor
c3	marcus	exon	8	8	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ValueError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        # skbio throws this error -- it knows that GFF3 files use one-indexing,
        # so a feature that has a coordinate of zero (or of -1, -2, ...) will
        # cause it to throw an error. Less work for us!
        assert str(ei.value) == (
            "Cannot set `bounds` ([(-1, 10)]) with coordinate smaller than "
            "lower bound (0)."
        )


def test_hotspot_feature_fractional_coord():
    gff_text = """##gff-version 3
c1	marcus	cds	1.5	10	.	+	0	ID=impostor
c3	marcus	exon	8	8	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ValueError) as ei:
            su.run_hotspot_feature_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        # Another error caused by skbio that we verify happens
        assert str(ei.value) == "invalid literal for int() with base 10: '1.5'"


def test_hotspot_good_min_num_1(capsys):
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=2
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	8	8	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        # A feature is a "hotspot" if it has at least 1 mutation
        # (translator's note: this is silly lol)
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            1,
            None,
            out_fh.name,
            mock_log,
        )
        # Notably, we don't use index_col=0 here -- because the first column
        # (contig) can be the same for multiple features, and the second column
        # (feature ID) can be the same across different contigs. So, we just
        # fall back on pandas' default of a RangeIndex (as of writing).
        # If this is a huge deal, we can modify the hotspot code to enforce
        # that feature IDs must be unique across ALL contigs in the file. But I
        # don't think anyone will really care, tbh.
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1", "c1", "c3"],
                "FeatureID": [
                    "first_feature_that_i_made_up",
                    "2",
                    "single_nt_feature",
                ],
                "FeatureStart_1IndexedInclusive": [5, 1, 8],
                "FeatureEnd_1IndexedInclusive": [19, 5, 8],
                "NumberMutatedPositions": [2, 1, 1],
                "PercentMutatedPositions": ["13.33%", "20.00%", "100.00%"],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df)

    # The test BCF file we use includes all naive r-mutations for r = 3.
    # We should thus see three hotspots (defined here as features with >= 1
    # mutation):
    # - c1, "first_feature_that_i_made_up" (has two mutations, at pos 11 & 13)
    # - c1, "2" (has one mutation, at pos 4)
    # - c3, "single_nt_feature" (has one mutation, at pos 8)
    #
    # ... All those coordinates use 1-indexing, fyi. See the comments in
    # test_run_small_dataset_r() in test_call_utils.py for details on the
    # mutations in this test dataset.
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Warning: feature single_nt_feature on contig c3 has equal "
        "start and end coordinates. We assume this refers to a feature of "
        "length 1 spanning this single position, rather than a feature of "
        "length 0.\n"
        "MockLog: Identified 3 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )

    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_hotspot_good_min_num_10_no_hotspots(capsys):
    # (so, no hotspots)
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=2
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	8	8	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            10,
            None,
            out_fh.name,
            mock_log,
        )
        # output file should just be the TSV header and no actual hotspots,
        # since no hotspots should have been detected
        with open(out_fh.name, "r") as of:
            text = of.read()
            assert text == (
                "Contig\tFeatureID\tFeatureStart_1IndexedInclusive\t"
                "FeatureEnd_1IndexedInclusive\tNumberMutatedPositions\t"
                "PercentMutatedPositions\n"
            )
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Warning: feature single_nt_feature on contig c3 has equal "
        "start and end coordinates. We assume this refers to a feature of "
        "length 1 spanning this single position, rather than a feature of "
        "length 0.\n"
        "MockLog: Identified 0 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_hotspot_good_min_perc_20(capsys):
    # (so, no hotspots)
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=2
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	8	8	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            None,
            20,
            out_fh.name,
            mock_log,
        )
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1", "c3"],
                "FeatureID": [
                    "2",
                    "single_nt_feature",
                ],
                "FeatureStart_1IndexedInclusive": [1, 8],
                "FeatureEnd_1IndexedInclusive": [5, 8],
                "NumberMutatedPositions": [1, 1],
                "PercentMutatedPositions": ["20.00%", "100.00%"],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Warning: feature single_nt_feature on contig c3 has equal "
        "start and end coordinates. We assume this refers to a feature of "
        "length 1 spanning this single position, rather than a feature of "
        "length 0.\n"
        "MockLog: Identified 2 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_hotspot_good_min_num_2_min_perc_20_no_hotspots(capsys):
    # (so, no hotspots, again)
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=2
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c3	marcus	exon	8	8	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            2,
            20,
            out_fh.name,
            mock_log,
        )
        # output file should just be the TSV header and no actual hotspots,
        # since no hotspots should have been detected
        with open(out_fh.name, "r") as of:
            text = of.read()
            assert text == (
                "Contig\tFeatureID\tFeatureStart_1IndexedInclusive\t"
                "FeatureEnd_1IndexedInclusive\tNumberMutatedPositions\t"
                "PercentMutatedPositions\n"
            )
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Warning: feature single_nt_feature on contig c3 has equal "
        "start and end coordinates. We assume this refers to a feature of "
        "length 1 spanning this single position, rather than a feature of "
        "length 0.\n"
        "MockLog: Identified 0 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_hotspot_good_same_fid_diff_contigs(capsys):
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=doppelganger
c3	marcus	exon	8	8	.	+	.	ID=doppelganger"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            1,
            None,
            out_fh.name,
            mock_log,
        )
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1", "c3"],
                "FeatureID": ["doppelganger", "doppelganger"],
                "FeatureStart_1IndexedInclusive": [5, 8],
                "FeatureEnd_1IndexedInclusive": [19, 8],
                "NumberMutatedPositions": [2, 1],
                "PercentMutatedPositions": ["13.33%", "100.00%"],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Warning: feature doppelganger on contig c3 has equal "
        "start and end coordinates. We assume this refers to a feature of "
        "length 1 spanning this single position, rather than a feature of "
        "length 0.\n"
        "MockLog: Identified 2 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_hotspot_good_feature_starts_and_ends_on_last_position_of_contig(
    capsys,
):
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=f1
c3	marcus	exon	16	16	.	+	.	ID=f2"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            1,
            None,
            out_fh.name,
            mock_log,
        )
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1"],
                "FeatureID": ["f1"],
                "FeatureStart_1IndexedInclusive": [5],
                "FeatureEnd_1IndexedInclusive": [19],
                "NumberMutatedPositions": [2],
                "PercentMutatedPositions": ["13.33%"],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Warning: feature f2 on contig c3 has equal "
        "start and end coordinates. We assume this refers to a feature of "
        "length 1 spanning this single position, rather than a feature of "
        "length 0.\n"
        "MockLog: Identified 1 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_hotspot_good_feature_ends_on_last_position_of_contig(capsys):
    # Basically, just decreasing the end coordinate of "worlds_shittiest_gene"
    # from 13 (see test_hotspot_feature_ends_after_contig() above) to a
    # 12 causes this to not fail.
    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=f1
c2	marcus	gene	2	12	100	-	.	ID=worlds_shittiest_gene"""
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_hotspot_feature_detection(
            TEST_BCF_PATH,
            io.StringIO(gff_text),
            1,
            None,
            out_fh.name,
            mock_log,
        )
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1"],
                "FeatureID": ["f1"],
                "FeatureStart_1IndexedInclusive": [5],
                "FeatureEnd_1IndexedInclusive": [19],
                "NumberMutatedPositions": [2],
                "PercentMutatedPositions": ["13.33%"],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through features in the GFF3 file and "
        "identifying hotspot features...\n"
        "MockLog: Identified 1 hotspot feature(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_coldspot_good_nocircular(capsys):
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_coldspot_gap_detection(
            TEST_BCF_PATH,
            5,
            False,
            out_fh.name,
            mock_log,
        )
        # For reference:
        # c1 has len 23 and mutations at (1-indexed) pos 4, 11, 13
        # c2 has len 12 and no mutations
        # c3 has len 16 and mutations at (1-indexed) pos 7, 8
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1", "c1", "c2", "c3", "c3"],
                "Start_1IndexedInclusive": [5, 14, 1, 1, 9],
                "End_1IndexedInclusive": [10, 23, 12, 6, 16],
                "Length": [6, 10, 12, 6, 8],
                "P_Value": [
                    np.nan,
                    approx(0.6392966, abs=1e-6),
                    np.nan,
                    np.nan,
                    approx(0.6872178, abs=1e-6),
                ],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df, check_dtype=False)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through contigs and identifying coldspot "
        "gaps...\n"
        "MockLog: Identified 5 coldspot gap(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_coldspot_good_circular(capsys):
    with tempfile.NamedTemporaryFile() as out_fh:
        su.run_coldspot_gap_detection(
            TEST_BCF_PATH,
            5,
            True,
            out_fh.name,
            mock_log,
        )
        # For reference:
        # c1 has len 23 and mutations at (1-indexed) pos 4, 11, 13
        # c2 has len 12 and no mutations
        # c3 has len 16 and mutations at (1-indexed) pos 7, 8
        obs_df = pd.read_csv(out_fh, sep="\t")
        exp_df = pd.DataFrame(
            {
                "Contig": ["c1", "c1", "c2", "c3"],
                "Start_1IndexedInclusive": [5, 14, 1, 9],
                "End_1IndexedInclusive": [10, 3, 12, 6],
                "Length": [6, 13, 12, 14],
                "P_Value": [
                    np.nan,
                    approx(0.3745209, abs=1e-6),
                    np.nan,
                    approx(0.19276259, abs=1e-6),
                ],
            }
        )
        pd.testing.assert_frame_equal(obs_df, exp_df, check_dtype=False)
    exp_out = (
        "PREFIX\nMockLog: Loading and checking the BCF file...\n"
        "MockLog: Looks good so far.\n"
        "PREFIX\nMockLog: Going through contigs and identifying coldspot "
        "gaps...\n"
        "MockLog: Identified 4 coldspot gap(s) across all 3 contigs in the "
        "BCF file.\n"
        "PREFIX\nMockLog: Writing out this information to a TSV file...\n"
        "MockLog: Done.\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out


def test_coldspot_zero_min_length():
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_coldspot_gap_detection(
                TEST_BCF_PATH,
                0,
                True,
                out_fh.name,
                mock_log,
            )

        assert str(ei.value) == (
            "Minimum coldspot gap length must be at least 1."
        )


def test_get_coldspot_gaps_len1():
    # If there's just one position, and it is a mutation, then there are no
    # gaps -- even of length 1.
    assert su.get_coldspot_gaps_in_contig([1], 1, 1, False) == []
    assert su.get_coldspot_gaps_in_contig([1], 1, 1, True) == []
    # ... or of length 2.
    assert su.get_coldspot_gaps_in_contig([1], 1, 2, False) == []
    assert su.get_coldspot_gaps_in_contig([1], 1, 2, True) == []

    # However, if this position isn't a mutation, it counts as a gap of len 1.
    assert su.get_coldspot_gaps_in_contig([], 1, 1, False) == [(1, 1, 1)]
    assert su.get_coldspot_gaps_in_contig([], 1, 1, True) == [(1, 1, 1)]


def test_get_coldspot_gaps_1mut_at_2len_end():
    # Contigs of length 2, where there is a single mutation at the left end...
    assert su.get_coldspot_gaps_in_contig([1], 2, 1, True) == [(2, 2, 1)]
    # ...or at the right end.
    assert su.get_coldspot_gaps_in_contig([2], 2, 1, True) == [(1, 1, 1)]

    # Increasing the minimum gap length to 2 causes us to not find anything.
    assert su.get_coldspot_gaps_in_contig([1], 2, 2, True) == []
    assert su.get_coldspot_gaps_in_contig([2], 2, 2, True) == []


def test_get_coldspot_gaps_weird_circular_3len():

    # Zero mutations:
    # NO NO NO
    assert su.get_coldspot_gaps_in_contig([], 3, 1, True) == [(1, 3, 3)]

    # One mutation:
    # MUT NO NO
    assert su.get_coldspot_gaps_in_contig([1], 3, 1, True) == [(2, 3, 2)]
    # NO MUT NO
    assert su.get_coldspot_gaps_in_contig([2], 3, 1, True) == [(3, 1, 2)]
    # NO NO MUT
    assert su.get_coldspot_gaps_in_contig([3], 3, 1, True) == [(1, 2, 2)]

    # Two mutations:
    # NO MUT MUT
    assert su.get_coldspot_gaps_in_contig([2, 3], 3, 1, True) == [(1, 1, 1)]
    # MUT NO MUT
    assert su.get_coldspot_gaps_in_contig([1, 3], 3, 1, True) == [(2, 2, 1)]
    # MUT MUT NO
    assert su.get_coldspot_gaps_in_contig([1, 2], 3, 1, True) == [(3, 3, 1)]

    # Three mutations:
    # MUT MUT MUT
    assert su.get_coldspot_gaps_in_contig([1, 2, 3], 3, 1, True) == []


def test_get_coldspot_gaps_no_mutations():
    # Cases where, despite there being no muts, the contig isn't long enough to
    # qualify as a gap
    for min_len in range(101, 201):
        assert su.get_coldspot_gaps_in_contig([], 100, min_len, True) == []

    # Cases where the contig is now long enough
    for min_len in range(1, 101):
        assert su.get_coldspot_gaps_in_contig([], 100, min_len, True) == [
            (1, 100, 100)
        ]


def test_get_coldspot_gaps_one_mutation_no_circular():
    # A contig has positions 1, 2, ..., 99, 100, with a mutation at pos 50.
    #
    # There are thus two gaps: [1, 49] (length 49) and [51, 100] (length 50).

    # Both gaps are long enough
    assert su.get_coldspot_gaps_in_contig([50], 100, 5, False) == [
        (1, 49, 49),
        (51, 100, 50),
    ]

    # Only one is long enough
    assert su.get_coldspot_gaps_in_contig([50], 100, 50, False) == [
        (51, 100, 50)
    ]

    # Neither is long enough
    assert su.get_coldspot_gaps_in_contig([50], 100, 51, False) == []


def test_get_coldspot_gaps_many_mutations():

    assert su.get_coldspot_gaps_in_contig([50, 55, 60], 100, 5, False) == [
        (1, 49, 49),
        (61, 100, 40),
    ]

    assert su.get_coldspot_gaps_in_contig([50, 55, 60], 100, 49, False) == [
        (1, 49, 49)
    ]

    assert su.get_coldspot_gaps_in_contig([50, 55, 60], 100, 5, True) == [
        (61, 49, 89)
    ]

    # positions are [1, 2, 3, ..., 98, 99, 100].
    # If each position is mutated, there are no coldspots.
    assert su.get_coldspot_gaps_in_contig(range(1, 101), 100, 1, True) == []

    # Try removing just the first mutated position -- creates a "gap" at this
    # one position
    assert su.get_coldspot_gaps_in_contig(range(2, 101), 100, 1, True) == [
        (1, 1, 1)
    ]

    # same deal with the last position
    assert su.get_coldspot_gaps_in_contig(range(1, 100), 100, 1, True) == [
        (100, 100, 1)
    ]

    # circular of last + first
    assert su.get_coldspot_gaps_in_contig(range(2, 100), 100, 1, True) == [
        (100, 1, 2)
    ]


def test_get_coldspot_gaps_literal_docs_example():
    # "better safe than sorry" -- not the guy who designed the titanic

    assert su.get_coldspot_gaps_in_contig([4, 6], 9, 1, True) == [
        (5, 5, 1),
        (7, 3, 6),
    ]

    assert su.get_coldspot_gaps_in_contig([4, 6], 9, 1, False) == [
        (5, 5, 1),
        (1, 3, 3),
        (7, 9, 3),
    ]

    # other stuff, just to be paranoid
    assert su.get_coldspot_gaps_in_contig([4, 6], 9, 2, False) == [
        (1, 3, 3),
        (7, 9, 3),
    ]


def test_longest_success_run_pvalue_bad():
    for m in range(7, 20):
        with pytest.raises(WeirdError) as ei:
            su.longest_success_run_pvalue(m, 7, 0.5)
        assert str(ei.value) == "n must be greater than m."

    for m in [-100, -2, -1, -0.5, 0, 0.5]:
        with pytest.raises(WeirdError) as ei:
            su.longest_success_run_pvalue(m, 7, 0.5)
        assert str(ei.value) == "m must be at least 1."

    for p in [-101, -100, -1, -0.5, 1.1, 100, 101]:
        with pytest.raises(WeirdError) as ei:
            su.longest_success_run_pvalue(7, 50, p)
        assert str(ei.value) == "p must be in the range [0, 1]."


def test_longest_success_run_pvalue_naus1982():
    # Values taken from the table on page 179 in Naus 1982:
    # https://www.tandfonline.com/doi/abs/10.1080/01621459.1982.10477783
    assert su.longest_success_run_pvalue(7, 50, 0.5) == approx(
        0.1653, abs=1e-4
    )
    assert su.longest_success_run_pvalue(10, 50, 0.5) == approx(
        0.0204, abs=1e-4
    )
    assert su.longest_success_run_pvalue(5, 50, 1 / 3) == approx(
        0.1214, abs=1e-4
    )
    assert su.longest_success_run_pvalue(4, 40, 1 / 3) == approx(
        0.2739, abs=1e-4
    )
    assert su.longest_success_run_pvalue(2, 16, 0.2) == approx(
        0.4107, abs=1e-4
    )


def test_get_coldspot_gap_pvalues_naus1982():
    # Using the first two values from the aforementioned table in Naus 1982.
    # We use num_muts == 25 in order to represent the mutation rate (p) being
    # 1/2, given n == 50.
    # (We could in theory use floating-point num_muts values, but that wouldn't
    # make sense -- and these cases are already tested above anyway.)
    assert su.get_coldspot_gap_pvalues(25, 50, [7]) == [
        approx(0.1653, abs=1e-4)
    ]
    assert su.get_coldspot_gap_pvalues(25, 50, [10]) == [
        approx(0.0204, abs=1e-4)
    ]


def test_get_coldspot_gap_pvalues_more():
    assert su.get_coldspot_gap_pvalues(5, 100, [17, 21, 19, 3]) == [
        "NA",
        approx(0.97316991876, abs=1e-6),
        "NA",
        "NA",
    ]
    assert su.get_coldspot_gap_pvalues(6, 100, [1, 1, 1, 1, 1]) == [
        approx(1),
        "NA",
        "NA",
        "NA",
        "NA",
    ]
    assert su.get_coldspot_gap_pvalues(98, 100, [2]) == [
        approx(0.0381085137, abs=1e-6)
    ]


def test_get_coldspot_gap_pvalues_zero_muts():
    assert su.get_coldspot_gap_pvalues(0, 100, [100]) == ["NA"]

    # test case where there's exactly 1 coldspot but the length is weird
    exp_err_msg = (
        "A contig with 0 mutations must have exactly 1 coldspot covering the "
        "entire contig."
    )
    for cl in list(range(-100, 100)) + list(range(101, 200)):
        with pytest.raises(WeirdError) as ei:
            su.get_coldspot_gap_pvalues(0, 100, [cl])
        assert str(ei.value) == exp_err_msg

    # test case where there's != 1 coldspot
    for coldspot_list in ([100, 1], [99, 1], [5, 5, 5], [], [50]):
        with pytest.raises(WeirdError) as ei:
            su.get_coldspot_gap_pvalues(0, 100, coldspot_list)
        assert str(ei.value) == exp_err_msg


def test_get_coldspot_gap_pvalues_all_muts():
    assert su.get_coldspot_gap_pvalues(100, 100, []) == []

    with pytest.raises(WeirdError) as ei:
        su.get_coldspot_gap_pvalues(100, 100, [1])
    assert str(ei.value) == (
        "There can't be any gaps if every position is mutated."
    )
