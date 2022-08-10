import io
import os
import tempfile
import pytest
import strainflye.spot_utils as su
from strainflye.errors import ParameterError
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
            su.run_hotspot_detection(
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
            su.run_hotspot_detection(
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
            su.run_hotspot_detection(
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


def test_hotspot_contig_in_gff3_but_not_bcf():
    # Here is our "table" of what's ok and what isn't here:
    # (I put too much effort into formatting this)
    #
    # In BCF? | In GFF3? | Result
    # --------|----------|----------------------------------------------
    #   Yes   |   Yes    | OK -- this contig could have features defined
    # --------|----------|----------------------------------------------
    #   Yes   |   No     | OK -- this contig just doesn't have any features
    # --------|----------|----------------------------------------------
    #         |          | strainFlye should throw an error! Something's wrong.
    #   No    |   Yes    | (Or we could just ignore this contig, I guess, but
    #         |          | whatever.)
    # --------|----------|----------------------------------------------
    #   No    |   No     | Literally doesn't matter LOL
    # --------|----------|----------------------------------------------
    #
    # ... Anyway, we check row 3 (the throw-an-error case) in this test.

    gff_text = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=first_feature_that_i_made_up
c1	marcus	gene	1	5	.	+	0	ID=split_feature
c2	marcus	polyA_site	1	6	50	+	.	ID=another_thing
c2	marcus	gene	6	7	100	-	.	ID=worlds_shittiest_gene
c9	marcus	gene	6	7	100	-	.	ID=sus_feature
c3	marcus	exon	9	9	.	+	.	ID=single_nt_feature"""
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_detection(
                TEST_BCF_PATH,
                io.StringIO(gff_text),
                1,
                None,
                out_fh.name,
                mock_log,
            )
        assert str(ei.value) == (
            "The GFF3 file describes feature(s) located on contig c9, but "
            "this contig is not described in the BCF file."
        )
