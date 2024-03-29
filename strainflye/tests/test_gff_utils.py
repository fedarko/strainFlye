# Note that a fair amount of stuff in gff_utils is tested indirectly through
# test_spot_utils, so these tests aren't as woefully inadequate as they might
# seem at first. (The reason for this is that gff_utils' code was abstracted
# from code that was initially unique to spot_utils.)


import skbio
import pytest
import strainflye.gff_utils as gu
from io import StringIO as sio
from strainflye import config
from strainflye.errors import ParameterError, WeirdError
from strainflye.tests.utils_for_testing import mock_log


def test_validate_basic_good():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        sfi = set()
        fid, frange = gu.validate_basic(feature, "c1", 20, sfi, mock_log)
        assert fid == "hi"
        assert frange == range(4, 19)
        assert sfi == set(["hi"])
        break
    # this check mandates that we pass through the body of the for loop at
    # least once. it's a safeguard against scikit-bio's GFF3 parser suddenly
    # malfunctioning, causing this test to trivially "pass" due to cim_tuples
    # being incorrectly empty.
    if no_loops:
        raise WeirdError


def test_validate_basic_good_zero_indexed():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        sfi = set(["bye"])
        fid, frange = gu.validate_basic(
            feature, "c1", 20, sfi, mock_log, zero_indexed_range=False
        )
        assert fid == "hi"
        assert frange == range(5, 20)
        assert sfi == set(["hi", "bye"])
        break
    if no_loops:
        raise WeirdError


def test_validate_basic_multi_bound_feature():
    silly_feature = skbio.metadata.Interval(
        interval_metadata=skbio.metadata.IntervalMetadata(100),
        bounds=[(5, 10), (11, 20)],
    )
    with pytest.raises(ParameterError) as ei:
        gu.validate_basic(silly_feature, "c1", 100, set(), mock_log)
    assert str(ei.value) == (
        "A feature in the GFF3 file on contig c1 exists without exactly "
        f"one set of bounds: {str(silly_feature)}"
    )


def test_validate_if_cds_good():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        is_cds, strand = gu.validate_if_cds(feature, "c1", mock_log)
        assert is_cds
        assert strand == "+"
        break
    if no_loops:
        raise WeirdError


def test_validate_if_cds_noncds(capsys):
    gff = "##gff-version 3\nc1	marcus	gene	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        is_cds, strand = gu.validate_if_cds(feature, "c1", mock_log)
        assert not is_cds
        assert strand is None
        assert capsys.readouterr().out == (
            "MockLog: Feature hi on contig c1 has a type that is not in "
            f"{config.CDS_TYPES}; ignoring this feature.\n"
        )
        break
    if no_loops:
        raise WeirdError


def test_validate_if_cds_badphase():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	1	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        with pytest.raises(ParameterError) as ei:
            gu.validate_if_cds(feature, "c1", mock_log)
        assert str(ei.value) == (
            "Feature hi on contig c1 has a phase of 1. This command "
            "does not support features with non-zero phases, at least for now."
        )
        break
    if no_loops:
        raise WeirdError


def test_validate_if_cds_nophase():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	.	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        with pytest.raises(ParameterError) as ei:
            gu.validate_if_cds(feature, "c1", mock_log)
        assert str(ei.value) == (
            "Feature hi on contig c1 does not have a phase "
            "attribute; this is required (more specifically, required to be "
            "0) here for all CDS features."
        )
    if no_loops:
        raise WeirdError


def test_validate_if_cds_badstrand():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	.	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        with pytest.raises(ParameterError) as ei:
            gu.validate_if_cds(feature, "c1", mock_log)
        assert str(ei.value) == (
            'Feature hi on contig c1 has a strand of ".". We require here '
            'that all CDS features have a strand of "+" or "-".'
        )
    if no_loops:
        raise WeirdError


def test_validate_if_cds_badlength():
    gff = "##gff-version 3\nc1	marcus	cds	5	20	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        feature = next(im.query(metadata={}))
        with pytest.raises(ParameterError) as ei:
            gu.validate_if_cds(feature, "c1", mock_log)
        assert str(ei.value) == (
            "Feature hi on contig c1 has a length of 16 bp. We require here "
            "that all CDS features have lengths divisible by 3."
        )
    if no_loops:
        raise WeirdError
