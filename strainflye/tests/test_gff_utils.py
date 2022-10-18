import skbio
import pytest
import strainflye.gff_utils as gu
from io import StringIO as sio
from strainflye.errors import ParameterError
from strainflye.tests.utils_for_testing import mock_log


def test_validate_basic_good():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    for contig, im in cim_tuples:
        feature = next(im.query(metadata={}))
        sfi = set()
        fid, frange = gu.validate_basic(feature, "c1", 20, sfi, mock_log)
        assert fid == "hi"
        assert frange == range(4, 19)
        assert sfi == set(["hi"])
        break


def test_validate_basic_good_zero_indexed():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    for contig, im in cim_tuples:
        feature = next(im.query(metadata={}))
        sfi = set(["bye"])
        fid, frange = gu.validate_basic(
            feature, "c1", 20, sfi, mock_log, zero_indexed_range=False
        )
        assert fid == "hi"
        assert frange == range(5, 20)
        assert sfi == set(["hi", "bye"])
        break


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
