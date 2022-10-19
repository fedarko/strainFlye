import skbio
import pytest
import pandas as pd
import strainflye.matrix_utils as mu
from io import StringIO as sio
from collections import defaultdict
from strainflye import config
from strainflye.errors import ParameterError
from strainflye.tests.utils_for_testing import mock_log, mock_log_2


def test_get_contig_cds_info_good():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    for contig, im in cim_tuples:
        cds_df, fid2codon2alignedcodons = mu.get_contig_cds_info(
            im, contig, {"c1": 23}, mock_log, mock_log_2
        )
        pd.testing.assert_frame_equal(
            cds_df,
            pd.DataFrame(
                {"LeftEnd": [5], "RightEnd": [19], "Strand": "+"},
                index=pd.Index(["hi"]),
            ),
        )
        assert fid2codon2alignedcodons == {
            "hi": {
                5: defaultdict(int),
                8: defaultdict(int),
                11: defaultdict(int),
                14: defaultdict(int),
                17: defaultdict(int),
            },
        }
        break
