import io
import os
import tempfile
import pytest
import strainflye.spot_utils as su
from strainflye.errors import ParameterError
from .utils_for_testing import mock_log


GFF_BASE = """##gff-version 3
c1  marcus  cds 5   19  .   +   0   ID=first_feature_that_i_made_up
c1  marcus  gene    1   5   .   +   0   ID=split_feature
c1  marcus  gene    12  20  .   +   0   ID=split_feature
c2  marcus  polyA_site  1   6   50  +   .   ID=another_thing
c2  marcus  gene    6   7   100 -   .   ID=worlds_shittiest_gene
c3  marcus  exon    9   9   .   +   .   ID=single_nt_feature"""


def test_hotspot_no_mins():
    ti_dir = os.path.join("strainflye", "tests", "inputs", "small")
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.run_hotspot_detection(
                os.path.join(ti_dir, "call-r-min3-di12345", "naive-calls.bcf"),
                io.StringIO(GFF_BASE),
                None,
                None,
                out_fh.name,
                mock_log,
            )

        assert str(ei.value) == (
            "At least one of (--min-num-mutations, --min-perc-mutations) must "
            "be specified."
        )
