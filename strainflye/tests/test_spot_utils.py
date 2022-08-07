import os
import tempfile
import pytest
import strainflye.spot_utils as su
from strainflye.errors import ParameterError


def test_hotspot_no_mins():
    ti_dir = os.path.join("strainflye", "tests", "inputs", "small")
    with tempfile.NamedTemporaryFile() as out_fh:
        with pytest.raises(ParameterError) as ei:
            su.find_hotspot_features(
                os.path.join(ti_dir, "call-r-min3-di12345", "naive-calls.bcf"),
                os.path.join(ti_dir, "features.gff"),
                None,
                None,
                out_fh.name,
            )

        assert str(ei.value) == (
            "At least one of (--min-num-mutations, --min-perc-mutations) must "
            "be specified."
        )
