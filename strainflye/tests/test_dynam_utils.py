from strainflye import dynam_utils as du


def test_skew():
    assert du.skew("ATATATATT") == 0
    assert du.skew("") == 0

    assert du.skew("CCCCC") == -1
    assert du.skew("GGGGG") == 1
    assert du.skew("GC") == 0
    assert du.skew("CG") == 0
    assert du.skew("GGGGC") == 3 / 5
    assert du.skew("CCCCG") == -3 / 5
    assert du.skew("CGCCC") == -3 / 5
