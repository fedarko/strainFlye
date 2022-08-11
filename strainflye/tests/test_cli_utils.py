from strainflye.cli_utils import fancystart, b2y


def test_fancylogging(capsys):
    # Tests fancystart() and fancylog() with quiet off, as is the default.
    # Uses pytest's capsys functionality to capture output -- see
    # https://docs.pytest.org/en/6.2.x/capture.html
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancylog = fancystart(*test_params, prefix="PREFIX ")
    exp_out = (
        "PREFIX strainFlye testing @ 0.00 sec: Starting...\n"
        "Input in1: Input #1\nInput in2: Input #2\nOutput out1: Output #1\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out

    fancylog("howdy", prefix="")
    captured = capsys.readouterr()
    assert captured.out == "strainFlye testing @ 0.00 sec: howdy\n"


def test_fancylogging_quiet(capsys):
    # Tests fancystart() and fancylog() with quiet on.
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancystart(*test_params, quiet=True, prefix="PREFIX ")
    captured = capsys.readouterr()
    # Changing quiet to True means that nothing should be logged by us
    assert captured.out == ""


def test_fancylogging_extrainfo(capsys):
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancylog = fancystart(
        *test_params,
        prefix="PREFIX ",
        extra_info=("Verbose?: Yes", "Sus?: Yeet")
    )
    exp_out = (
        "PREFIX strainFlye testing @ 0.00 sec: Starting...\n"
        "Input in1: Input #1\nInput in2: Input #2\n"
        "Verbose?: Yes\nSus?: Yeet\n"
        "Output out1: Output #1\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out

    fancylog("howdy", prefix="")
    captured = capsys.readouterr()
    assert captured.out == "strainFlye testing @ 0.00 sec: howdy\n"


def test_b2y():
    assert b2y(True) == "Yes"
    assert b2y(1) == "Yes"
    assert b2y(False) == "No"
    assert b2y(None) == "No"
    assert b2y(0) == "No"
