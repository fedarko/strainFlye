from strainflye.cli_utils import fancystart


def test_fancylogging_verbose(capsys):
    # Tests fancystart() and fancylog() with verbose on.
    # Uses pytest's capsys functionality to capture output -- see
    # https://docs.pytest.org/en/6.2.x/capture.html
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancylog = fancystart(*test_params, True, prefix="PREFIX ")
    exp_out = (
        "PREFIX strainFlye testing @ 0.00 sec: Starting...\n"
        "Input in1: Input #1\nInput in2: Input #2\nOutput out1: Output #1\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out

    fancylog("howdy", prefix="")
    captured = capsys.readouterr()
    assert captured.out == "strainFlye testing @ 0.00 sec: howdy\n"


def test_fancylogging_not_verbose(capsys):
    # Tests fancystart() and fancylog() with verbose off.
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancylog = fancystart(*test_params, False, prefix="PREFIX ")
    captured = capsys.readouterr()
    # Changing verbose to False means that nothing should be logged by us
    assert captured.out == ""
