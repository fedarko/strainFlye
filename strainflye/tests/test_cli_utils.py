from strainflye import __version__
from strainflye.cli_utils import fancystart, b2y, proglog, get_verboselog
from .utils_for_testing import mock_log


def test_fancylogging(capsys):
    # Uses pytest's capsys functionality to capture output -- see
    # https://docs.pytest.org/en/6.2.x/capture.html
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancylog = fancystart(*test_params, prefix="PREFIX ", version=False)
    exp_out = (
        "PREFIX strainFlye testing @ 0.00s: Starting...\n"
        "Input in1: Input #1\nInput in2: Input #2\nOutput out1: Output #1\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out

    fancylog("howdy", prefix="")
    captured = capsys.readouterr()
    assert captured.out == "strainFlye testing @ 0.00s: howdy\n"


def test_fancylogging_extrainfo(capsys):
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2")),
        (("out1", "Output #1"),),
    )
    fancylog = fancystart(
        *test_params,
        prefix="PREFIX ",
        extra_info=("Verbose?: Yes", "Sus?: Yeet"),
        version=False,
    )
    exp_out = (
        "PREFIX strainFlye testing @ 0.00s: Starting...\n"
        "Input in1: Input #1\nInput in2: Input #2\n"
        "Verbose?: Yes\nSus?: Yeet\n"
        "Output out1: Output #1\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out

    fancylog("howdy", prefix="")
    captured = capsys.readouterr()
    assert captured.out == "strainFlye testing @ 0.00s: howdy\n"


def test_fancylogging_version(capsys):
    test_params = (
        "strainFlye testing",
        (("in1", "Input #1"), ("in2", "Input #2"), ("in3", "Input #3")),
        (("out1", "Output #1"), ("out2", "Output #2")),
    )
    fancylog = fancystart(*test_params, prefix="PREFIX ")
    exp_out = (
        f'Using strainFlye version "{__version__}".\n'
        "PREFIX strainFlye testing @ 0.00s: Starting...\n"
        "Input in1: Input #1\n"
        "Input in2: Input #2\n"
        "Input in3: Input #3\n"
        "Output out1: Output #1\n"
        "Output out2: Output #2\n"
    )
    captured = capsys.readouterr()
    assert captured.out == exp_out
    fancylog("howdy", prefix="")
    captured = capsys.readouterr()
    assert captured.out == "strainFlye testing @ 0.00s: howdy\n"


def test_b2y():
    assert b2y(True) == "Yes"
    assert b2y(1) == "Yes"
    assert b2y(False) == "No"
    assert b2y(None) == "No"
    assert b2y(0) == "No"


def test_proglog(capsys):
    proglog("c1", 5, 10000, mock_log)
    captured = capsys.readouterr()
    assert (
        captured.out == "MockLog: On contig c1 (5 / 10,000 contigs = 0.05%).\n"
    )

    proglog("c1", 5, 10000, mock_log, contig_len=1234567890)
    captured = capsys.readouterr()
    assert captured.out == (
        "MockLog: On contig c1 (1,234,567,890 bp) (5 / 10,000 contigs = "
        "0.05%).\n"
    )

    proglog("c1", 5, 10000, mock_log, contig_len=1234567890, prefix="LOL")
    captured = capsys.readouterr()
    assert captured.out == (
        "MockLog: LOLcontig c1 (1,234,567,890 bp) (5 / 10,000 contigs = "
        "0.05%).\n"
    )


def test_get_verboselog(capsys):
    vl = get_verboselog(mock_log, False)
    vl("HI")
    vl("I SHOULDN'T BE DISPLAYED", prefix="butts lol")
    assert capsys.readouterr().out == ""
    vl = get_verboselog(mock_log, True)
    vl("OK SO LIKE")
    vl("NOW I SHOULD SHOW UP", prefix="")
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: OK SO LIKE\n" "MockLog: NOW I SHOULD SHOW UP\n"
    )
