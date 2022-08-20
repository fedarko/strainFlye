import strainflye.smooth_utils as su


def test_convert_to_runs():
    assert su.convert_to_runs([]) == []
    assert su.convert_to_runs([3]) == [(3, 3)]
    assert su.convert_to_runs([3, 4, 9, 10, 12, 22]) == [
        (3, 4),
        (9, 10),
        (12, 12),
        (22, 22),
    ]
    assert su.convert_to_runs(
        [1, 2, 3, 5, 6, 7, 10, 20, 31, 32, 33, 34, 35, 40]
    ) == [(1, 3), (5, 7), (10, 10), (20, 20), (31, 35), (40, 40)]
