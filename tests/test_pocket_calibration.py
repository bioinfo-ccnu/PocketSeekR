"""Current quality/novelty selection and hard exclusions."""

import pytest

from pocketseekr.pocket_calibration import SelectionCalibration, select_calibrated


def pocket(name, x, score, atoms):
    return {
        "candidate_id": name,
        "center": [x, 0, 0],
        "score": score,
        "atom_indices": atoms,
    }


def test_zero_exponent_preserves_quality_order_but_not_duplicates():
    pool = [
        pocket("a", 0, 10, [0, 1, 2, 3]),
        pocket("b", 6, 9, [0, 1, 2, 4]),
        pocket("c", 12, 8, [5, 6, 7, 8]),
        pocket("duplicate", 18, 7, [0, 1, 2, 3]),
    ]
    assert [
        p["candidate_id"]
        for p in select_calibrated(pool, SelectionCalibration(overlap_exponent=0))
    ] == ["a", "b", "c"]
    assert [
        p["candidate_id"] for p in select_calibrated(pool, SelectionCalibration())
    ] == ["a", "c", "b"]


def test_truncation_applies_to_geometry_order_before_scoring_and_never_refills():
    pool = [pocket(f"r{i}", 0, i + 1, [0, 1]) for i in range(6)]
    selected = select_calibrated(pool, SelectionCalibration(region_budget=5))
    assert [p["candidate_id"] for p in selected] == ["r4"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"region_budget": 4},
        {"region_budget": 5.5},
        {"overlap_exponent": -1},
        {"overlap_exponent": float("nan")},
        {"overlap_exponent": float("inf")},
    ],
)
def test_invalid_calibration_parameters_raise(kwargs):
    with pytest.raises(ValueError):
        SelectionCalibration(**kwargs)
