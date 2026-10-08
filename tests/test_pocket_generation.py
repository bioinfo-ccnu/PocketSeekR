"""Current region generation, geometry and chemistry invariants."""

from dataclasses import replace
from itertools import product

import numpy as np
import pytest
import torch

from pocketseekr.config import StageAConfig
from pocketseekr.pocket_directional import exposed_base_potential
from pocketseekr.pocket_generation import (
    GenerationCalibration,
    RegionalPocketSearch,
    select_generation,
)
from pocketseekr.utils.io import _extract_rna_atoms, parse_pdb_atoms


@pytest.fixture
def search(tmp_path):
    xyz = torch.tensor(
        [v for v in product([-1.0, 0, 1.0], repeat=3) if any(v)], dtype=torch.float64
    )
    xyz = 5 * xyz / torch.linalg.norm(xyz, dim=1)[:, None]
    source = tmp_path / "shell.pdb"
    source.write_text(
        "\n".join(
            [
                f"ATOM  {i:5d}  C8    A A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
                for i, (x, y, z) in enumerate(xyz.tolist(), 1)
            ]
            + ["END"]
        )
        + "\n"
    )
    atoms = _extract_rna_atoms(parse_pdb_atoms(source), source)
    stage = StageAConfig()
    return RegionalPocketSearch(atoms, torch.stack([a.position for a in atoms]), stage)


def test_spatial_chemical_screening_preserves_fields_and_site_order(search):
    for i in range(min(15, len(search.grid[1]))):
        p = search.grid[1][i]
        expected = exposed_base_potential(p, search.atoms, search.geometry, search.tree)
        assert search.point_chemistry(i) == (
            expected["face_support"],
            expected["polar_support"],
        )
    # Two visible rings with different centers, plus a distant excluded ring.
    from pocketseekr.pocket_physics import Ring
    from scipy.spatial import cKDTree

    q = search.grid[1][0]
    rings = [
        Ring(q - np.array([0.0, 0, 3.5]), np.array([0.0, 0, 1.0])),
        Ring(q - np.array([1.0, 0, 3.6]), np.array([0.0, 0, 1.0])),
        Ring(q - np.array([0.0, 0, 10.0]), np.array([0.0, 0, 1.0])),
    ]
    search.geometry = replace(search.geometry, rings=rings)
    search.site_trees[0] = (rings, cKDTree([r.center for r in rings]))
    search.point_cache.clear()
    expected = exposed_base_potential(q, search.atoms, search.geometry, search.tree)
    assert search.point_chemistry(0) == (
        expected["face_support"],
        expected["polar_support"],
    )


def test_explicit_growth_does_not_mutate_parent_protocol_or_shared_grid(search):
    from pocketseekr.pocket_regions import REGION_PROTOCOL

    old = dict(REGION_PROTOCOL)
    points = search.grid[1].copy()
    search.pool(GenerationCalibration(growth_radius_A=4, growth_min_fraction=0.3))
    assert REGION_PROTOCOL == old
    assert np.array_equal(search.grid[1], points)
    assert len(search.region_cache) == 1
    search.pool(GenerationCalibration())
    assert len(search.region_cache) == 2


def test_chemical_weight_changes_quality_without_changing_centers_or_crops():
    pool = [
        {
            "candidate_id": "a",
            "center": [0.0, 0, 0],
            "atom_indices": [0, 1],
            "geometry_score": 1.0,
            "chemical_feature": 0.0,
            "score": 1.0,
        },
        {
            "candidate_id": "b",
            "center": [6.0, 0, 0],
            "atom_indices": [2, 3],
            "geometry_score": 0.6,
            "chemical_feature": 1.0,
            "score": 1.2,
        },
    ]
    zero = select_generation(pool, GenerationCalibration(chemical_weight=0))
    strong = select_generation(pool, GenerationCalibration(chemical_weight=2))
    assert zero[0]["candidate_id"] == "a" and strong[0]["candidate_id"] == "b"
    for p in zero + strong:
        original = next(x for x in pool if x["candidate_id"] == p["candidate_id"])
        assert (
            p["center"] == original["center"]
            and p["atom_indices"] == original["atom_indices"]
        )


@pytest.mark.parametrize(
    "fields",
    [
        {"center_separation_A": 0},
        {"growth_radius_A": -1},
        {"growth_min_fraction": 0},
        {"growth_min_fraction": 1.1},
        {"chemical_weight": -1},
        {"chemical_weight": float("nan")},
        {"region_budget": 4},
    ],
)
def test_invalid_generation_settings_fail_explicitly(fields):
    with pytest.raises(ValueError):
        GenerationCalibration(**fields)
