from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import eye
import torch

from pocketseekr.config import StageAConfig
from pocketseekr.pocket_chemistry_calibration import (
    ChemistryCalibration,
    ChemicalRegionalSearch,
    bounded_fields,
    chemical_descriptor,
    representative_index,
)
from pocketseekr.pocket_generation import RegionalPocketSearch
from pocketseekr.pocket_region_chemistry import regional_support_summary
from pocketseekr.utils.io import _extract_rna_atoms, parse_pdb_atoms


def test_saturation_uses_raw_sums_and_monotone_scale():
    raw = np.array([[0, 0], [1, 2], [40, 50.0]])
    assert np.all(bounded_fields(raw, 0.5) >= bounded_fields(raw, 2))
    assert bounded_fields(raw, 1)[0].tolist() == [0, 0]
    assert bounded_fields(raw, 2)[2, 0] < 1


def test_internal_weights_are_independent_and_bounded():
    fields = np.array([[1.0, 0], [1, 0]])
    base = ChemistryCalibration()
    assert (
        chemical_descriptor(
            fields, eye(2), replace(base, face_weight=1, cooperation_weight=0)
        )[0]
        == 1
    )
    assert (
        chemical_descriptor(
            fields, eye(2), replace(base, face_weight=0, cooperation_weight=0)
        )[0]
        == 0
    )
    assert (
        chemical_descriptor(fields, eye(2), replace(base, cooperation_weight=1))[0] == 0
    )


def test_center_is_a_member_and_chemical_support_moves_it():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0]])
    cfg = ChemistryCalibration()
    assert representative_index(points, np.ones(3), np.array([0, 0, 1.0]), cfg) == 1
    assert (
        representative_index(
            points,
            np.ones(3),
            np.array([0, 0, 1.0]),
            replace(cfg, center_chemical_weight=10),
        )
        == 2
    )


@pytest.mark.parametrize("pdb_id", ["1F1T", "1NTA"])
def test_real_rna_reproduces_parent_and_default_descriptor(pdb_id):
    path = f"tests/data/{pdb_id}.pdb"
    from pathlib import Path

    if not Path(path).exists():
        pytest.skip("local benchmark coordinates unavailable")
    atoms = _extract_rna_atoms(parse_pdb_atoms(path), path)
    xyz = torch.stack([a.position for a in atoms])
    stage = StageAConfig()
    search = ChemicalRegionalSearch(atoms, xyz, stage)
    cfg = ChemistryCalibration(face_weight=0.5, cooperation_radius_A=3.0)
    pool = search.pool(cfg)
    old = RegionalPocketSearch(atoms, xyz, stage).pool(cfg.generation())
    assert len(pool) == len(old)
    for p, q in zip(pool, old):
        for k in ("candidate_id", "center", "atom_indices", "voxel_id"):
            assert p[k] == q[k]
        assert p["chemical_feature"] == pytest.approx(q["chemical_feature"], abs=2e-15)
    c = search.regions(cfg)[0]
    members = c["_members"]
    fields = np.array([search.point_chemistry(i) for i in members])
    old_summary = regional_support_summary(
        search.grid[1][members],
        fields[:, 0],
        fields[:, 1],
        search.geometry,
        search.tree,
        stage.cavity_min_clearance,
    )
    feature, _ = search.descriptor(c, cfg)
    assert feature == pytest.approx(old_summary["cooperative_feature"], abs=2e-15)
    points = search.grid[1][members]
    kernel = search.kernel(members, 4.5)
    assert np.asarray(kernel.sum(1)).flatten() == pytest.approx(np.ones(len(points)))
    assert (kernel - kernel.T).shape == (len(points), len(points))
    for p in search.pool(replace(cfg, center_chemical_weight=2)):
        region = search.prepared_regions[
            replace(cfg.generation(), chemical_weight=1.0)
        ][p["parent_region_id"]]
        assert p["voxel_id"] in search.grid[0][region["_members"]]


@pytest.mark.parametrize(
    "fields",
    [
        {"face_weight": -1},
        {"cooperation_weight": 2},
        {"cooperation_radius_A": 0},
        {"saturation_scale": float("nan")},
        {"center_geometry_power": -1},
        {"center_chemical_weight": -1},
    ],
)
def test_invalid_settings_fail(fields):
    with pytest.raises(ValueError):
        ChemistryCalibration(**fields)


def test_cooperative_kernel_cannot_cross_rna_wall():
    from types import SimpleNamespace
    from scipy.spatial import cKDTree

    search = object.__new__(ChemicalRegionalSearch)
    search.grid = (None, np.array([[-2.5, 0.0, 0.0], [2.5, 0.0, 0.0]]))
    search.xyz = np.array([[0.0, 0.0, 0.0]])
    search.radii = np.array([1.0])
    search.tree = cKDTree(search.xyz)
    search.stage = SimpleNamespace(cavity_min_clearance=1.2)
    search.kernel_cache = {}
    blocked = search.kernel(np.array([0, 1]), 6.0)
    assert blocked.toarray() == pytest.approx(np.eye(2))
    search.xyz = np.array([[0.0, 10.0, 0.0]])
    search.tree = cKDTree(search.xyz)
    search.kernel_cache.clear()
    accessible = search.kernel(np.array([0, 1]), 6.0).toarray()
    assert accessible[0, 1] > 0
    assert accessible.sum(1) == pytest.approx(np.ones(2))
