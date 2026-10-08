"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

from dataclasses import dataclass, replace

import math

from pathlib import Path

import numpy as np

from scipy.spatial import cKDTree

import torch

import yaml

from .pocket_calibration import SelectionCalibration, select_calibrated

from .pocket_directional import DIRECTIONAL_PROTOCOL, exposed_base_potential

from .pocket_physics import rna_geometry

from .pocket_region_chemistry import regional_support_summary

from .pocket_regions import REGION_PROTOCOL, face_neighbors, grow_regions, scan_region_grid

from .pocket_search import _rna_frame

from .pocket_selection import select_pockets

from .stage_a import pocket_detector_config

from .utils.io import RNA_BASE_ATOMS, _atom_radii, _rna_atom_key

@dataclass(frozen=True)
class GenerationCalibration:
    center_separation_A: float = 4.0
    growth_radius_A: float = 6.0
    growth_min_fraction: float = 0.5
    region_budget: int = 40
    chemical_weight: float = 2.0
    overlap_exponent: float = 0.25

    def __post_init__(self):
        SelectionCalibration(self.region_budget, self.overlap_exponent)
        for value in (self.center_separation_A, self.growth_radius_A):
            if not math.isfinite(value) or value <= 0:
                raise ValueError("generation lengths must be positive and finite")
        if (
            not math.isfinite(self.growth_min_fraction)
            or not 0 < self.growth_min_fraction <= 1
        ):
            raise ValueError("growth_min_fraction must be in (0,1]")
        if not math.isfinite(self.chemical_weight) or self.chemical_weight < 0:
            raise ValueError("chemical_weight must be finite and nonnegative")

    @classmethod
    def load(cls, path: str | Path):
        return cls(**yaml.safe_load(Path(path).read_text()))

class RegionalPocketSearch:
    """One RNA grid and accessible adjacency, shared by explicit parameter sets."""

    def __init__(self, atoms, positions, stage):
        if (
            stage.pocket_detector != "cavity_grid"
            or stage.cavity_detector_version != 3
            or stage.pocket_ranker != "region_cooperative_v1"
            or stage.cavity_region_field != "opposed_rays"
            or stage.max_pockets != 5
        ):
            raise ValueError(
                "generation calibration requires the fixed v4 regional parent"
            )
        pocket_detector_config(stage)
        self.atoms, self.positions, self.stage = atoms, positions, stage
        self.xyz = positions.detach().cpu().double().numpy()
        self.radii = _atom_radii(atoms).double().numpy()
        frame = _rna_frame(self.xyz)
        if frame is None:
            raise ValueError("regional search requires a noncollinear RNA frame")
        self.grid = scan_region_grid(
            self.xyz,
            self.radii,
            np.array([a.atom_name in RNA_BASE_ATOMS for a in atoms]),
            stage,
            *frame,
        )
        self.adjacency = (
            None
            if self.grid is None
            else face_neighbors(
                self.grid[0],
                self.grid[1],
                self.grid[-1],
                self.xyz,
                self.radii,
                stage.cavity_min_clearance,
            )
        )
        self.tree = cKDTree(self.xyz)
        self.geometry = rna_geometry(atoms, self.xyz, self.radii)
        # Only sites outside the unchanged 7.5 A cutoff are omitted. Original
        # order and full-RNA occlusion are retained, including at boundaries.
        self.site_trees = []
        for sites, coordinates in (
            (self.geometry.rings, [r.center for r in self.geometry.rings]),
            (self.geometry.donors, [self.xyz[s.atom] for s in self.geometry.donors]),
            (
                self.geometry.acceptors,
                [self.xyz[s.atom] for s in self.geometry.acceptors],
            ),
        ):
            self.site_trees.append(
                (sites, cKDTree(np.asarray(coordinates).reshape(-1, 3)))
            )
        self.region_cache, self.crop_cache, self.point_cache, self.chemical_cache = (
            {},
            {},
            {},
            {},
        )

    def point_chemistry(self, index):
        if index not in self.point_cache:
            point = self.grid[1][index]
            local = [
                [
                    sites[i]
                    for i in tree.query_ball_point(
                        point,
                        DIRECTIONAL_PROTOCOL["chemical_cutoff_A"] + 1e-8,
                        return_sorted=True,
                    )
                ]
                for sites, tree in self.site_trees
            ]
            geometry = replace(
                self.geometry, rings=local[0], donors=local[1], acceptors=local[2]
            )
            f = exposed_base_potential(point, self.atoms, geometry, self.tree)
            self.point_cache[index] = (f["face_support"], f["polar_support"])
        return self.point_cache[index]

    def regions(self, config):
        key = (config.growth_radius_A, config.growth_min_fraction)
        if key in self.region_cache:
            return self.region_cache[key]
        ids, points, scores, enclosure, opposed, clearance, _ = self.grid
        grown = grow_regions(
            points,
            scores,
            self.adjacency,
            radius=config.growth_radius_A,
            min_fraction=config.growth_min_fraction,
            min_voxels=REGION_PROTOCOL["minimum_voxel_count"],
        )
        result = []
        volume_scale = 4 * math.pi * REGION_PROTOCOL["volume_support_radius_A"] ** 3 / 3
        for seed, members in grown:
            mean = np.average(points[members], axis=0, weights=scores[members])
            distance = np.linalg.norm(points[members] - mean, axis=1)
            best = members[np.flatnonzero(distance <= distance.min() + 1e-8)[0]]
            volume = len(members) * self.stage.cavity_grid_spacing**3
            peak = float(scores[members].max())
            fraction = float(enclosure[members].mean())
            result.append(
                {
                    "center": points[best],
                    "score": peak
                    * (0.5 + 0.5 * (1 - math.exp(-volume / volume_scale))),
                    "_members": members,
                    "region": {
                        "region_id": int(ids[seed]),
                        "center_voxel_id": int(ids[best]),
                        "type": "high_enclosure"
                        if fraction
                        >= REGION_PROTOCOL["high_enclosure_min_mean_ray_fraction"]
                        else "groove_like",
                        "voxel_count": int(len(members)),
                        "probe_center_volume_A3": float(volume),
                        "peak_point_score": peak,
                        "mean_point_score": float(scores[members].mean()),
                        "mean_ray_enclosure": fraction,
                        "mean_opposed_ray_fraction": float(opposed[members].mean()),
                        "center_clearance_A": float(clearance[best]),
                        "seed_center": points[seed].tolist(),
                        "score_weighted_mean": mean.tolist(),
                        "field": self.stage.cavity_region_field,
                    },
                }
            )
        result.sort(key=lambda p: (-round(p["score"], 12), p["region"]["region_id"]))
        self.region_cache[key] = result
        return result

    def crop(self, candidate):
        key = candidate["region"]["center_voxel_id"]
        if key not in self.crop_cache:
            center = torch.as_tensor(candidate["center"], dtype=self.positions.dtype)
            inside = torch.linalg.norm(self.positions - center, dim=-1) <= 10.0
            residues = {
                _rna_atom_key(a)[:3] for a, keep in zip(self.atoms, inside) if keep
            }
            indices = torch.tensor(
                [
                    i
                    for i, a in enumerate(self.atoms)
                    if _rna_atom_key(a)[:3] in residues
                ],
                dtype=torch.long,
            )
            self.crop_cache[key] = (center, indices)
        return self.crop_cache[key]

    def pool(self, config):
        if self.grid is None:
            return []

        def proposals():
            separated = []
            for c in self.regions(config):
                if any(
                    np.linalg.norm(c["center"] - p) < config.center_separation_A - 1e-8
                    for p in separated
                ):
                    continue
                separated.append(c["center"])
                center, indices = self.crop(c)
                if len(indices) >= self.stage.min_pocket_atoms:
                    yield c | {"center": center.tolist(), "atom_indices": indices}

        cropped = select_pockets(
            proposals(),
            config.region_budget,
            crop_jaccard_limit=REGION_PROTOCOL["crop_jaccard_limit"],
        )
        pool = []
        for rank, c in enumerate(cropped, 1):
            key = tuple(c["_members"].tolist())
            if key not in self.chemical_cache:
                values = np.array([self.point_chemistry(i) for i in c["_members"]])
                fields = regional_support_summary(
                    self.grid[1][c["_members"]],
                    values[:, 0],
                    values[:, 1],
                    self.geometry,
                    self.tree,
                    self.stage.cavity_min_clearance,
                )
                self.chemical_cache[key] = fields["cooperative_feature"]
            chemical = self.chemical_cache[key]
            pool.append(
                {
                    "candidate_id": f"r{rank:03d}_anchor",
                    "parent_geometry_rank": rank,
                    "parent_region_id": c["region"]["region_id"],
                    "voxel_id": c["region"]["center_voxel_id"],
                    "center": c["center"],
                    "atom_indices": c["atom_indices"].tolist(),
                    "score": c["score"] * (1 + config.chemical_weight * chemical),
                    "kind": "original_anchor",
                    "geometry_score": c["score"],
                    "chemical_feature": chemical,
                }
            )
        return pool

def select_generation(pool, config):
    weighted = [
        p
        | {
            "score": p["geometry_score"]
            * (1 + config.chemical_weight * p["chemical_feature"])
        }
        for p in pool
    ]
    return select_calibrated(
        weighted, SelectionCalibration(config.region_budget, config.overlap_exponent)
    )
