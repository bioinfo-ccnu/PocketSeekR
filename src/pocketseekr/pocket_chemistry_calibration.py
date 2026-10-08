"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace

import math

import numpy as np

from scipy.sparse import coo_matrix

from scipy.spatial import cKDTree

from .pocket_directional import DIRECTIONAL_PROTOCOL, exposed_base_potential

from .pocket_generation import GenerationCalibration, RegionalPocketSearch, select_generation

from .pocket_regions import probe_segments_clear

from .pocket_selection import select_pockets

@dataclass(frozen=True)
class ChemistryCalibration(GenerationCalibration):
    chemical_weight: float = 2.0
    face_weight: float = 0.75
    cooperation_weight: float = 0.5
    cooperation_radius_A: float = 1.5
    saturation_scale: float = 1.0
    center_geometry_power: float = 1.0
    center_chemical_weight: float = 0.0

    def __post_init__(self):
        super().__post_init__()
        for value in (self.face_weight, self.cooperation_weight):
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("internal mixing weights must be in [0,1]")
        for value in (self.cooperation_radius_A, self.saturation_scale):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(
                    "cooperation radius and saturation scale must be positive"
                )
        for value in (self.center_geometry_power, self.center_chemical_weight):
            if not math.isfinite(value) or value < 0:
                raise ValueError("center weights must be finite and nonnegative")

    def generation(self):
        return GenerationCalibration(
            **{
                k: v
                for k, v in asdict(self).items()
                if k in GenerationCalibration.__dataclass_fields__
            }
        )

def bounded_fields(raw, scale):
    raw = np.asarray(raw, dtype=float)
    if (
        raw.ndim != 2
        or raw.shape[1] != 2
        or not np.isfinite(raw).all()
        or (raw < 0).any()
    ):
        raise ValueError("raw chemistry must be a finite nonnegative N x 2 array")
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("saturation scale must be positive")
    # Use raw sums; inversion of an already saturated field loses information.
    return -np.expm1(-raw / scale)

def chemical_descriptor(fields, kernel, config):
    face, polar = fields.T
    local_face, local_polar = (kernel @ fields).T
    mean = config.face_weight * face + (1 - config.face_weight) * polar
    cross = 0.5 * (np.sqrt(face * local_polar) + np.sqrt(polar * local_face))
    feature = (1 - config.cooperation_weight) * float(
        mean.mean()
    ) + config.cooperation_weight * float(cross.mean())
    return feature, mean

def representative_index(points, geometry_scores, mean_support, config):
    weights = np.asarray(geometry_scores) ** config.center_geometry_power
    weights = weights * (1 + config.center_chemical_weight * mean_support)
    if not np.isfinite(weights).all() or not (weights > 0).all():
        raise ValueError("representative center requires positive finite voxel weights")
    mean = np.average(points, axis=0, weights=weights)
    distance = np.linalg.norm(points - mean, axis=1)
    return int(np.flatnonzero(distance <= distance.min() + 1e-8)[0])

class ChemicalRegionalSearch(RegionalPocketSearch):
    """Keep the original geometry pool, then relocate/filter its representatives.

    All chemistry uses full original connected-region voxels. Re-centering does
    not search for replacement regions, increase the pool, or use ligand labels.
    """

    def __init__(self, atoms, positions, stage):
        super().__init__(atoms, positions, stage)
        self.raw_cache, self.kernel_cache, self.descriptor_cache = {}, {}, {}
        self.relocated_cache = {}
        self.base_pool_cache, self.prepared_regions, self.center_cache = {}, {}, {}

    def point_chemistry(self, index):
        if index not in self.raw_cache:
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
            values = exposed_base_potential(point, self.atoms, geometry, self.tree)
            self.raw_cache[index] = (values["face_sum"], values["polar_sum"])
            self.point_cache[index] = (values["face_support"], values["polar_support"])
        return self.point_cache[index]

    def kernel(self, members, radius):
        key = (tuple(members.tolist()), radius)
        if key not in self.kernel_cache:
            points = self.grid[1][members]
            pairs = cKDTree(points).query_pairs(radius + 1e-8, output_type="ndarray")
            if len(pairs):
                clear = probe_segments_clear(
                    points[pairs[:, 0]],
                    points[pairs[:, 1]],
                    self.xyz,
                    self.radii,
                    self.tree,
                    self.stage.cavity_min_clearance,
                )
                pairs = pairs[clear]
            n = len(points)
            if len(pairs):
                a, b = pairs.T
                weights = np.exp(
                    -0.5 * ((points[a] - points[b]) ** 2).sum(1) / radius**2
                )
                matrix = coo_matrix(
                    (
                        np.r_[np.ones(n), weights, weights],
                        (np.r_[np.arange(n), a, b], np.r_[np.arange(n), b, a]),
                    ),
                    shape=(n, n),
                ).tocsr()
            else:
                matrix = coo_matrix(
                    (np.ones(n), (np.arange(n), np.arange(n))), shape=(n, n)
                ).tocsr()
            self.kernel_cache[key] = matrix.multiply(
                1 / np.asarray(matrix.sum(1))
            ).tocsr()
        return self.kernel_cache[key]

    def descriptor(self, region, config):
        members = region["_members"]
        key = (
            tuple(members.tolist()),
            config.face_weight,
            config.cooperation_weight,
            config.cooperation_radius_A,
            config.saturation_scale,
        )
        if key not in self.descriptor_cache:
            for i in members:
                self.point_chemistry(i)
            fields = bounded_fields(
                [self.raw_cache[i] for i in members], config.saturation_scale
            )
            feature, mean = chemical_descriptor(
                fields, self.kernel(members, config.cooperation_radius_A), config
            )
            self.descriptor_cache[key] = (feature, mean)
        return self.descriptor_cache[key]

    def pool(self, config):
        if self.grid is None:
            return []
        generation = replace(config.generation(), chemical_weight=1.0)
        if generation not in self.base_pool_cache:
            # Frozen generator defines eligibility before center relocation.
            self.base_pool_cache[generation] = super().pool(generation)
            self.prepared_regions[generation] = {
                c["region"]["region_id"]: c for c in self.regions(generation)
            }
        pool = self.base_pool_cache[generation]
        regions = self.prepared_regions[generation]

        relocation_key = (
            generation,
            config.center_geometry_power,
            config.center_chemical_weight,
            config.face_weight if config.center_chemical_weight else None,
            config.saturation_scale if config.center_chemical_weight else None,
        )
        if relocation_key not in self.relocated_cache:

            def proposals():
                separated = []
                for p in pool:
                    region = regions[p["parent_region_id"]]
                    _, mean = self.descriptor(region, config)
                    center_key = (
                        p["parent_region_id"],
                        config.growth_radius_A,
                        config.growth_min_fraction,
                        config.center_geometry_power,
                        config.center_chemical_weight,
                        config.face_weight if config.center_chemical_weight else None,
                        config.saturation_scale
                        if config.center_chemical_weight
                        else None,
                    )
                    if center_key not in self.center_cache:
                        members = region["_members"]
                        j = representative_index(
                            self.grid[1][members], self.grid[2][members], mean, config
                        )
                        best = int(members[j])
                        c = region | {
                            "center": self.grid[1][best],
                            "region": region["region"]
                            | {"center_voxel_id": int(self.grid[0][best])},
                        }
                        self.center_cache[center_key] = (
                            int(self.grid[0][best]),
                            self.crop(c),
                            self.grid[1][best],
                        )
                    voxel, (center, indices), precise_center = self.center_cache[
                        center_key
                    ]
                    if any(
                        np.linalg.norm(precise_center - other)
                        < config.center_separation_A - 1e-8
                        for other in separated
                    ):
                        continue
                    separated.append(precise_center)
                    if len(indices) < self.stage.min_pocket_atoms:
                        continue
                    yield p | {
                        "center": center.tolist(),
                        "voxel_id": voxel,
                        "atom_indices": indices,
                        "kind": "regional_representative",
                    }

            relocated = [
                p | {"atom_indices": p["atom_indices"].tolist()}
                for p in select_pockets(
                    proposals(), config.region_budget, crop_jaccard_limit=0.8
                )
            ]
            self.relocated_cache[relocation_key] = relocated
        result = []
        for p in self.relocated_cache[relocation_key]:
            chemical, _ = self.descriptor(regions[p["parent_region_id"]], config)
            result.append(
                p
                | {
                    "chemical_feature": chemical,
                    "score": p["geometry_score"]
                    * (1 + config.chemical_weight * chemical),
                }
            )
        return result

def detect_chemistry_pockets(atoms, xyz, stage, config):
    pool = ChemicalRegionalSearch(atoms, xyz, stage).pool(config)
    return {
        "rna_atom_count": len(atoms),
        "region_pool": pool,
        "pockets": select_generation(pool, config),
    }
