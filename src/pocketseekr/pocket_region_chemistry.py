"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

import numpy as np

from scipy.spatial import cKDTree

from .pocket_regions import probe_segments_clear

REGION_CHEMISTRY_PROTOCOL = {
    "version": "RNA_region_chemistry_v1",
    "domain": "all_original_v3_region_voxels_no_resampling_or_center_replacement",
    "measure": "uniform_voxel_volume_normalized_by_original_region_volume",
    "point_fields": "existing_visible_base_face_and_base_polar_saturated_support",
    "cooperation_radius_A": 3.0,
    "cooperation_kernel": "exp_minus_half_squared_distance_over_radius_within_radius",
    "cooperation_visibility": "full_segment_clear_of_Bondi_plus_v3_probe_spheres",
    "kernel_normalization": "row_stochastic_with_unit_self_weight",
    "mean_feature": "half_mean_face_plus_mean_polar",
    "co_support": "mean_half_sqrt_face_times_local_polar_plus_sqrt_polar_times_local_face",
    "cooperative_feature": "half_sum_mean_feature_and_co_support",
    "score": "v3_score_times_1_plus_selected_region_feature_optional_1_plus_center_tensor",
    "missing_region": "error_no_center_scoring_substitution",
}

def regional_support_summary(points, face, polar, geometry, tree, probe):
    """Integrate bounded fields; local cross-support cannot jump an RNA wall."""
    points = np.asarray(points, dtype=float)
    face, polar = np.asarray(face, dtype=float), np.asarray(polar, dtype=float)
    count = len(points)
    if (
        points.shape != (count, 3)
        or count == 0
        or face.shape != (count,)
        or polar.shape != (count,)
        or not np.isfinite(points).all()
        or not np.isfinite(face).all()
        or not np.isfinite(polar).all()
        or np.any((face < 0) | (face > 1) | (polar < 0) | (polar > 1))
        or len(np.unique(points, axis=0)) != count
    ):
        raise ValueError(
            "regional chemistry requires nonempty finite points and bounded fields"
        )
    radius = REGION_CHEMISTRY_PROTOCOL["cooperation_radius_A"]
    pairs = cKDTree(points).query_pairs(radius + 1e-8, output_type="ndarray")
    if len(pairs):
        clear = probe_segments_clear(
            points[pairs[:, 0]],
            points[pairs[:, 1]],
            geometry.xyz,
            geometry.radii,
            tree,
            probe,
        )
        pairs = pairs[clear]
    # Unit diagonal is part of the kernel, including for a voxel with no other
    # directly accessible neighbor; it is not a fallback score.
    row_sum = np.ones(count)
    local_face, local_polar = face.copy(), polar.copy()
    if len(pairs):
        a, b = pairs.T
        distance2 = ((points[a] - points[b]) ** 2).sum(1)
        weights = np.exp(-0.5 * distance2 / radius**2)
        np.add.at(row_sum, a, weights)
        np.add.at(row_sum, b, weights)
        np.add.at(local_face, a, weights * face[b])
        np.add.at(local_face, b, weights * face[a])
        np.add.at(local_polar, a, weights * polar[b])
        np.add.at(local_polar, b, weights * polar[a])
    local_face /= row_sum
    local_polar /= row_sum
    co_support = 0.5 * (np.sqrt(face * local_polar) + np.sqrt(polar * local_face))
    base_support = 0.5 * (face + polar)
    mean_feature = float(base_support.mean())
    cooperative_feature = 0.5 * (mean_feature + float(co_support.mean()))

    def participation(values):
        squared_sum = float(values @ values)
        return float(values.sum() ** 2 / (count * squared_sum)) if squared_sum else 0.0

    support_total = float(base_support.sum())
    centroid, covariance = None, None
    if support_total > 0:
        centroid = np.sum(points * base_support[:, None], axis=0) / support_total
        delta = points - centroid
        covariance = (
            np.einsum("i,ij,ik->jk", base_support, delta, delta) / support_total
        )
    return {
        "face_support": float(face.mean()),
        "polar_support": float(polar.mean()),
        "mean_feature": mean_feature,
        "co_support": float(co_support.mean()),
        "cooperative_feature": cooperative_feature,
        "face_effective_volume_fraction": participation(face),
        "polar_effective_volume_fraction": participation(polar),
        "co_support_effective_volume_fraction": participation(co_support),
        "support_centroid_A": None if centroid is None else centroid.tolist(),
        "support_covariance_A2": None if covariance is None else covariance.tolist(),
        "accessible_neighbor_pair_count": int(len(pairs)),
        "local_face_support": local_face.tolist(),
        "local_polar_support": local_polar.tolist(),
        "point_co_support": co_support.tolist(),
    }
