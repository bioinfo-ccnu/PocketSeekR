"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

from collections import deque

import numpy as np

from scipy.spatial import cKDTree

from .pocket_search import _ray_hits

REGION_PROTOCOL = {
    "version": "connected_regions_v1",
    "connectivity": "six_faces_with_probe_clear_segments",
    "minimum_opposed_ray_pairs": 1,
    "ray_pair_count": 13,
    "growth_radius_A": 6.0,
    "growth_minimum_seed_score_fraction": 0.5,
    "minimum_voxel_count": 3,
    "volume_support_radius_A": 3.0,
    "region_score": "peak_times_half_plus_half_saturating_probe_center_volume",
    "center": "free_voxel_nearest_score_weighted_region_mean",
    "crop_jaccard_limit": 0.8,
    "high_enclosure_min_mean_ray_fraction": 0.75,
}

def point_field(points, xyz, radii, bases, tree, config, frame):
    neighbors = tree.query_ball_point(points, config.shell_cutoff, return_sorted=True)
    counts = np.array([len(n) for n in neighbors])
    atoms = np.concatenate(neighbors).astype(int)
    rows = np.repeat(np.arange(len(points)), counts)
    vectors = xyz[atoms] - points[rows]
    distances = np.linalg.norm(vectors, axis=1)
    clearance = np.full(len(points), np.inf)
    np.minimum.at(clearance, rows, distances - radii[atoms])
    weight = np.exp(-0.5 * ((distances - 4.0) / 2.0) ** 2)
    total = np.bincount(rows, weights=weight, minlength=len(points))
    unit = vectors / distances[:, None]

    def mean(values):
        return np.bincount(rows, weights=weight * values, minlength=len(points)) / total

    balance = (
        1
        - np.linalg.norm(np.stack([mean(unit[:, i]) for i in range(3)], axis=1), axis=1)
    ).clip(0, 1)
    hits = _ray_hits(
        vectors, radii[atoms] + config.cavity_min_clearance, rows, len(points), frame
    )
    enclosure = hits.mean(1)
    pairs = (hits[:, :13] & hits[:, ::-1][:, :13]).sum(1)
    opposed = pairs / 13
    support = (0.75 + 0.25 * mean(bases[atoms])) * (1 - np.exp(-total / 12.0))
    score = (
        enclosure
        * support
        * (balance if config.cavity_region_field == "v2" else np.sqrt(opposed))
    )
    valid = (
        (clearance >= config.cavity_min_clearance)
        & (clearance <= config.cavity_max_clearance)
        & (counts >= 3)
        & (pairs >= REGION_PROTOCOL["minimum_opposed_ray_pairs"])
        & (score > 0)
    )
    return valid, score, enclosure, opposed, clearance

def scan_region_grid(xyz, radii, bases, config, origin, frame):
    spacing = config.cavity_grid_spacing
    local = (xyz - origin) @ frame
    margin = config.cavity_max_clearance
    low = np.floor((local.min(0) - margin) / spacing + 1e-8) * spacing
    high = np.ceil((local.max(0) + margin) / spacing - 1e-8) * spacing
    dimensions = np.rint((high - low) / spacing).astype(np.int64) + 1
    count = int(np.prod(dimensions))
    tree = cKDTree(xyz)
    chunks = []
    for start in range(0, count, 4096):
        ids = np.arange(start, min(start + 4096, count))
        coordinates = np.stack(
            [
                ids // (dimensions[1] * dimensions[2]),
                (ids // dimensions[2]) % dimensions[1],
                ids % dimensions[2],
            ],
            axis=1,
        )
        points = (low + coordinates * spacing) @ frame.T + origin
        nearest, _ = tree.query(points)
        near = (nearest >= radii.min() + config.cavity_min_clearance) & (
            nearest <= radii.max() + config.cavity_max_clearance
        )
        points, ids = points[near], ids[near]
        if not len(points):
            continue
        valid, score, enclosure, opposed, clearance = point_field(
            points, xyz, radii, bases, tree, config, frame
        )
        if valid.any():
            chunks.append(
                (
                    ids[valid],
                    points[valid],
                    score[valid],
                    enclosure[valid],
                    opposed[valid],
                    clearance[valid],
                )
            )
    if not chunks:
        return None
    return (*[np.concatenate([c[i] for c in chunks]) for i in range(6)], dimensions)

def probe_segments_clear(a, b, xyz, radii, tree, probe):
    """Check exact closest distance along each edge to probe-expanded RNA atoms."""
    segment = b - a
    squared_length = np.sum(segment * segment, axis=1)
    if np.any(squared_length <= 0):
        raise ValueError("probe graph edges must have positive length")
    half_length = np.sqrt(squared_length).max() / 2
    neighbors = tree.query_ball_point(
        (a + b) / 2, float(radii.max() + probe + half_length), return_sorted=True
    )
    counts = np.array([len(n) for n in neighbors])
    atoms = np.concatenate(neighbors).astype(int)
    rows = np.repeat(np.arange(len(a)), counts)
    vectors = xyz[atoms] - a[rows]
    projection = (np.sum(vectors * segment[rows], axis=1) / squared_length[rows]).clip(
        0, 1
    )
    distances = np.linalg.norm(vectors - projection[:, None] * segment[rows], axis=1)
    blocked = distances - radii[atoms] < probe - 1e-8
    return np.bincount(rows, weights=blocked, minlength=len(a)) == 0

def face_neighbors(ids, points, dimensions, xyz, radii, probe):
    """Six face neighbors, rejecting edges that cross an expanded RNA sphere."""
    adjacency = np.full((len(ids), 6), -1, dtype=np.int64)
    coordinates = np.stack(
        [
            ids // (dimensions[1] * dimensions[2]),
            (ids // dimensions[2]) % dimensions[1],
            ids % dimensions[2],
        ],
        axis=1,
    )
    strides = (dimensions[1] * dimensions[2], dimensions[2], 1)
    tree = cKDTree(xyz)
    for axis, stride in enumerate(strides):
        a = np.flatnonzero(coordinates[:, axis] + 1 < dimensions[axis])
        indices = np.searchsorted(ids, ids[a] + stride)
        inside = indices < len(ids)
        a, indices = a[inside], indices[inside]
        exists = ids[indices] == ids[a] + stride
        a, b = a[exists], indices[exists]
        if not len(a):
            continue
        # Endpoints and midpoints alone do not guarantee a collision-free edge.
        clear = probe_segments_clear(points[a], points[b], xyz, radii, tree, probe)
        a, b = a[clear], b[clear]
        adjacency[a, 2 * axis + 1] = b
        adjacency[b, 2 * axis] = a
    return adjacency

def grow_regions(points, scores, adjacency, *, radius, min_fraction, min_voxels):
    """Score-seeded connected regions with bounded spatial and saddle growth."""
    assigned = np.zeros(len(points), dtype=bool)
    regions = []
    for seed in np.argsort(-scores.round(12), kind="stable"):
        if assigned[seed]:
            continue
        queue = deque([int(seed)])
        assigned[seed] = True
        members = []
        threshold = scores[seed] * min_fraction
        while queue:
            point = queue.popleft()
            members.append(point)
            for neighbor in adjacency[point]:
                if neighbor < 0 or assigned[neighbor]:
                    continue
                if (
                    scores[neighbor] < threshold - 1e-12
                    or np.linalg.norm(points[neighbor] - points[seed]) > radius + 1e-8
                ):
                    continue
                assigned[neighbor] = True
                queue.append(int(neighbor))
        if len(members) >= min_voxels:
            regions.append((int(seed), np.array(sorted(members), dtype=int)))
    return regions
