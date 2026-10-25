"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

import numpy as np

def _rna_frame(xyz):
    """Deterministic atom-index tie breaking in an RNA-defined rigid frame."""
    origin = xyz.mean(0)
    centered = xyz - origin
    norms = np.linalg.norm(centered, axis=1)
    if norms.max() < 1e-8:
        return None
    first = np.flatnonzero(norms >= norms.max() - 1e-8)[0]
    axis1 = centered[first] / norms[first]
    transverse = centered - np.outer(centered @ axis1, axis1)
    norms = np.linalg.norm(transverse, axis=1)
    if norms.max() < 1e-8:
        # A collinear receptor has no intrinsic 3-D search frame.
        return None
    second = np.flatnonzero(norms >= norms.max() - 1e-8)[0]
    axis2 = transverse[second] / norms[second]
    return origin, np.stack([axis1, axis2, np.cross(axis1, axis2)], axis=1)

def _ray_hits(vectors, expanded_radii, rows, point_count, frame):
    """26 blocked-ray flags in fixed order; opposite directions have reversed IDs."""
    axes = np.array(
        [
            [a, b, c]
            for a in (-1, 0, 1)
            for b in (-1, 0, 1)
            for c in (-1, 0, 1)
            if a or b or c
        ],
        dtype=float,
    )
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    axes = axes @ frame.T
    squared_distance = (vectors * vectors).sum(-1)
    blocked = np.zeros((point_count, len(axes)), dtype=bool)
    for i, axis in enumerate(axes):
        projection = vectors @ axis
        perpendicular = (squared_distance - projection**2).clip(0)
        hit = (projection > 0) & (perpendicular <= expanded_radii**2)
        blocked[:, i] = np.bincount(rows, weights=hit, minlength=point_count) > 0
    return blocked
