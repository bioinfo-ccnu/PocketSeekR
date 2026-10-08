"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

import math

def select_pockets(pockets, limit, *, crop_jaccard_limit=None):
    """Keep a ranked prefix, optionally suppressing strongly overlapping RNA crops.

    Crop overlap is a region-diversity approximation, not connected-cavity
    segmentation. Suppressed candidates are replaced only from the supplied
    pool; no new center is invented, and a short pool remains short.
    """
    if limit < 1:
        raise ValueError("final pocket limit must be positive")
    if crop_jaccard_limit is not None and (
        not math.isfinite(crop_jaccard_limit) or not 0 < crop_jaccard_limit <= 1
    ):
        raise ValueError("crop Jaccard limit must be in (0, 1]")
    selected, selected_sets = [], []
    for pocket in pockets:
        indices = pocket["atom_indices"].tolist()
        atoms = set(indices)
        if not atoms or len(atoms) != len(indices):
            raise ValueError("pocket crop must contain nonempty unique RNA indices")
        if crop_jaccard_limit is not None and any(
            len(atoms & other) / len(atoms | other) >= crop_jaccard_limit
            for other in selected_sets
        ):
            continue
        selected.append(pocket)
        selected_sets.append(atoms)
        if len(selected) == limit:
            break
    return selected
