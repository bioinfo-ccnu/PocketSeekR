"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

from dataclasses import dataclass

import math

from pathlib import Path

import numpy as np

import yaml

@dataclass(frozen=True)
class SelectionCalibration:
    region_budget: int = 40
    overlap_exponent: float = 0.25

    def __post_init__(self):
        if type(self.region_budget) is not int or self.region_budget < 5:
            raise ValueError("region_budget must be an integer at least five")
        if not math.isfinite(self.overlap_exponent) or self.overlap_exponent < 0:
            raise ValueError("overlap_exponent must be finite and nonnegative")

    @classmethod
    def load(cls, path: str | Path):
        return cls(**yaml.safe_load(Path(path).read_text()))

def select_calibrated(pool, config):
    """Q(r) * (1 - max crop Jaccard)**beta; retain v5 hard exclusions."""
    remaining = list(pool[: config.region_budget])
    if len({p["candidate_id"] for p in remaining}) != len(remaining):
        raise ValueError("candidate IDs must be unique")
    selected, selected_sets = [], []
    while remaining and len(selected) < 5:
        choices = []
        for pocket in remaining:
            atoms = set(pocket["atom_indices"])
            center = np.asarray(pocket["center"], dtype=float)
            if (
                not atoms
                or center.shape != (3,)
                or not np.isfinite(center).all()
                or not math.isfinite(pocket["score"])
                or pocket["score"] <= 0
            ):
                raise ValueError(
                    "selection requires finite centers and positive scores/crops"
                )
            overlap = max(
                (len(atoms & other) / len(atoms | other) for other in selected_sets),
                default=0.0,
            )
            if overlap >= 0.8 or any(
                np.linalg.norm(center - np.asarray(p["center"])) < 3.0 - 1e-8
                for p in selected
            ):
                continue
            gain = pocket["score"] * (1 - overlap) ** config.overlap_exponent
            choices.append((gain, pocket, atoms, overlap))
        if not choices:
            break
        gain, pocket, atoms, overlap = min(
            choices, key=lambda x: (-round(x[0], 12), x[1]["candidate_id"])
        )
        selected.append(
            pocket | {"selection_gain": gain, "selected_crop_overlap": overlap}
        )
        selected_sets.append(atoms)
        remaining = [
            p for p in remaining if p["candidate_id"] != pocket["candidate_id"]
        ]
    return selected
