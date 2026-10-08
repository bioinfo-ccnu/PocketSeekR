"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

from dataclasses import dataclass

import math

import numpy as np

PROTOCOL = {
    "version": "physics_v1",
    "score_units": "empirical_score_not_free_energy",
    "orientation_count": 48,
    "initial_translation_A": 1.5,
    "translation_bound_A": 2.0,
    "refinement_translation_steps_A": [0.5, 0.25],
    "refinement_rotation_steps_deg": [15.0, 7.5],
    "pair_cutoff_A": 8.0,
    "gauss1_weight": -0.035579,
    "gauss2_weight": -0.005156,
    "repulsion_weight": 0.840245,
    "hbond_weight": -0.587439,
    "stacking_weight": -0.587439,
    "hbond_heavy_distance_max_A": 3.5,
    "hbond_heavy_distance_min_A": 2.2,
    "donor_hydrogen_distance_A": 1.0,
    "hbond_angle_min_deg": 120.0,
    "acceptor_direction_max_deviation_deg": 60.0,
    "stacking_normal_distance_range_A": [2.8, 4.5],
    "stacking_normal_distance_optimum_A": 3.5,
    "stacking_normal_distance_width_A": 0.6,
    "stacking_plane_angle_max_deg": 30.0,
    "stacking_lateral_offset_max_A": 2.5,
    "stacking_lateral_width_A": 2.0,
    "ring_planarity_rms_max_A": 0.25,
    "polar_matching": "one_pair_per_donor_heavy_atom_and_acceptor",
    "ligand_polar_typing": "RDKit_BaseFeatures_actual_H_donors_NOS_acceptors",
    "candidate_pool": "unchanged_geometry_candidates_no_refill",
}

@dataclass
class Site:
    atom: int
    directions: np.ndarray
    # An unresolved single-bond azimuth is represented as a chemical cone;
    # never fabricate a laboratory-axis hydrogen or lone pair.
    cone_deg: float | None = None

@dataclass
class Ring:
    center: np.ndarray
    normal: np.ndarray

@dataclass
class Geometry:
    xyz: np.ndarray
    radii: np.ndarray
    donors: list[Site]
    acceptors: list[Site]
    rings: list[Ring]
    unavailable_polar_sites: int = 0

def _unit(vector):
    length = np.linalg.norm(vector)
    return None if length < 1e-8 else vector / length

def _ring(points):
    if len(points) < 3:
        return None
    center = points.mean(0)
    _, singular, axes = np.linalg.svd(points - center, full_matrices=False)
    if (
        singular[1] < 1e-8
        or singular[-1] / math.sqrt(len(points)) > PROTOCOL["ring_planarity_rms_max_A"]
    ):
        return None
    return Ring(center, axes[-1])

def _outward(xyz, index, neighbors):
    if not neighbors:
        return None
    directions = [_unit(xyz[index] - xyz[j]) for j in neighbors]
    if any(d is None for d in directions):
        return None
    return _unit(np.sum(directions, axis=0))

def _lone_pairs(xyz, index, neighbors, sp2, normal):
    axis = _outward(xyz, index, neighbors)
    if axis is None:
        return None
    if sp2 and len(neighbors) >= 2:
        return Site(index, axis[None])
    if sp2 and normal is not None:
        side = _unit(np.cross(normal, axis))
        if side is not None:
            return Site(
                index,
                np.stack(
                    [0.5 * axis + sign * math.sqrt(0.75) * side for sign in (-1, 1)]
                ),
            )
    if len(neighbors) >= 2:
        side = _unit(
            np.cross(xyz[neighbors[0]] - xyz[index], xyz[neighbors[1]] - xyz[index])
        )
        if side is not None:
            angle = math.radians(54.75)
            return Site(
                index,
                np.stack(
                    [
                        math.cos(angle) * axis + sign * math.sin(angle) * side
                        for sign in (-1, 1)
                    ]
                ),
            )
    return Site(index, axis[None], cone_deg=60.0 if sp2 else 54.75)

_PURINE = [
    ("N9", "C8"),
    ("C8", "N7"),
    ("N7", "C5"),
    ("C5", "C4"),
    ("C4", "N9"),
    ("C4", "N3"),
    ("N3", "C2"),
    ("C2", "N1"),
    ("N1", "C6"),
    ("C6", "C5"),
]

_PYRIMIDINE = [
    ("N1", "C2"),
    ("C2", "N3"),
    ("N3", "C4"),
    ("C4", "C5"),
    ("C5", "C6"),
    ("C6", "N1"),
]

_SUGAR = [
    ("C1'", "O4'"),
    ("O4'", "C4'"),
    ("C4'", "C3'"),
    ("C3'", "C2'"),
    ("C2'", "C1'"),
    ("C2'", "O2'"),
    ("C3'", "O3'"),
    ("C4'", "C5'"),
    ("C5'", "O5'"),
    ("O5'", "P"),
    ("P", "OP1"),
    ("P", "OP2"),
    ("P", "O1P"),
    ("P", "O2P"),
]

_BASE = {
    "A": (_PURINE + [("C6", "N6")], {"N6"}, {"N1", "N3", "N7"}),
    "G": (_PURINE + [("C6", "O6"), ("C2", "N2")], {"N1", "N2"}, {"O6", "N3", "N7"}),
    "I": (_PURINE + [("C6", "O6")], {"N1"}, {"O6", "N3", "N7"}),
    "C": (_PYRIMIDINE + [("C2", "O2"), ("C4", "N4")], {"N4"}, {"N3", "O2"}),
    "U": (_PYRIMIDINE + [("C2", "O2"), ("C4", "O4")], {"N3"}, {"O2", "O4"}),
    "PSU": (_PYRIMIDINE + [("C2", "O2"), ("C4", "O4")], {"N1", "N3"}, {"O2", "O4"}),
    "H2U": (_PYRIMIDINE + [("C2", "O2"), ("C4", "O4")], {"N3"}, {"O2", "O4"}),
}

def rna_geometry(atoms, positions, radii):
    xyz = np.asarray(positions, dtype=float)
    groups = {}
    for i, atom in enumerate(atoms):
        key = atom.chain_id, atom.residue_id, atom.residue_name
        names = groups.setdefault(key, {})
        if atom.atom_name in names:
            raise ValueError(
                f"duplicate RNA atom in chemical typing: {key}/{atom.atom_name}"
            )
        names[atom.atom_name] = i
    donors, acceptors, rings, unavailable = [], [], [], 0
    for (_, _, resname), names in groups.items():
        base = resname[1:] if resname in {"RA", "RU", "RG", "RC"} else resname
        if base not in _BASE:
            raise ValueError(f"unsupported RNA chemistry for physics_v1: {resname}")
        bonds, donor_names, acceptor_names = _BASE[base]
        glycoside = "C5" if base == "PSU" else "N9" if base in {"A", "G", "I"} else "N1"
        bonds = bonds + _SUGAR + [("C1'", glycoside)]
        neighbors = {i: [] for i in names.values()}
        expected = {name: set() for name in names}
        for a, b in bonds:
            if a in expected:
                expected[a].add(b)
            if b in expected:
                expected[b].add(a)
            if a in names and b in names:
                neighbors[names[a]].append(names[b])
                neighbors[names[b]].append(names[a])
        ring_names = {
            name
            for pair in (_PURINE if base in {"A", "G", "I"} else _PYRIMIDINE)
            for name in pair
        }
        # Missing ring atoms must not create a fabricated aromatic plane.
        plane = (
            _ring(xyz[[names[n] for n in sorted(ring_names)]])
            if ring_names <= names.keys()
            else None
        )
        if plane is not None and base != "H2U":
            rings.append(plane)
        normal = plane.normal if plane is not None else None
        for name, i in names.items():
            bonded = neighbors[i]
            polar = (
                name in donor_names
                or name in acceptor_names
                or name in {"O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O1P", "O2P"}
            )
            if polar and not expected[name] <= names.keys():
                unavailable += int(name in donor_names or name == "O2'") + int(
                    name in acceptor_names or name.startswith("O")
                )
                continue
            if name in donor_names or name == "O2'":
                axis = _outward(xyz, i, bonded)
                site = None
                if axis is not None:
                    if name == "O2'":
                        site = Site(i, axis[None], cone_deg=70.5)
                    elif name in {"N2", "N4", "N6"}:
                        side = (
                            _unit(np.cross(normal, axis))
                            if normal is not None
                            else None
                        )
                        if side is not None:
                            site = Site(
                                i,
                                np.stack(
                                    [
                                        0.5 * axis + s * math.sqrt(0.75) * side
                                        for s in (-1, 1)
                                    ]
                                ),
                            )
                    elif len(bonded) == 2:
                        site = Site(i, axis[None])
                if site is None:
                    unavailable += 1
                else:
                    donors.append(site)
            if name in acceptor_names or name in {
                "O2'",
                "O3'",
                "O4'",
                "O5'",
                "OP1",
                "OP2",
                "O1P",
                "O2P",
            }:
                # For 3'-phosphodiester O, include the next residue's covalent P.
                if name == "O3'":
                    bonded = bonded + [
                        j
                        for j, a in enumerate(atoms)
                        if a.atom_name == "P"
                        and a.chain_id == atoms[i].chain_id
                        and j not in bonded
                        and np.linalg.norm(xyz[i] - xyz[j]) < 1.9
                    ]
                site = _lone_pairs(xyz, i, bonded, name in acceptor_names, normal)
                if site is None:
                    unavailable += 1
                else:
                    acceptors.append(site)
    return Geometry(
        xyz, np.asarray(radii, dtype=float), donors, acceptors, rings, unavailable
    )

def _angular_cos(site, rotation, approach):
    if site.cone_deg is None:
        return float(np.max((site.directions @ rotation.T) @ approach))
    cosine = np.clip(np.dot(site.directions[0] @ rotation.T, approach), -1, 1)
    return math.cos(math.acos(cosine) - math.radians(site.cone_deg))
