"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

import math

import numpy as np

from .pocket_physics import _angular_cos, PROTOCOL

from .pocket_regions import probe_segments_clear

DIRECTIONAL_PROTOCOL = {
    "version": "RNA_directional_v1",
    "candidate_pool": "fixed_v3_after_crop_filter_and_deduplication",
    "directions": "26_RNA_frame_rays_13_antipodal_pairs",
    "ray_spheres": "Bondi_plus_v3_probe_radius",
    "ray_cutoff": "v3_shell_cutoff_atom_center_distance",
    "pair_weight": "exp_minus_mean_first_intersection_distance_over_length",
    "confinement_length_A": 3.0,
    "tensor": "3_over_13_sum_pair_weight_outer_axis",
    "tensor_feature": "sqrt_product_largest_two_eigenvalues",
    "plane_fit": "existing_RNA_complete_ring_SVD",
    "ring_planarity_rms_max_A": PROTOCOL["ring_planarity_rms_max_A"],
    "face_height_A": 3.5,
    "face_height_width_A": 1.0,
    "face_lateral_width_A": 2.0,
    "face_visibility_launch_height_A": 2.0,
    "polar_distance_A": 3.0,
    "polar_distance_width_A": 1.0,
    "visibility_surface_offset_A": 0.05,
    "chemical_cutoff_A": 7.5,
    "polar_sites": "base_atoms_only_existing_neutral_RNA_typing",
    "polar_alignment": "positive_best_direction_cosine_squared_including_chemical_cones",
    "visibility": "entire_centerline_segment_clear_of_unexpanded_RNA_spheres",
    "chemical_feature": "mean_of_one_minus_exp_minus_face_sum_and_polar_sum",
    "score": "v3_score_times_1_plus_tensor_times_1_plus_chemical_for_enabled_terms",
    "missing_chemistry": "omit_unresolved_sites_record_counts_no_fabricated_plane_or_direction",
}

def line_visible(start, end, geometry, tree):
    """Actual segment collision, not just sampling endpoints or midpoint."""
    if np.linalg.norm(end - start) < 1e-8:
        return False
    return bool(
        probe_segments_clear(
            start[None], end[None], geometry.xyz, geometry.radii, tree, 0.0
        )[0]
    )

def exposed_base_potential(center, atoms, geometry, tree):
    face_sum, polar_sum = 0.0, 0.0
    visible_faces, visible_sites = 0, 0
    cutoff = DIRECTIONAL_PROTOCOL["chemical_cutoff_A"]
    for ring in geometry.rings:
        delta = center - ring.center
        distance = float(np.linalg.norm(delta))
        height = float(delta @ ring.normal)
        if (
            distance > cutoff
            or abs(height) < DIRECTIONAL_PROTOCOL["face_visibility_launch_height_A"]
        ):
            continue
        lateral2 = max(0.0, distance**2 - height**2)
        start = (
            ring.center
            + math.copysign(
                DIRECTIONAL_PROTOCOL["face_visibility_launch_height_A"], height
            )
            * ring.normal
        )
        if not line_visible(start, center, geometry, tree):
            continue
        value = math.exp(
            -0.5
            * (
                (abs(height) - DIRECTIONAL_PROTOCOL["face_height_A"])
                / DIRECTIONAL_PROTOCOL["face_height_width_A"]
            )
            ** 2
            - 0.5 * lateral2 / DIRECTIONAL_PROTOCOL["face_lateral_width_A"] ** 2
        )
        face_sum += value
        visible_faces += 1
    donors = [site for site in geometry.donors if "'" not in atoms[site.atom].atom_name]
    acceptors = [
        site
        for site in geometry.acceptors
        if "'" not in atoms[site.atom].atom_name
        and atoms[site.atom].atom_name not in {"OP1", "OP2", "O1P", "O2P"}
    ]
    for site in donors + acceptors:
        delta = center - geometry.xyz[site.atom]
        distance = float(np.linalg.norm(delta))
        if (
            distance
            <= geometry.radii[site.atom]
            + DIRECTIONAL_PROTOCOL["visibility_surface_offset_A"]
            or distance > cutoff
        ):
            continue
        approach = delta / distance
        alignment = max(0.0, _angular_cos(site, np.eye(3), approach)) ** 2
        if alignment == 0:
            continue
        start = (
            geometry.xyz[site.atom]
            + (
                geometry.radii[site.atom]
                + DIRECTIONAL_PROTOCOL["visibility_surface_offset_A"]
            )
            * approach
        )
        if not line_visible(start, center, geometry, tree):
            continue
        polar_sum += alignment * math.exp(
            -0.5
            * (
                (distance - DIRECTIONAL_PROTOCOL["polar_distance_A"])
                / DIRECTIONAL_PROTOCOL["polar_distance_width_A"]
            )
            ** 2
        )
        visible_sites += 1
    face_support = -math.expm1(-face_sum)
    polar_support = -math.expm1(-polar_sum)
    return {
        "face_sum": face_sum,
        "polar_sum": polar_sum,
        "face_support": face_support,
        "polar_support": polar_support,
        "base_feature": 0.5 * (face_support + polar_support),
        "visible_base_faces": visible_faces,
        "visible_base_polar_sites": visible_sites,
        "typed_base_planes": len(geometry.rings),
        "typed_base_donors": len(donors),
        "typed_base_acceptors": len(acceptors),
        "unavailable_polar_sites_all_RNA": geometry.unavailable_polar_sites,
    }
