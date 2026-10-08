"""Frozen RNA geometry settings used by the adopted PocketSeekR method."""

from dataclasses import dataclass


@dataclass(frozen=True)
class StageAConfig:
    pocket_detector: str = 'cavity_grid'
    cavity_detector_version: int = 3
    pocket_ranker: str = 'region_cooperative_v1'
    cavity_region_field: str = 'opposed_rays'
    max_pockets: int = 5
    min_pocket_atoms: int = 12
    close_contact_cutoff: float = 3.25
    shell_cutoff: float = 7.5
    cavity_grid_spacing: float = 1.5
    cavity_min_clearance: float = 1.2
    cavity_max_clearance: float = 4.0
    cavity_center_separation: float = 6.0
