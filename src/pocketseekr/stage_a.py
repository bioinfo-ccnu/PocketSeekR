"""Describe the fixed RNA-only regional detector, without docking variants."""

from copy import deepcopy

from .config import StageAConfig
from .pocket_directional import DIRECTIONAL_PROTOCOL
from .pocket_region_chemistry import REGION_CHEMISTRY_PROTOCOL
from .pocket_regions import REGION_PROTOCOL


def pocket_detector_config(config: StageAConfig) -> dict:
    if config != StageAConfig():
        raise ValueError('PocketSeekR requires the fixed RNA-only regional detector')
    directional = deepcopy(DIRECTIONAL_PROTOCOL)
    directional.update({
        'version': 'RNA_directional_region_v1',
        'chemical_feature': 'selected_feature_from_region_chemistry_protocol',
        'score': REGION_CHEMISTRY_PROTOCOL['score'],
    })
    return {
        'pocket_detector': 'cavity_grid',
        'cavity_detector_version': 3,
        'pocket_ranker': 'region_cooperative_v1',
        'directional_protocol': directional,
        'region_chemistry_protocol': deepcopy(REGION_CHEMISTRY_PROTOCOL),
        'cavity_region_field': 'opposed_rays',
        'region_protocol': deepcopy(REGION_PROTOCOL),
        **{name: getattr(config, name) for name in (
            'max_pockets', 'min_pocket_atoms', 'close_contact_cutoff',
            'shell_cutoff', 'cavity_grid_spacing', 'cavity_min_clearance',
            'cavity_max_clearance', 'cavity_center_separation',
        )},
    }
