"""Public RNA-only prediction API."""

from dataclasses import asdict
from pathlib import Path

import torch

from .config import StageAConfig
from .pocket_chemistry_calibration import ChemistryCalibration, detect_chemistry_pockets
from .stage_a import pocket_detector_config
from .utils.files import sha256_file
from .utils.io import _extract_rna_atoms, parse_pdb_atoms


def predict(rna_pdb: str | Path, settings: ChemistryCalibration | None = None) -> dict:
    """Return an RNA-only report; empty searches remain empty.

    The first structural model and atom-wise occupancy/altloc selection follow
    the source parser. Only supported RNA residue names contribute to detection.
    No ligand coordinates, trained checkpoint or fpocket output is consumed.
    """
    source = Path(rna_pdb).resolve()
    settings = ChemistryCalibration() if settings is None else settings
    atoms = _extract_rna_atoms(parse_pdb_atoms(source), source)
    xyz = torch.stack([atom.position for atom in atoms])
    stage = StageAConfig()
    result = detect_chemistry_pockets(atoms, xyz, stage, settings)
    pockets = [
        pocket | {
            'rank': rank,
            'crop_atom_identities': [
                {
                    'chain_id': atoms[i].chain_id,
                    'residue_id': atoms[i].residue_id,
                    'residue_name': atoms[i].residue_name,
                    'atom_name': atoms[i].atom_name,
                }
                for i in pocket['atom_indices']
            ],
        }
        for rank, pocket in enumerate(result['pockets'], 1)
    ]
    package = Path(__file__).resolve().parent
    return {
        'method': 'PocketSeekR',
        'algorithm': 'v5_selection40',
        'protocol': {
            'inputs': 'RNA_only_no_native_ligand_no_learned_model',
            'crop': '10_A_center_sphere_expanded_to_complete_observed_nucleotides',
            'selection': 'greedy_Q_times_one_minus_max_crop_Jaccard_to_beta',
            'calibration': 'global_parameters_selected_offline_no_ligand_labels_at_prediction',
        },
        'settings': asdict(settings),
        'source_pdb': str(source),
        'source_sha256': sha256_file(source),
        'parent_detector': pocket_detector_config(stage),
        'implementation_sha256': {
            str(path.relative_to(package)): sha256_file(path)
            for path in sorted(package.rglob('*.py'))
        },
        'rna_atom_count': len(atoms),
        'pockets': pockets,
        'region_pool': result['region_pool'],
    }
