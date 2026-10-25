from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import sys

import pytest
import torch
import yaml

from pocketseekr import predict
from pocketseekr.cli import main
from pocketseekr.pocket_chemistry_calibration import ChemistryCalibration
from pocketseekr.utils.io import parse_pdb_atoms

ROOT = Path(__file__).resolve().parents[1]


def test_cli_and_api_preserve_pool_and_crop_identity(tmp_path):
    source = ROOT/'examples/synthetic_shell.pdb'
    api = predict(source)
    assert main(['--rna-pdb', str(source), '--output-dir', str(tmp_path)]) == 0
    report = json.loads((tmp_path/'pockets.json').read_text())
    candidates = json.loads((tmp_path/'candidates.json').read_text())
    assert report['pockets'] == api['pockets']
    assert candidates['region_pool'] == api['region_pool']
    assert report['settings'] == asdict(ChemistryCalibration())
    assert len(report['pockets']) == 1
    assert len(candidates['region_pool']) <= 40
    for p in report['pockets']:
        assert len(p['crop_atom_identities']) == len(p['atom_indices'])
    assert yaml.safe_load((ROOT/'configs/default.yaml').read_text()) == report['settings']


def test_module_entrypoint_runs_outside_repository(tmp_path):
    result = subprocess.run([
        sys.executable, '-m', 'pocketseekr', '--rna-pdb',
        str(ROOT/'examples/synthetic_shell.pdb'), '--output-dir', str(tmp_path/'out'),
    ], cwd=tmp_path, capture_output=True, text=True, check=True)
    assert 'pockets: 1' in result.stdout


def test_legacy_settings_are_rejected(tmp_path):
    path = tmp_path/'legacy.yaml'
    path.write_text('region_budget: 20\nfinal_budget: 5\n')
    with pytest.raises(TypeError, match='final_budget'):
        ChemistryCalibration.load(path)


def test_no_supported_rna_is_rejected(tmp_path):
    path = tmp_path/'empty.pdb'
    path.write_text('END\n')
    with pytest.raises(ValueError, match='no RNA atoms'):
        predict(path)


def test_parser_preserves_first_model_and_occupancy_selection(tmp_path):
    def atom(x, alt, occupancy):
        return f'ATOM      1  C8 {alt}  A A   1    {x:8.3f}{0:8.3f}{0:8.3f}{occupancy:6.2f} 20.00           C'
    path = tmp_path/'models.pdb'
    path.write_text('\n'.join([
        'MODEL        1', atom(1, 'A', 0.4), atom(2, 'B', 0.6), 'ENDMDL',
        'MODEL        2', atom(99, ' ', 1), 'ENDMDL', 'END',
    ])+'\n')
    atoms = parse_pdb_atoms(path)
    assert len(atoms) == 1
    assert torch.equal(atoms[0].position, torch.tensor([2., 0., 0.]))
