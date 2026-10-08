from pathlib import Path
import numpy as np

from pocketseekr import predict
from pocketseekr.browser import predict_text
from pocketseekr.utils.io import _extract_rna_atoms, _rna_atom_key, parse_pdb_atoms

ROOT = Path(__file__).resolve().parents[1]


def test_browser_bridge_reuses_core_and_exports_exact_crop_atoms(tmp_path):
    source = ROOT/'tests/data/1F1T.pdb'
    native = predict(source)
    browser = predict_text(source.read_text(), '1F1T.pdb')
    assert browser['region_pool'] == native['region_pool']
    assert browser['source_pdb'] == '1F1T.pdb'
    atoms = _extract_rna_atoms(parse_pdb_atoms(source), source)
    full = tmp_path/'rna.pdb'
    full.write_text(browser['rna_pdb'])
    exported = parse_pdb_atoms(full)
    assert [_rna_atom_key(a) for a in exported] == [_rna_atom_key(a) for a in atoms]
    assert np.array_equal(np.stack([a.position for a in exported]), np.stack([a.position for a in atoms]))
    for original, pocket in zip(native['pockets'], browser['pockets']):
        assert {k:v for k,v in pocket.items() if k != 'crop_pdb'} == original
        crop = tmp_path/'crop.pdb'
        crop.write_text(pocket['crop_pdb'])
        assert [_rna_atom_key(a) for a in parse_pdb_atoms(crop)] == [
            _rna_atom_key(atoms[i]) for i in pocket['atom_indices']
        ]
