"""Browser bridge: the same prediction core with RNA-only PDB exports."""

from pathlib import Path
import re
from tempfile import TemporaryDirectory

from .api import predict
from .utils.io import _extract_rna_atoms, parse_pdb_atoms


def atoms_to_pdb(atoms) -> str:
    """Write selected observed atoms without ligands or synthetic crop boundaries."""
    lines = []
    for serial, atom in enumerate(atoms, 1):
        residue = re.fullmatch(r'(-?\d+)([A-Za-z]?)', atom.residue_id)
        if residue is None:
            raise ValueError(f'Cannot export PDB residue identity: {atom.residue_id}')
        number, insertion = int(residue[1]), residue[2]
        if serial > 99999 or not -999 <= number <= 9999:
            raise ValueError('Structure exceeds PDB serial/residue field limits')
        x, y, z = atom.position
        lines.append(
            f'ATOM  {serial:5d} {atom.atom_name:>4} {atom.residue_name:>3} '
            f'{atom.chain_id:1}{number:4d}{insertion:1}   '
            f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom.element:>2}'
        )
    return '\n'.join(lines + ['END']) + '\n'


def predict_text(pdb_text: str, filename: str = 'rna.pdb') -> dict:
    """Read text locally in Python/WASM; no network operation is performed."""
    with TemporaryDirectory(prefix='pocketseekr-') as directory:
        source = Path(directory)/'input.pdb'
        source.write_text(pdb_text, encoding='utf-8')
        report = predict(source)
        atoms = _extract_rna_atoms(parse_pdb_atoms(source), source)
        report['source_pdb'] = filename
        report['rna_pdb'] = atoms_to_pdb(atoms)
        for pocket in report['pockets']:
            pocket['crop_pdb'] = atoms_to_pdb([atoms[i] for i in pocket['atom_indices']])
        return report
