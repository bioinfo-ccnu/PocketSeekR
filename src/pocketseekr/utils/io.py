"""RNA-only pocket prediction components extracted from RiboPoseDiff."""

from __future__ import annotations

from dataclasses import dataclass

from pathlib import Path

import torch

from torch import Tensor

from .chemistry import BONDI_RADII

from .errors import SampleValidationError

RNA_RESIDUES = {"A", "U", "G", "C", "I", "RA", "RU", "RG", "RC", "PSU", "H2U"}

RNA_BASE_ATOMS = {
    "N1",
    "C2",
    "N3",
    "C4",
    "C5",
    "C6",
    "N7",
    "C8",
    "N9",
    "O2",
    "O4",
    "O6",
    "N2",
    "N4",
    "N6",
}

@dataclass
class ParsedAtom:
    record_type: str
    element: str
    atom_name: str
    residue_name: str
    chain_id: str
    residue_id: str
    position: Tensor

def _element_radii(elements: list[str]) -> Tensor:
    radii = BONDI_RADII
    unsupported = {element.upper() for element in elements} - radii.keys()
    if unsupported:
        raise SampleValidationError(f"unsupported vdW elements: {sorted(unsupported)}")
    return torch.tensor([radii[e.upper()] for e in elements])

def _atom_radii(atoms: list[ParsedAtom]) -> Tensor:
    return _element_radii([atom.element for atom in atoms])

def _rna_atom_key(atom: ParsedAtom) -> tuple[str, str, str, str]:
    return atom.chain_id, atom.residue_id, atom.residue_name, atom.atom_name

def parse_pdb_atoms(path: str | Path) -> list[ParsedAtom]:
    # One structural model and one alternative location per atom. Insertion codes
    # are part of residue identity, preventing silent collisions during alignment.
    selected: dict[tuple[str, str, str, str], tuple[float, str, ParsedAtom]] = {}
    model_seen = False
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if model_seen:
                    break
                model_seen = True
            if record == "ENDMDL":
                break
            if record not in {"ATOM", "HETATM"}:
                continue
            atom_name = line[12:16].strip()
            element = (
                line[76:78].strip().upper()
                or atom_name.lstrip("0123456789")[:1].upper()
            )
            if element in {"H", "D"}:
                continue
            try:
                position = torch.tensor(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
            except ValueError as exc:
                raise SampleValidationError(
                    f"invalid PDB coordinates in {path}: {line.rstrip()}"
                ) from exc
            if not torch.isfinite(position).all():
                raise SampleValidationError(f"non-finite coordinates in {path}")
            atom = ParsedAtom(
                record,
                element,
                atom_name,
                line[17:20].strip().upper(),
                line[21].strip(),
                line[22:26].strip() + line[26:27].strip(),
                position,
            )
            key = _rna_atom_key(atom)
            occupancy = float(line[54:60].strip() or "0")
            alt = line[16:17].strip()
            previous = selected.get(key)
            priority = (occupancy, alt in {"", "A"})
            if previous is None or priority > (previous[0], previous[1] in {"", "A"}):
                selected[key] = occupancy, alt, atom
    return [entry[2] for entry in selected.values()]

def _extract_rna_atoms(atoms: list[ParsedAtom], path: str | Path) -> list[ParsedAtom]:
    rna_atoms = [atom for atom in atoms if atom.residue_name in RNA_RESIDUES]
    if not rna_atoms:
        raise SampleValidationError(f"no RNA atoms found in {path}")
    return rna_atoms
