"""Command-line interface for the standalone PocketSeekR tool."""

import argparse
from pathlib import Path

import torch

from .api import predict
from .pocket_chemistry_calibration import ChemistryCalibration
from .utils.files import sha256_file, write_json_atomic


def build_parser():
    parser = argparse.ArgumentParser(description='Identify RNA ligand-binding pockets with PocketSeekR.')
    parser.add_argument('--rna-pdb', required=True, help='Input structure in PDB format.')
    parser.add_argument('--output-dir', required=True, help='Directory for JSON reports.')
    parser.add_argument('--config', help='Explicit regional chemistry YAML; defaults to the adopted 40-region settings.')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    torch.set_num_threads(1)
    settings = ChemistryCalibration.load(args.config) if args.config else ChemistryCalibration()
    report = predict(args.rna_pdb, settings)
    if args.config:
        report['config_sha256'] = sha256_file(args.config)
    candidates = {
        'rna_atom_count': report['rna_atom_count'],
        'region_pool': report.pop('region_pool'),
        'pockets': [
            {key: value for key, value in pocket.items() if key not in {'rank', 'crop_atom_identities'}}
            for pocket in report['pockets']
        ],
    }
    output = Path(args.output_dir)
    write_json_atomic(report, output/'pockets.json')
    write_json_atomic(candidates, output/'candidates.json')
    print(f"RNA atoms: {report['rna_atom_count']}; pockets: {len(report['pockets'])}; saved {output/'pockets.json'}")
    return 0
