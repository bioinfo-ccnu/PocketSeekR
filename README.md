# PocketSeekR

**Identifying RNA ligand-binding pockets through geometry and regional chemical support.**

PocketSeekR is a standalone, non-neural tool for identifying candidate ligand-binding pockets from an RNA structure. It takes observed RNA coordinates, searches probe-accessible free-space regions, evaluates regional chemical support, and selects up to five complementary pocket centers. No ligand coordinates or trained checkpoint are required.

## Install

Python 3.10 or later is required. A CPU environment is sufficient.

```bash
git clone https://github.com/bioinfo-ccnu/PocketSeekR.git
cd PocketSeekR
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Dependencies are NumPy, SciPy, PyYAML and PyTorch. PyTorch currently provides coordinate tensors and cropping operations; it is not used to load a trained model. RDKit and the original docking/training dependencies are not required.

## Identify pockets

```bash
pocketseekr --rna-pdb path/to/rna.pdb --output-dir outputs/my_rna

# Equivalent module invocation
python -m pocketseekr --rna-pdb path/to/rna.pdb --output-dir outputs/my_rna

# Explicit settings; omitted settings use the same adopted defaults
pocketseekr --rna-pdb path/to/rna.pdb --config configs/default.yaml --output-dir outputs/my_rna
```

The input must be a PDB file. The first structural model is used. Alternative locations are selected atom-wise by occupancy, with blank/A preference at equal occupancy. Supported residue names are A, U, G, C, I, RA, RU, RG, RC, PSU and H2U; unsupported residue names are excluded, and an input without supported RNA is rejected. Hydrogen/deuterium atoms are omitted. Supported modified residues use the explicit chemical typing inherited from the source algorithm. Direct mmCIF input is not implemented in this standalone release.

Two reports are saved:

- `pockets.json`: ordered final pockets with centers in Å, quality scores, selection gains, RNA crop atom identities, parameters and implementation hashes.
- `candidates.json`: the internal candidate pool and selected pockets. Atom indices refer to the parsed RNA atom order, not PDB atom serial numbers.

RNA crops use a 10 Å neighborhood expanded to complete observed nucleotides. They are context crops, not cavity boundaries. The CLI exports JSON; coordinate-file export and a browser viewer are not implemented in this release.

## Python API

```python
from pocketseekr import predict

report = predict("path/to/rna.pdb")
for pocket in report["pockets"]:
    print(pocket["rank"], pocket["center"], pocket["score"])
```

The API also returns `region_pool`. Empty searches return an empty pool and pocket list; candidates are never filled using another detector.

## Method

1. Build a 1.5 Å grid in a deterministic RNA frame and search probe-accessible free space using local clearance and 26 directional rays.
2. Grow connected regions through probe-clear adjacency. After geometric ranking, separation and crop filtering, retain at most **40** regions.
3. Evaluate visible aromatic-face and directional base-polar support throughout each region, including local chemical co-support.
4. Score regions using `Q(r) = G(r) × [1 + 2c(r)]`.
5. Greedily select up to **5** pockets with `gain(r | S) = Q(r) × [1 − Jmax(r,S)]^0.25`, subject to crop-overlap and center-distance exclusions.

The internal 40-region budget and final five-pocket limit are different. The default representative center remains determined by geometry. Chemical support is a structural heuristic, not a binding free energy or proof that a particular ligand can satisfy all proposed interactions.

See [method details](docs/method.md), [algorithm tables](docs/pocket_v5_algorithms.tex), [benchmark scope](docs/benchmark.md) and [source provenance](docs/provenance.md).

## Example and tests

```bash
pocketseekr --rna-pdb examples/synthetic_shell.pdb --output-dir outputs/example
python -m pip install -e '.[test]'
python -m pytest -q
```

The shell example is synthetic and tests geometry/integration; it is not biological validation. Real-structure fixtures from PDB entries 1F1T and 1NTA are used for inherited chemistry tests. Only supported RNA atoms contribute to prediction, even when an input file contains other molecules.

## Benchmark boundary

The source developmental comparison includes 940 ligand instances from 554 PDB entries (553 independent RNA inputs). The adopted method has Top1/3/5 localization counts of **314/565/652**, compared with **238/387/428** for the fixed fpocketR localization core. These are results on data used during method development and calibration, not a new independent test of this repository. RBind/Rsite comparisons include explicitly documented spatial adapters; GHECOM uses one fixed configuration.

## Deployment status

This repository provides a local Python CLI and API. A GitHub Pages/Pyodide implementation is planned but is not included or deployed. A Pages website alone does not execute the current Python package.

## Provenance

PocketSeekR was extracted from the current RNA-only prediction path of [RiboPoseDiff](https://github.com/wangleiofficial/RiboPoseDiff). The extraction preserves computational function bodies, coordinate precision, ordering and default parameters while removing unrelated docking and training modules. See [the extraction manifest](docs/source_extraction.json) for the source-file hashes and selected symbols.
