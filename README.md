# PocketSeekR

**Find RNA small-molecule binding pockets from structure alone.**

[![CI](https://github.com/bioinfo-ccnu/PocketSeekR/actions/workflows/ci.yml/badge.svg)](https://github.com/bioinfo-ccnu/PocketSeekR/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/bioinfo-ccnu/PocketSeekR?color=416788)](https://github.com/bioinfo-ccnu/PocketSeekR/releases/latest)
[![Python](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![pip](https://img.shields.io/badge/pip-installable-3776AB?logo=pypi&logoColor=white)](#installation)
[![CPU](https://img.shields.io/badge/CPU-supported-52796F)](#quick-start)
[![Input](https://img.shields.io/badge/Input-RNA%20structure-52796F)](#input-and-output)

English | [简体中文](README.zh-CN.md)

PocketSeekR identifies and ranks candidate small-molecule binding pockets from an RNA PDB structure. It combines probe-accessible geometry, regional aromatic and directional polar support, and overlap-aware selection to return up to **five pockets** from an internal pool of at most **40 regions**. Prediction runs on CPU and requires neither a ligand nor a trained model.

```text
RNA structure  →  Accessible regions  →  Regional chemical support  →  Ranked pockets
```

## Installation

Install the tagged version directly from GitHub:

```bash
python -m pip install "git+https://github.com/bioinfo-ccnu/PocketSeekR.git@v0.1.0"
pocketseekr --version
```

Alternatively, download the wheel from [Releases](https://github.com/bioinfo-ccnu/PocketSeekR/releases/latest) and install it with `python -m pip install pocketseekr-0.1.0-py3-none-any.whl`.

Python 3.10 or newer is required. Runtime dependencies are NumPy, SciPy, PyYAML and PyTorch. PyTorch supplies coordinate tensors; no CUDA device, RDKit, model checkpoint or docking package is required. On Linux, a CPU-only PyTorch build can be installed before PocketSeekR:

```bash
python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
```

## Quick start

```bash
pocketseekr --rna-pdb path/to/rna.pdb --output-dir outputs/my_rna
```

The equivalent module command is:

```bash
python -m pocketseekr --rna-pdb path/to/rna.pdb --output-dir outputs/my_rna
```

To run the included example, clone the repository:

```bash
git clone https://github.com/bioinfo-ccnu/PocketSeekR.git
cd PocketSeekR
pocketseekr --rna-pdb examples/1F1T_rna.pdb --output-dir outputs/1F1T
```

## Python API

```python
from pocketseekr import predict

result = predict("path/to/rna.pdb")
for pocket in result["pockets"]:
    print(pocket["rank"], pocket["center"], pocket["score"])
```

The API returns the ranked pockets and the internal `region_pool`. Coordinates are in ångströms in the input structure's coordinate frame. Settings can also be passed explicitly:

```python
from pocketseekr import PocketSettings, predict

settings = PocketSettings()
result = predict("path/to/rna.pdb", settings=settings)
```

The CLI accepts the same settings as YAML through `--config`; [configs/default.yaml](configs/default.yaml) records the shipped defaults. Defaults are part of the installed package, so prediction works outside the repository directory.

## Input and output

| Item | Description |
| --- | --- |
| Input | RNA coordinates in PDB format; the first structural model is used. |
| `pockets.json` | Ranked pocket centers, scores, selection gains, RNA crop atom identities, settings and source hashes. |
| `candidates.json` | Internal candidate regions and the selected pockets. |
| Pocket context | A 10 Å center-based crop expanded to complete observed nucleotides. |

Only supported RNA heavy atoms contribute to prediction, even if other molecules occur in the input. Supported residue names are `A`, `U`, `G`, `C`, `I`, `RA`, `RU`, `RG`, `RC`, `PSU` and `H2U`; other residue names are excluded. Alternative locations are selected atom-wise by occupancy, preferring blank/A labels at equal occupancy. Direct mmCIF input is not supported in this release; convert it to a valid PDB structure before running the tool. Atom indices in JSON refer to the parsed RNA order, not PDB serial numbers. A context crop is distinct from the geometric cavity boundary.

## How it works

1. Search a 1.5 Å grid in an RNA-defined frame and connect probe-accessible voxels.
2. Grow and geometrically rank candidate regions, retaining at most 40 after separation and crop filtering.
3. Evaluate aromatic-face and directional polar support across each region, including local cross-support through accessible neighbors.
4. Combine geometric and chemical quality as `Q(r) = G(r) × [1 + 2c(r)]`.
5. Select up to five pockets with `gain(r | S) = Q(r) × [1 − Jmax(r,S)]^0.25`, limiting redundant crops.

See [Method](docs/method.md) for the geometry, scoring and selection definitions. Chemical support is a structural descriptor, not a calculated binding free energy.

## Development

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python -m build
python -m twine check --strict dist/*
```

CI runs the tests on Python 3.10, 3.12 and 3.13, builds a source distribution and wheel, and exercises the installed wheel outside the checkout. [Source provenance](docs/provenance.md) records the extraction from [RiboPoseDiff](https://github.com/wangleiofficial/RiboPoseDiff); [the source manifest](docs/source_extraction.json) identifies the computational definitions. Large datasets, historical evaluation outputs and docking/training code are outside this package.
