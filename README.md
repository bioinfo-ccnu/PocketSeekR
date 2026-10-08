# PocketSeekR

**Identifying RNA ligand-binding pockets through geometry and regional chemical support.**

PocketSeekR is a standalone, non-neural tool for identifying candidate ligand-binding pockets from an RNA structure. It takes observed RNA coordinates, searches probe-accessible free-space regions, evaluates regional chemical support, and selects up to five complementary pocket centers. No ligand coordinates or trained checkpoint are required.

## Calculate online

Open **[PocketSeekR on GitHub Pages](https://bioinfo-ccnu.github.io/PocketSeekR/)**, choose a local PDB file or try the 1F1T RNA example, then click **Identify pockets**. Inspect the ranked centers in the interactive RNA viewer and download the full JSON report, RNA-only PDB or selected RNA context crop.

Calculation runs locally in a Web Worker using Pyodide. The website downloads its Python runtime and numerical libraries from the same Pages site; it does not upload your structure or send it to a computation service. No server, account or trained checkpoint is needed. The first calculation loads about 30 MB of static assets and depends on network speed; later calculations reuse the initialized worker. Cancel stops the worker, and running again initializes a new one.

The browser accepts PDB files up to 20 MB and requires WebAssembly, module workers and WebGL. Only the first structural model is used, including for NMR ensembles. Larger structures may need more time and memory; use the Python CLI for inputs beyond the browser limit. Direct mmCIF input is not supported. The browser and CLI use the same prediction code and adopted settings.

## Install

Python 3.10 or later is required. A CPU environment is sufficient.

```bash
git clone https://github.com/bioinfo-ccnu/PocketSeekR.git
cd PocketSeekR
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Dependencies are NumPy, SciPy and PyYAML. PyTorch, RDKit and the original docking/training dependencies are not required.

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

RNA crops use a 10 Å neighborhood expanded to complete observed nucleotides. They are context crops, not cavity boundaries. The CLI exports JSON. The browser additionally exports RNA-only and crop PDB files with consecutive atom serial numbers; residue identities and observed coordinates are preserved.

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

The static browser application is deployed through GitHub Actions to [GitHub Pages](https://bioinfo-ccnu.github.io/PocketSeekR/). Each publication first runs real Chromium/WASM predictions and compares reports with native Python. See [build, test and deployment instructions](docs/browser_deployment.md).

## Provenance

PocketSeekR was extracted from the current RNA-only prediction path of [RiboPoseDiff](https://github.com/wangleiofficial/RiboPoseDiff). The original extraction preserved computational function bodies while removing unrelated docking and training modules. Version 0.2 replaces coordinate storage and cropping operations with NumPy to run the same core in native Python and WebAssembly, retaining float32 input coordinates, float64 geometry and adopted parameters. See [the extraction manifest](docs/source_extraction.json) for the source-file hashes and selected symbols.
