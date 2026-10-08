# Source and standalone extraction

The standalone package was extracted on 2026-10-08 from the working tree of RiboPoseDiff. The source HEAD was `fad5008b2ca6b8fcb70d6d843de0d957a9e43241`, but the source working tree contained uncommitted changes; that commit alone does not identify the complete exported implementation.

`source_extraction.json` identifies the exact source files using SHA256 and lists the extracted symbols. The computational definitions are copied verbatim, with only the package namespace/import dependency closure changed. The standalone `config.py` fixes the original default geometry. `stage_a.py` describes only the adopted regional detector, with a signature checked against the source default.

The independent CLI/API are new wrappers. They use the same first-model PDB parser, occupancy/altloc selection, supported RNA residues, coordinate dtype, geometry, chemistry, crop rules, ranking and default calibration. No docking model, training loop, native ligand label pipeline or legacy prediction entry point is included.

The inherited internal names containing `calibration` refer to parameter containers and fixed prediction calculations. This repository does not fit new parameters. `configs/default.yaml` records the adopted settings; the installed CLI defaults are provided by the same parameter dataclass and do not depend on a repository-relative config path.

The archived protocols document the original developmental experiments. Their RiboPoseDiff paths and commands describe the source repository and are not commands provided by this standalone package. Follow this repository's README for prediction.

Fixture sources: PDB entries [1F1T](https://www.rcsb.org/structure/1F1T) and [1NTA](https://www.rcsb.org/structure/1NTA). The fixtures are used as known structural inputs for numerical and chemical invariants, not as a held-out benchmark.

## Extraction verification

The source suite passed 330 tests, and the installed standalone wheel passed 34 tests. Full candidate pools and selected pockets match the source implementation exactly for 1F1T, 1NTA, 2GDI and 8K7W. All 51 extracted computational definitions have identical syntax trees to their source definitions, and the frozen detector signature matches the original default. This representative check does not substitute for a new independent benchmark. See [verification records](extraction_verification.json).
