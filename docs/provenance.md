# Source provenance

PocketSeekR's computational code was extracted on 2026-10-25 from the RNA-only prediction path in [RiboPoseDiff](https://github.com/wangleiofficial/RiboPoseDiff), source commit `26d6c5715b6792726b792dbe962f08044b724021`. The extraction includes only the dependencies reached by the current regional geometry, chemical scoring, pocket selection and RNA parsing functions. It does not copy docking, training, browser interfaces or historical evaluation outputs.

[source_extraction.json](source_extraction.json) lists the exact source SHA256 values and selected symbols for each computational module. The extracted function and class definitions preserve their source syntax trees. The standalone geometry configuration and detector signature match `configs/pocketseekr.yaml` in the source repository; `configs/default.yaml` matches the adopted chemical and selection settings.

The public API, CLI, package metadata and installation workflow are standalone wrappers. All prediction defaults are installed with the Python package and are independent of repository-relative paths. The internal name `calibration` denotes the source parameter containers; prediction does not fit a model or update parameters.

The full candidate regions and final selected pockets were compared against the source implementation on PDB structures 1F1T, 1NTA, 2GDI and 8K7W. The per-structure checks and computational-definition checks are recorded in [extraction_verification.json](extraction_verification.json).

The PDB fixtures are structural inputs for chemistry and integration tests. `examples/1F1T_rna.pdb` contains only the supported RNA heavy atoms of PDB [1F1T](https://www.rcsb.org/structure/1F1T), exported using the prediction parser. `examples/synthetic_shell.pdb` is an artificial geometric test input.
