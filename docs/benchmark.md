# Benchmark scope

The developmental comparison uses the fixed HARIBOSS snapshot and structural processing documented by the source project: RNA-only ligand contacts, `our_interest=1`, X-ray/cryo-EM structures at resolution ≤5 Å, and the first model of NMR structures after the existing quality filters. It comprises 940 ligand instances, 554 PDB entries and 553 independent RNA inputs.

## Localization results

| Method | Top1 | Top3 | Top5 | Top5 contact coverage ≥80% |
|---|---:|---:|---:|---:|
| Adopted PocketSeekR source method | 314 (33.4%) | 565 (60.1%) | 652 (69.4%) | 647 (68.8%) |
| Fixed fpocketR localization core | 238 (25.3%) | 387 (41.2%) | 428 (45.5%) | 397 (42.2%) |
| GHECOM, fixed multiscale configuration | 102 (10.9%) | 207 (22.0%) | 211 (22.4%) | 280 (29.8%) |
| RBind formula reproduction + spatial adapter | 49 (5.2%) | 101 (10.7%) | 124 (13.2%) | 236 (25.1%) |
| Official Rsite numeric algorithm + spatial adapter | 50 (5.3%) | 158 (16.8%) | 217 (23.1%) | 326 (34.7%) |

The denominator is 940 for each method, including instances with no candidate. A localization hit requires at least one of the first k centers to be within 4.5 Å of any observed native ligand heavy atom. Coverage is measured against RNA atoms within 4.5 Å of the observed ligand, using each center's 10 Å whole-observed-nucleotide crop.

These numbers describe the original developmental cohort, which was used during method development and parameter calibration. They are not results from a new independent test, and this export has not rerun all 940 instances. Standalone extraction equivalence is checked separately on representative structures.

## Comparator boundaries

fpocketR uses the fixed fpocket 4.0.3 localization core (`-m 3 -M 5.7 -i 42 -D 1.65 -A 3 -p 0`), not the complete PyMOL/PyProDy wrapper. GHECOM uses one fixed multiscale configuration. RBind is a formula reproduction checked against an official example. Rsite uses the original R script's numeric output.

RBind/Rsite predict nucleotide sites, so a shared additional spatial adapter groups positive nucleotide heavy-atom centers with edges at ≤8 Å, takes a mean regional center and ranks by the highest node score. This adapter is not the authors' native pocket-generation output. Its effect is part of the reported localization performance. See the [archived comparison protocol](source_classical_pocket_comparison.md).

The full experimental caches, downloaded databases and comparator binaries are retained in the source workspace and are not part of this small prediction repository.
