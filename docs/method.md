# Method

PocketSeekR takes observed RNA heavy-atom coordinates and produces up to five ranked pocket centers with RNA context crops. The fixed default method searches at most 40 candidate regions before final selection. All prediction steps use RNA geometry and chemical typing; ligand coordinates and trained model weights are not inputs.

## Accessible regions

A deterministic coordinate frame is constructed from the RNA coordinates. A 1.5 Å grid is searched for points with 1.2–4.0 Å surface clearance and at least three nearby atoms. Twenty-six ray directions describe local enclosure and opposing support. Six-face grid neighbors are connected only when the entire segment is clear of the RNA atoms expanded by a 1.2 Å probe radius. Connected regions grow within 6 Å of a seed and above half its geometric score, with a minimum of three voxels.

The regional geometric score combines the maximum voxel score with a saturating volume factor. A representative center is the member voxel nearest the geometry-weighted mean. Regions are filtered with a 4 Å center-separation rule and RNA crop overlap, retaining at most 40. Context crops extend 10 Å from each center and include complete observed nucleotides, with at least 12 RNA heavy atoms.

## Regional chemical support

Aromatic faces and polar directions are derived from observed RNA base coordinates. Ring planes require complete observed rings with plane-fit RMS deviation at most 0.25 Å. Directions and visibility use local bond geometry and atom occlusion. Missing structural information is not replaced with invented coordinates.

Aromatic and polar potentials are transformed into bounded fields `F` and `P`. A row-normalized Gaussian kernel `K` connects neighbors within 1.5 Å whose full probe segment is accessible, including an identity self-term. The descriptor is evaluated on every original region voxel:

```text
F_i = 1 - exp(-U_i)
P_i = 1 - exp(-H_i)
B_i = 0.75 F_i + 0.25 P_i
C_i = [sqrt(F_i (K P)_i) + sqrt(P_i (K F)_i)] / 2
c(r) = 0.5 mean(B_i) + 0.5 mean(C_i)
Q(r) = G(r) [1 + 2 c(r)]
```

The default representative centers remain determined by geometry. The bounded chemical descriptor measures regional support and local cross-support; it is not a force-field energy or a ligand-specific affinity prediction.

## Joint selection

Let `Jmax(r,S)` be the maximum Jaccard overlap between the RNA atom crop of a candidate and those already selected. Candidates are chosen greedily by:

```text
gain(r | S) = Q(r) [1 - Jmax(r,S)]^0.25
```

Candidates within 3 Å of an already selected center, or with crop Jaccard overlap at least 0.8, are excluded. Ties use fixed candidate identifiers. The output preserves the actual selection order and contains at most five pockets; searches with fewer candidates keep their actual size.

## Implementation

The installed defaults are defined by `StageAConfig` and `PocketSettings`; the matching YAML is included in the wheel. The original computational definitions and their source hashes are listed in [source_extraction.json](source_extraction.json). No model fitting occurs during prediction.
