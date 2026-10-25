# Example structures

`1F1T_rna.pdb` contains the supported RNA heavy atoms from PDB entry [1F1T](https://www.rcsb.org/structure/1F1T), exported with the same first-model and alternative-location rules as prediction. Ligand, solvent and other non-RNA records are omitted. The original coordinates are unchanged.

```bash
pocketseekr --rna-pdb examples/1F1T_rna.pdb --output-dir outputs/1F1T
```

`synthetic_shell.pdb` is an artificial shell used to check geometry and integration. It is not an experimental RNA structure.
