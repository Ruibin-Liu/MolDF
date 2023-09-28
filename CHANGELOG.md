# MolDf v0.7.2

To install this version

```
pip install moldf
```

## Fixes

- The `covalent_bonds` folder is now included in the PyPI dist so that the `get_bonds_by_distance` method is usable if the package is pip-installed.

- Covalent bonds between hetero ligands and normal residues are included as well.

## Break changes

- The internal helper function `get_covalent_radii` is replaced by the `get_covalent_bond_cutoffs` function so that the code is cleaner and more maintainable.
