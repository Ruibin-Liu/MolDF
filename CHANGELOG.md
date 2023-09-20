# pdbx2df v0.6.6 (in-dev)

From this version on, changelog is kept for each new version, compared to the previous version.

## New features
- Parsing `mol2` file header lines starting with `#` in `read_mol2`.
- `radius_of_gyration` can be accessed as a property in `PDBDataFrame`.
- `rmsd` method added in `PDBDataFrame`.
- Can feed `PDBDataFrame` with the PDBx `_atom_site` DataFrame now.

## Fixes
- Corrected `center_of_mass` calculation in `PDBDataFrame`.

## Docs
- Cleaned up empty lines in all function and method's docstrings.
