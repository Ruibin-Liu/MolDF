# MolDf v0.7.5 (dev)

To install this version

```
pip install moldf -U
```

## New features

- Added `write_jcsv` function. Any dict of Pandas Frames can be writen to the [JCSV format](https://github.com/Ruibin-Liu/JCSV).
- Added `read_jcsv` function. The `JCSV` files written by the above `write_jcsv` function can be read into dicts of DataFrames.
- The above two close [#20](https://github.com/Ruibin-Liu/MolDF/issues/20#issue-1905974450).

## Fixes

- `get_residue_template` now returns both (A, B) and (B, A) for A-B bonding information.
- For `mmcif`, columns of the `_atom_site` category but not in the provided file are skipped in building `PDBDataFrame`.
 - `label_entity_id` is char in `mmcif` specification.
 - The newline character in windows os is respected like in CSV and Pandas.
