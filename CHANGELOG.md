# MolDf v0.7.5 (dev)

To install this version

```
pip install moldf -U
```

## New features

- Added `write_jcsv` function. Any dict of Pandas Frames can be writen to the [JCSV format](https://github.com/Ruibin-Liu/JCSV).

## Fixes

- `get_residue_template` now returns both (A, B) and (B, A) for A-B bonding information.
