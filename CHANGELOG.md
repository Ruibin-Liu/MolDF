# MolDf v0.7.6 (dev)

To install this version

```
pip install moldf -U
```

## New features


## Fixes

- `label_seq_id` dtype is not set to ``int`` because many cif files don't respect that specification by [RCSB](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_atom_site.label_seq_id.html).
- Check local directories for files first before attempting to download them from online resources for ``read_pdb`` and ``read_pdbx``.
