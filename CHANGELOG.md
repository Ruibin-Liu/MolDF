# MolDf v0.7.6 (dev)

To install this version

```
pip install moldf -U
```

## New features
- Can get a list of chain ids through `PDBDataFrame.chain_list` property, even when two chains are identical as to `chaid_id`s.
- Can read the `HET`, `HETNAM`, `HETSYN`, and `FORMUL` lines into a `_chem_comp` category now for `read_pdb`.
The `_chem_comp` category mimics the same category in the `mmCIF` format but is slightly different in that the column names are mostly taken from the `FIELD` specification in the [RCSB's doc](https://www.wwpdb.org/documentation/file-format-content/format33/sect4.html).

## Fixes

- By default, those numerical columns are not converted to numerical types because many cif files don't respect those type specifications by [RCSB](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_atom_site.label_seq_id.html).
- Check local directories for files first before attempting to download them from online resources for ``read_pdb`` and ``read_pdbx``.
- Residues in different chains but with identical (`chain_id`, `residue_name`, `residue_number`) are now listed separately in `residue_list`. A fixed file by pdbfixer if the original file is a bioassembly file containing same asym chains can have such a problem.
