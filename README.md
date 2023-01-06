# pdbx2df

Parse a PDBx file (mmCIF file: pdb_id.cif) into a python dict with PDBx category names as keys and contents belonging to the category as the corresponding values. Each category content is parsed as a Pandas DataFrame whose columns are the attribute names.

## Requirements

. Pandas (>=1.0)

## Install

```
pip install pdbx2df
```

## Usage examples

1. If you want to read the 3D coordinates for PDB `1vii` into a Pandas DataFrame, and you have downloaded the `1vii.cif` file to your current working directory `./`, you can:

```python
from pdbx2df import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['_atom_site'])
atoms_df = pdbx['_atom_site']
# 'atoms_df' is a Pandas DataFrame containing the '_atom_site' category which has the detailed 3D coordinates for each atom.
```

2. If you want to read the FASTA sequence of `1vii`, you can:

```python
from pdbx2df import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['_entity_poly'])
fasta_df = pdbx['_entity_poly']
fasta = fasta_df['pdbx_seq_one_letter_code_can'].to_list()[0]  # 1vii only has one sequence
# fasta == 'MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF'
```

3. You can read them simutanously:

```python
from pdbx2df import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['_entity_poly', '_atom_site'])
atoms_df = pdbx['_atom_site']
fasta_df = pdbx['_entity_poly']
```

Putting a list of category names to `category_names`, you will get them if they are in the PDBx file.

4. You can parse the whole file by using 'all':

```python
from pdbx2df import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['all'])
atoms_df = pdbx['_atom_site']
fasta_df = pdbx['_entity_poly']
# and more
```
