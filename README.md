# pdbx2df

Parse a PDBx file (mmCIF file: pdb_id.cif) into a python dict with PDBx category names as keys and contents belonging to the category as the corresponding values. Each category content is parsed as a Pandas DataFrame whose columns are the attribute names. On the other hand, we can write a dict of Pandas DataFrame(s) into a PDBx format in which the dict key(s) are used as category names, the DataFrame column names as attribute names, and the DataFrame row(s) as the corresponding record(s).

Also supports parsing a PDB file (pdb_id.pdb) into a python dict of Pandas DataFrames. Currently, only the lines starting with 'ATOM', 'HETATM', and 'TER' are read into a category named '_atom_site' which corresponds to the same category in a mmCIF file.

## Requirements

- Pandas (>=1.0)

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

5. Write back to a PDBx file:

```python
from pdbx2df import read_pdbx, write_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['all'])
keep = ['_atom_site', '_entity_poly']  # suppose we only want to keep the FASTA sequence and 3D coordinates.
pdbx_keep = {k: v for k, v in pdbx.items() if k in keep}
write_pdbx(pdbx_keep, '1vii_save.cif')
```

6. For reading the atomic information in a PDB file `1vii.pdb`:

```python
from pdbx2df import read_pdb
pdb_file = './1vii.pdb'
pdb = read_pdb(pdb_file, category_names=['_atom_site'])  # We use '_atom_site' here to mirror the mmCIF format
atoms_df = pdb['_atom_site']
# 'atoms_df' is a Pandas DataFrame containing the '_atom_site' category which has the detailed 3D coordinates for each atom.
```

7. Suppose we only want to keep the residue atoms in `5u8l.pdb`:

```python
from pdbx2df import read_pdb, write_pdb
pdb_file = './5u8l.pdb'
pdb = read_pdb(pdb_file, category_names=['_atom_site'])
df = pdb['_atom_site']
df = df[df.record_name == 'ATOM']
pdb['_atom_site'] = df
write_pdb(pdb, '5u8l_nohetero.pdb')
# The '5u8l_nohetero.pdb' file contains only the protein residues.
```
