<div align="center">
  <img src="https://raw.githubusercontent.com/Ruibin-Liu/MolDF/main/assets/logo_name.svg" width="400"><br>
</div>

-----------------

# MolDF: molecular structure data processing and analysis with DataFrame


![Tests](https://github.com/Ruibin-Liu/MolDF/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/moldf/badge/?version=latest)](https://moldf.readthedocs.io/en/latest/?badge=latest)
![Python](https://img.shields.io/badge/python-3.7-blue.svg)
[![PyPI version](https://badge.fury.io/py/moldf.svg)](https://badge.fury.io/py/moldf)

Many file formats are about different ways of integrating structured data blocks into a single file in that those blocks are related to each other in some way. The `PDBx` or `mmCIF` file format organizes structural biology data into `categories` and each category contains a structured data block which includes several `attributes`, and each attribute contains the same number of elements within a category. Those characteristics make a PDBx/mmCIF file naturally to be representable as a `Python dict` of `Pandas DataFrames`.

Our `moldf` package primarily parses a PDBx file (mmCIF file: <pdb_id>.cif) into a Python dict with PDBx category names as keys and contents belonging to the category as the corresponding values. Each category content is parsed as a Pandas DataFrame whose columns are the attribute names. On the other hand, we can write a dict of Pandas DataFrame(s) into a PDBx format in which the dict key(s) are used as category names, the DataFrame column names as attribute names, and the DataFrame row(s) as the corresponding record(s).

The old style `PDB` file format is not very well structured compared to the new PDBx file format. However, we can make `moldf` support parsing a PDB file (pdb_id.pdb) into a Python dict of Pandas DataFrames similarly, although many 'blocks' need more post processing. As such, the lines starting with `ATOM`, `HETATM`, and `TER` are read into a category named `_atom_site` which corresponds to the same category in a mmCIF file. And for NMR models, all `ATOM`, `HETATM`, and `TER` lines are read into a single DataFrame but atoms in a NMR model has the same value in the `nmr_model` column which is determined by the number in the corresponding `MODEL` line. `SEQRES` lines are read into a category named `_seq_res`.

The `TRIPOS MOL2` file format is also supported for reading using the same keyword as the PDB and PDBx files about selecting file categories. Currently, the `ATOM`, `BOND`, and `MOLECULE` are supported.

For the `PDB` and `MOL2` formats, if any other categories are required in your workflow, please raise an issue or PR.

## Requirements

Only `Pandas` is required since we need to export Pandas DataFrames. Python versions >= 3.9 in Linux, Windows, and MacOS are tested.

## Install

```bash
pip install moldf
```

## Usage examples

### 1. PDBx file

#### 1.1 If you want to read the 3D coordinates for PDB `1vii` into a Pandas DataFrame, and you have downloaded the `1vii.cif` file to your current working directory `./`, you can:

```python
from moldf import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['_atom_site'])
atoms_df = pdbx['_atom_site']
# 'atoms_df' is a Pandas DataFrame containing the '_atom_site' category which has the detailed 3D coordinates for each atom.
```

#### 1.2. If you want to read the FASTA sequence of `1vii`, you can:

```python
from moldf import read_pdbx
pdb_id = '1vii'
pdbx = read_pdbx(pdb_id=pdb_id, category_names=['_entity_poly'])
fasta_df = pdbx['_entity_poly']
fasta = fasta_df['pdbx_seq_one_letter_code_can'].to_list()[0]  # 1vii only has one sequence
# fasta == 'MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF'
```

The file content of `1VII.cif` is fetched from [RCSB](https://files.rcsb.org/view/1VII.cif) if only `pdb_id` is provided.

You can also put a Uniprot ID in the `pdb_id` keyword to get the AlphaFold2 file from [alphafold.ebi.ac.uk](alphafold.ebi.ac.uk).

By default, the fetched content will be saved to a file named `<pdb_id>.cif` under the directory given by `pdbx_file_dir` (which by default is your current working directory). You can choose not to save the file by setting `save_pdbx_file=False` in the `read_pdbx` function calling.

#### 1.3. You can read them simutanously:

```python
from moldf import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['_entity_poly', '_atom_site'])
atoms_df = pdbx['_atom_site']
fasta_df = pdbx['_entity_poly']
```

Putting a list of category names to `category_names`, you will get them if they are in the PDBx file.

#### 1.4. You can parse the whole file by using 'all':

```python
from moldf import read_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['all'])
atoms_df = pdbx['_atom_site']
fasta_df = pdbx['_entity_poly']
# and more
```

#### 1.5. Write back to a PDBx file:

```python
from moldf import read_pdbx, write_pdbx
pdbx_file = './1vii.cif'
pdbx = read_pdbx(pdbx_file, category_names=['all'])
keep = ['_atom_site', '_entity_poly']  # suppose we only want to keep the FASTA sequence and 3D coordinates.
pdbx_keep = {k: v for k, v in pdbx.items() if k in keep}
write_pdbx(pdbx_keep, '1vii_save.cif')
```
### 2. PDB file
#### 2.1. For reading the atomic information in a PDB file `1vii.pdb`:

```python
from moldf import read_pdb
pdb_file = './1vii.pdb'
pdb = read_pdb(pdb_file=pdb_file, category_names=['_atom_site'])  # We use '_atom_site' here to mirror the mmCIF format and it is the default
atoms_df = pdb['_atom_site']
# 'atoms_df' is a Pandas DataFrame containing the '_atom_site' category which has the detailed 3D coordinates for each atom.
```

#### 2.2. Suppose we only want to keep the protein residue atoms in `5u8l.pdb`:

```python
from moldf import read_pdb, write_pdb
pdb_file = './5u8l.pdb'
pdb = read_pdb(pdb_file=pdb_file, category_names=['_atom_site'])
df = pdb['_atom_site']
df = df[df.record_name == 'ATOM']
pdb['_atom_site'] = df
write_pdb(pdb, '5u8l_nohetero.pdb')
# The '5u8l_nohetero.pdb' file contains only the protein residues.
```

The `read_pdb` function can parse PDB files generated by `Chimera` by default. You can set `allow_chimera=False` in its input to fully follow the standard PDB format (although I don't see a use case).

For using a `pdb_id` or Uniprot ID (as `pdb_id`) to fetch a RCSB or AlphaFold2 PDB file, it is almost the same as the above `read_pdbx` example, except that the keywords are `save_pdb_file` and `pdb_file_dir` for saving the fetched content.

The `write_pdb` function can write PDB files that can be parsed by `Chimera` by setting `allow_chimera=True`. `allow_chimera=False` by default so that the output PDB files follow the standard PDB format strictly.

Since our package can read from and write to PDB files containing NMR models, it is straightforward to read and write trajectory files saved as PDB files by molecular dynamics software, if different frames are surrounded by pairs of `MODEL` and `ENDMDL` lines.

### 3. MOL2 file
#### For example, to read the `test.mol2` file in the `tests/test_files` folder in this repository:
```python
from moldf import read_mol2
mol2_file = './tests/test_files/test.mol2'
mol2 = read_mol2(test_file)
atoms_df = mol2['ATOM']  # The 'ATOM' category as a DataFrame.
bonds_df = mol2['BOND']  # The 'BOND' category as a DataFrame.
```
