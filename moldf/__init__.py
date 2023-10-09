# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""
# MolDF - Super lightweight and fast mmCIF/PDB/MOL2 file parser into Pandas DataFrames and backwards writer.

### Key features:

### 1. Read a mmCIF/PDB/MOL2 file as a dictionary of Pandas DataFrame.

### 2. Atom selection using dot syntax enabled by a DataFrame subclass: PDBDataFrame

### 3. Write a dictionary of Pandas DataFrames into mmCIF/PDB files.
"""  # noqa
from moldf.constants import AMINO_ACIDS, ELEMENT_MASSES
from moldf.covalent_bond import get_covalent_bond_cutoffs, get_residue_template
from moldf.pdb_dataframe import PDBDataFrame
from moldf.read_jcsv import read_jcsv
from moldf.read_mol2 import read_mol2
from moldf.read_pdb import read_pdb
from moldf.read_pdbx import read_pdbx
from moldf.split_line import split_line
from moldf.write_jcsv import write_jcsv
from moldf.write_mol2 import write_mol2
from moldf.write_pdb import write_pdb
from moldf.write_pdbx import write_pdbx

from .version import __version__

__all__ = [
    "PDBDataFrame",
    "read_pdbx",
    "write_pdbx",
    "read_pdb",
    "write_pdb",
    "read_mol2",
    "write_mol2",
    "read_jcsv",
    "write_jcsv",
    "get_covalent_bond_cutoffs",
    "get_residue_template",
    "split_line",
    "AMINO_ACIDS",
    "ELEMENT_MASSES",
    "__version__",
]
