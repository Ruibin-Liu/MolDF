# pdbx2df
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/pdbx2df
"""
# pdbx2df - Super lightweight and fast mmCIF/PDB/MOL2 file parser into Pandas DataFrames and backwards writer.

### Key features:

### 1. Read a mmCIF/PDB/MOL2 file as a dictionary of Pandas DataFrame.

### 2. Atom selection using dot syntax enabled by a DataFrame subclass: PDBDataFrame

### 3. Write a dictionary of Pandas DataFrames into mmCIF/PDB files.
"""  # noqa
from pdbx2df.constants import AMINO_ACIDS, ELEMENT_MASSES
from pdbx2df.pdb_dataframe import PDBDataFrame
from pdbx2df.read_mol2 import read_mol2
from pdbx2df.read_pdb import read_pdb
from pdbx2df.read_pdbx import read_pdbx
from pdbx2df.split_line import split_line
from pdbx2df.write_mol2 import write_mol2
from pdbx2df.write_pdb import write_pdb
from pdbx2df.write_pdbx import write_pdbx

from .version import __version__

__all__ = [
    "PDBDataFrame",
    "read_pdbx",
    "write_pdbx",
    "read_pdb",
    "write_pdb",
    "read_mol2",
    "write_mol2",
    "split_line",
    "AMINO_ACIDS",
    "ELEMENT_MASSES",
    "__version__",
]
