# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for reading mol2 files."""
import os
import sys

import pytest

from moldf.read_mol2 import read_mol2

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_read_mol2():
    """
    Test read_mol2 function
    """
    # Correct ones
    file_path = [CFD, "test_files", "test.mol2"]
    test_file = f"{os.sep}".join(file_path)
    mol2 = read_mol2(test_file)
    mol2_atom_df = mol2["ATOM"]
    mol2_bond_df = mol2["BOND"]
    mol2_molecule_df = mol2["MOLECULE"]
    assert (mol2_atom_df.head(1).atom_id == 1).item(), "The first atom_id is not 1."
    assert (mol2_bond_df.tail(1).bond_id == 62).item(), "The last bond_id is not 62."
    assert (
        mol2_molecule_df.tail(1).mol_name == "KEGG_DURG-00000915-01"
    ).item(), "The MOLECULE mol_name is not KEGG_DURG-00000915-01."

    # Not implemented
    message = """Only ATOM, MOLECULE, BOND, HEADER categories are implemented for the MOL2 format.
                Create an issue at https://github.com/Ruibin-Liu/moldf if
                you want the ANCHOR_ATOM category to be implemented.
                """  # noqa
    with pytest.raises(NotImplementedError) as error_msg:
        _ = read_mol2(test_file, category_names=["ANCHOR_ATOM"])
    assert message in str(error_msg.value), "NotImplementedError not raised correctly."
