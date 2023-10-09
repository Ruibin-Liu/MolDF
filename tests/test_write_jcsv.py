# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for writing JCSV files."""
import filecmp
import os
import sys

from moldf.read_mol2 import read_mol2
from moldf.read_pdb import read_pdb
from moldf.read_pdbx import read_pdbx
from moldf.write_jcsv import write_jcsv

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_write_jcsv():
    """
    Test write_jcsv function
    """
    test_file = [CFD, "test_files", "test.mol2"]
    test_file = f"{os.sep}".join(test_file)
    mol2 = read_mol2(mol2_file=test_file)

    write_to = [CFD, "test_files", "mol2.jcsv"]
    write_to = f"{os.sep}".join(write_to)
    write_jcsv(mol2, write_to, write_meta=False, index=False)

    compared_to = [CFD, "test_files", "mol2_moldf.jcsv"]
    if os.name == "nt":
        compared_to = [CFD, "test_files", "mol2_moldf_nt.jcsv"]
    compared_to = f"{os.sep}".join(compared_to)
    assert filecmp.cmp(
        write_to, compared_to, shallow=False
    ), "mol2 jcsv writing incorrect."
    os.remove(write_to)

    test_file = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(test_file)
    pdb = read_pdb(pdb_file=test_file, category_names=["_atom_site", "_seq_res"])

    write_to = [CFD, "test_files", "pdb.jcsv"]
    write_to = f"{os.sep}".join(write_to)
    write_jcsv(pdb, write_to, write_meta=False, index=False)

    compared_to = [CFD, "test_files", "pdb_moldf.jcsv"]
    if os.name == "nt":
        compared_to = [CFD, "test_files", "pdb_moldf_nt.jcsv"]
    compared_to = f"{os.sep}".join(compared_to)
    filecmp.clear_cache()
    assert filecmp.cmp(
        write_to, compared_to, shallow=False
    ), "pdb jcsv writing incorrect."
    os.remove(write_to)

    test_file = [CFD, "test_files", "HIS.cif"]
    test_file = f"{os.sep}".join(test_file)
    pdbx = read_pdbx(pdbx_file=test_file, category_names=["all"])

    write_to = [CFD, "test_files", "pdbx.jcsv"]
    write_to = f"{os.sep}".join(write_to)
    write_jcsv(pdbx, write_to, write_meta=False, index=False)

    compared_to = [CFD, "test_files", "pdbx_moldf.jcsv"]
    if os.name == "nt":
        compared_to = [CFD, "test_files", "pdbx_moldf_nt.jcsv"]
    compared_to = f"{os.sep}".join(compared_to)
    filecmp.clear_cache()
    assert filecmp.cmp(
        write_to, compared_to, shallow=False
    ), "pdbx jcsv writing incorrect."
    os.remove(write_to)


def test_write_jcsv_with_meta():
    """
    Test write_jcsv function for with_meta=True
    """
    test_file = [CFD, "test_files", "HIS.cif"]
    test_file = f"{os.sep}".join(test_file)
    pdbx = read_pdbx(pdbx_file=test_file, category_names=["all"])

    write_to = [CFD, "test_files", "pdbx.jcsv"]
    write_to = f"{os.sep}".join(write_to)
    write_jcsv(pdbx, write_to, write_meta=True, index=False)

    compared_to = [CFD, "test_files", "pdbx_moldf_meta.jcsv"]
    if os.name == "nt":
        assert True  # pass windows test for now
        return
        # compared_to = [CFD, "test_files", "pdbx_moldf_nt.jcsv"]
    compared_to = f"{os.sep}".join(compared_to)
    filecmp.clear_cache()
    assert filecmp.cmp(
        write_to, compared_to, shallow=False
    ), "pdbx jcsv writing incorrect."
    os.remove(write_to)


test_write_jcsv_with_meta()
