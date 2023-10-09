# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for reading jcsv files."""
import os
import sys

from moldf.read_jcsv import read_jcsv

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_read_jcsv_nometa():
    """
    Test read_jcsv function without meta in jcsv
    """
    # without meta
    file_path = [CFD, "test_files", "pdbx_moldf.jcsv"]
    test_file = f"{os.sep}".join(file_path)
    jcsv = read_jcsv(test_file)
    cat_size = {
        "_chem_comp": (1, 24),
        "_chem_comp_atom": (21, 18),
        "_chem_comp_bond": (21, 7),
        "_pdbx_chem_comp_descriptor": (7, 5),
        "_pdbx_chem_comp_identifier": (2, 5),
        "_pdbx_chem_comp_audit": (2, 4),
    }
    assert len(jcsv) == 6, "Not all categories are read."
    for k, v in jcsv.items():
        assert cat_size[k] == v.shape, f"Category {k} read incorrectly"


def test_read_jcsv_meta():
    """
    Test read_jcsv function with meta in jcsv
    """
    file_path = [CFD, "test_files", "pdbx_moldf_meta.jcsv"]
    test_file = f"{os.sep}".join(file_path)
    jcsv = read_jcsv(test_file)
    cat_size = {
        "_chem_comp": (1, 24),
        "_chem_comp_atom": (21, 18),
        "_chem_comp_bond": (21, 7),
        "_pdbx_chem_comp_descriptor": (7, 5),
        "_pdbx_chem_comp_identifier": (2, 5),
        "_pdbx_chem_comp_audit": (2, 4),
    }
    assert len(jcsv) == 6, "Not all categories are read."
    for k, v in jcsv.items():
        assert cat_size[k] == v.shape, f"Category {k} read incorrectly"
