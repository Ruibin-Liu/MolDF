# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for writing MOL2 files."""
import os
import sys

from moldf.read_mol2 import read_mol2
from moldf.write_mol2 import write_mol2

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_write_mol2():
    """
    Test write_mol2 function
    """
    test_file = [CFD, "test_files", "test.mol2"]
    test_file = f"{os.sep}".join(test_file)
    mol2 = read_mol2(mol2_file=test_file)

    compare_to = [CFD, "test_files", "test_moldf.mol2"]
    compare_to = f"{os.sep}".join(compare_to)

    write_to = [CFD, "test_files", "test_test.mol2"]
    write_to = f"{os.sep}".join(write_to)

    write_mol2(mol2, file_name=write_to)
    with open(compare_to, "r", encoding="utf-8") as cf, open(
        write_to,
        "r",
        encoding="utf-8",
    ) as wf:
        for _ in range(4):
            next(cf)
            next(wf)
        cf_lines = cf.readlines()
        wf_lines = wf.readlines()
        assert cf_lines == wf_lines
    os.remove(write_to)
