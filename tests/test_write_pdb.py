# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for writing PDB files."""
import os
import sys

from moldf.read_pdb import read_pdb
from moldf.write_pdb import write_pdb

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_write_pdb():
    """
    Test write_pdb function
    """
    pdb_id = "1G03"
    test_file = [CFD, "test_files", f"{pdb_id}.pdb"]
    test_file = f"{os.sep}".join(test_file)
    pdb = read_pdb(
        pdb_file=test_file, category_names=["_atom_site"], allow_chimera=False
    )

    compare_to = [CFD, "test_files", f"{pdb_id}_moldf.pdb"]
    compare_to = f"{os.sep}".join(compare_to)

    write_to = [CFD, "test_files", f"{pdb_id}_test.pdb"]
    write_to = f"{os.sep}".join(write_to)

    write_pdb(pdb, file_name=write_to)
    with open(compare_to, "r", encoding="utf-8") as cf, open(
        write_to,
        "r",
        encoding="utf-8",
    ) as wf:
        next(cf)
        next(wf)
        cf_lines = cf.readlines()
        wf_lines = wf.readlines()
        assert cf_lines == wf_lines
    os.remove(write_to)
