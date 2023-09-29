# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for reading mol2 files."""
import os
import shutil
import sys

from moldf.covalent_bond import get_covalent_bond_cutoffs, get_residue_template

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_get_covalent_bond_cutoffs():
    """Tests for the 'get_covalent_bond_cutoffs' function."""

    single_bonds, double_bonds, triple_bonds = get_covalent_bond_cutoffs(
        ["C", "H", "O", "N", "S"]
    )
    assert (
        abs(single_bonds[(" O", " O")] - 2.3104) < 0.00001
    ), "O-O single bond incorrect."
    assert (
        abs(double_bonds[(" C", " C")] - 1.9321) < 0.00001
    ), "C-C double bond incorrect."
    assert (
        abs(triple_bonds[(" N", " C")] - 1.4161) < 0.00001
    ), "N-C triple bond incorrect."


def test_get_residue_template():
    """Tests for the 'get_residue_template' function."""
    file_path = [CFD, "test_files", "HIS.cif"]
    test_file = f"{os.sep}".join(file_path)
    his_bonds = get_residue_template(
        residue_name="HIS", residue_template_file=test_file
    )
    assert his_bonds[("OXT", "HXT")] == ("SING", True, "N"), "HIS bonds read wrongly."
    assert (
        his_bonds[("HXT", "OXT")] == his_bonds[("OXT", "HXT")]
    ), "HIS bonds info not symetrical"

    o44_bonds = get_residue_template(residue_name="O44")
    assert o44_bonds[("CAM", "CAL")] == (
        "SING",
        True,
        "N",
    ), "O44 bond information incorrect."
    shutil.rmtree("./template_files")
