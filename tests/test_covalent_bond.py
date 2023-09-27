# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for reading mol2 files."""
import os
import shutil
import sys

import pytest

from moldf.covalent_bond import get_covalent_radii, get_residue_template

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_get_covalent_radii():
    """Tests for the 'get_covalent_radii' function."""

    covalent_radii = get_covalent_radii(["C", "H", "O", "N", "S"])
    assert covalent_radii.shape == (5, 5), "C/H/O/N/S covalent radii read wrongly."
    assert covalent_radii.single_C.to_list() == [
        31,
        76,
        71,
        66,
        105,
    ], "C/H/O/N/S single covalent radii read wrongly for 'single_C'."

    covalent_radii = get_covalent_radii(
        ["C", "H", "O", "N", "S"], single_radii_set="single_PA"
    )
    assert covalent_radii.single_PA.to_list() == [
        32,
        75,
        71,
        63,
        103,
    ], "C/H/O/N/S single covalent radii read wrongly for 'single_PA'."
    with pytest.raises(AttributeError) as exception_info:
        covalent_radii.single_C
    assert (
        str(exception_info.value) == "'DataFrame' object has no attribute 'single_C'"
    ), "Column 'single_C' not dropped when 'single_PA' is used."
    assert covalent_radii.double.isna()[0], "'H' double radius is not NaN."


def test_get_residue_template():
    """Tests for the 'get_residue_template' function."""
    file_path = [CFD, "test_files", "HIS.cif"]
    test_file = f"{os.sep}".join(file_path)
    his_bonds = get_residue_template(
        residue_name="HIS", residue_template_file=test_file
    )
    expected = {
        ("N", "CA"): ("SING", True, "N"),
        ("N", "H"): ("SING", True, "N"),
        ("N", "H2"): ("SING", True, "N"),
        ("CA", "C"): ("SING", True, "N"),
        ("CA", "CB"): ("SING", True, "N"),
        ("CA", "HA"): ("SING", True, "N"),
        ("C", "O"): ("DOUB", True, "N"),
        ("C", "OXT"): ("SING", True, "N"),
        ("CB", "CG"): ("SING", True, "N"),
        ("CB", "HB2"): ("SING", True, "N"),
        ("CB", "HB3"): ("SING", True, "N"),
        ("CG", "ND1"): ("SING", False, "N"),
        ("CG", "CD2"): ("DOUB", False, "N"),
        ("ND1", "CE1"): ("DOUB", False, "N"),
        ("ND1", "HD1"): ("SING", True, "N"),
        ("CD2", "NE2"): ("SING", False, "N"),
        ("CD2", "HD2"): ("SING", True, "N"),
        ("CE1", "NE2"): ("SING", False, "N"),
        ("CE1", "HE1"): ("SING", True, "N"),
        ("NE2", "HE2"): ("SING", True, "N"),
        ("OXT", "HXT"): ("SING", True, "N"),
    }
    assert his_bonds == expected, "HIS bonds read wrongly."

    o44_bonds = get_residue_template(residue_name="O44")
    assert o44_bonds[("CAM", "CAL")] == (
        "SING",
        True,
        "N",
    ), "O44 bond information incorrect."
    shutil.rmtree("./template_files")
