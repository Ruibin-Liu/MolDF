# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for reading PDB files."""
import os
import sys

from moldf.read_pdb import read_pdb

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_read_pdb():
    """
    Test read_pdb function
    """
    # NMR PDB
    file_path = [CFD, "test_files", "1G03.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    pdb_df = pdb["_atom_site"]
    assert (pdb_df.tail(1).nmr_model == 20).item(), "NMR model index not read as 20."
    assert (
        pdb_df.tail(1).record_name == "TER   "
    ).item(), "Last record name not 'TER   ' when TER lines are required."
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=False,
    )
    pdb_df = pdb["_atom_site"]
    assert (
        pdb_df.tail(1).record_name == "ATOM  "
    ).item(), "Last record name not 'ATOM  ' when TER lines are not required."

    # From PDB ID
    pdb_id = "1VIi"
    pdb = read_pdb(
        pdb_id=pdb_id,
        category_names=["_atom_site"],
        save_pdb_file=True,
        allow_chimera=False,
    )
    pdb_df = pdb["_atom_site"]
    assert os.path.exists(
        f"./PDB_files/{pdb_id.upper()}.pdb"
    ), "File not saved to PDB_files if pdb_file_dir not provided."
    assert (
        pdb_df.head(1).residue_name == "MET"
    ).item(), "First 1VII residue is not MET."

    # From Uniprot ID
    uniprot_id = "P01116"
    pdb = read_pdb(
        pdb_id=uniprot_id,
        category_names=["_atom_site"],
        save_pdb_file=True,
        allow_chimera=False,
    )
    pdb_df = pdb["_atom_site"]
    assert os.path.exists(
        f"./PDB_files/{uniprot_id.upper()}.pdb"
    ), "File not saved to PDB_files if pdb_file_dir not provided."
    assert (
        pdb_df.head(1).residue_number == 1
    ).item(), "First P01116 residue number is not 1."

    # Non-NMR PDB
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(pdb_file=test_file, category_names=["_atom_site"])
    pdb_df = pdb["_atom_site"]
    assert (
        "nmr_model" not in pdb_df.columns
    ), "nmr_model column appears when the file is not NMR."

    # Chimera
    file_path = [CFD, "test_files", "5K9I_mimic_Chimera.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    pdb_df = pdb["_atom_site"]
    pdb_df = pdb_df[pdb_df.record_name == "TER   "]
    assert (
        pdb_df.tail(1).atom_number == 14209
    ).item(), "Chimera PDB didn't get larger than 9999 atoms."
    assert (
        pdb_df.tail(1).residue_name == "LEUX"
    ).item(), "Chimera PDB didn't get 4-letter residue names."


def test_read_seq_res():
    """Test reading the 'SEQRES' section."""
    pdb_id = "1VII"
    pdb = read_pdb(
        pdb_id=pdb_id,
        category_names=["_seq_res"],
        save_pdb_file=True,
        allow_chimera=False,
    )
    seq_df = pdb["_seq_res"]
    assert seq_df["chain_sequence"].to_list() == [
        "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF"
    ], "1VII sequence not read in correctly."
