import os
import sys

from pdbx2df.read_pdb import read_pdb

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
        test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    pdb_df = pdb["_atom_site"]
    assert (pdb_df.tail(1).nmr_model == 20).bool(), "NMR model index not read as 20."
    assert (
        pdb_df.tail(1).record_name == "TER   "
    ).bool(), "Last record name not 'TER   ' when TER lines are required."
    pdb = read_pdb(
        test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=False,
    )
    pdb_df = pdb["_atom_site"]
    assert (
        pdb_df.tail(1).record_name == "ATOM  "
    ).bool(), "Last record name not 'ATOM  ' when TER lines are not required."

    # Non-NMR PDB
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(test_file, category_names=["_atom_site"])
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
    print(pdb_df)
    assert (
        pdb_df.tail(1).atom_number == 14209
    ).bool(), "Chimera PDB didn't get larger than 9999 atoms."
    assert (
        pdb_df.tail(1).residue_name == "LEUX"
    ).bool(), "Chimera PDB didn't get 4-letter residue names."
