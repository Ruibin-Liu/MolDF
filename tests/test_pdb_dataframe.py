# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for the PDBDataFrame class."""
import os
import shutil
import sys
from pathlib import Path

import numpy as np  # type: ignore
import pytest
from scipy.spatial.transform import Rotation  # type: ignore

from moldf.constants import ELEMENT_MASSES
from moldf.pdb_dataframe import RESIDUE_CODES, PDBDataFrame
from moldf.read_pdb import read_pdb

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_construct_from_df():
    """Testing constructing from plain DataFrame."""
    file_path = [CFD, "test_files", "1G03.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    assert isinstance(pdb_df, PDBDataFrame), "Does not return a 'PDBDataFrame instance."
    assert (
        df.columns == pdb_df.columns
    ).all(), "Columns changed when constructing from a plain DataFrame."
    assert (
        df.shape == pdb_df.shape
    ), "DataFrame shape changed when constructing from a plain DataFrame."


def test_hash_eq():
    """Testing whether PDBDataFrame is hashable and comparable."""
    file_path = [CFD, "test_files", "1G03.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    with pytest.raises(TypeError) as exception_info:
        assert df == pdb_df
    assert (
        str(exception_info.value) == "'NoneType' object is not callable"
    ), "Compared to non-directly-comparable plain DataFrame."

    pdb_df2 = pdb_df.copy()
    assert pdb_df == pdb_df2, "Can not compare with same type."


def test_property_getter():
    """Testing getting properties."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)

    assert (
        len(list(pdb_df.RESIDUE_CODES.keys())[0]) == 4
    ), "Residue name length not 4 for Chimera PDBs."

    assert (
        len(list(pdb_df.ELEMENT_MASSES.keys())[0]) == 2
    ), "Element symbol length not 2 for Chimera PDBs."

    assert pdb_df.is_chimera, "'is_chimera' not correct."

    assert pdb_df.hash_random_state == 0, "Default hash random state is not 0."

    assert (
        pdb_df.use_squared_distance
    ), "Default not using R2 for distance calculations."

    assert (
        not pdb_df.use_square_form
    ), "Default not using condensed form for distance matrix."

    atoms = pdb_df[pdb_df.record_name.isin(["ATOM  ", "HETATM"])]
    assert pdb_df.atoms.shape == atoms.shape, "TER lines not removed correctly."

    coords = pdb_df.atoms[["x_coord", "y_coord", "z_coord"]]
    assert np.allclose(
        coords.values, pdb_df.coords.values
    ), "Coords not selected correctly."

    element_set = pdb_df.element_set
    assert element_set == {"C", "O", "N", "S"}

    sequences = {
        "A": "AWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGEMGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKDPEERPTFEYLQAFLEDYFTSTEPQYQPGENL",  # noqa
        "B": "AWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGEMGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKDPEERPTFEYLQAFLEDYFTSTEPQYQPGENL",  # noqa
    }
    assert pdb_df.sequences == sequences, "Chain sequences not right."

    assert len(pdb_df.residue_list) == 520 and pdb_df.residue_list[0] == (
        "A",
        "ALA ",
        259,
    ), "Residue list not recovered."

    backbone = pdb_df[
        (pdb_df.record_name == "ATOM  ")
        & (pdb_df.atom_name.isin([" C  ", " CA ", " N  ", " O  "]))
        & (pdb_df.element_symbol.isin([" C", " N", " O"]))
    ]
    assert pdb_df.backbone == backbone, "Backbone not correct."

    side_chain = pdb_df[
        (pdb_df.record_name == "ATOM  ")
        & (~pdb_df.atom_name.isin([" C  ", " CA ", " N  ", " O  "]))
        & (pdb_df.element_symbol.isin([" C", " N", " O", " S"]))
    ]
    assert pdb_df.side_chain == side_chain, "Side chain not correct."

    ca_atoms = pdb_df[
        (pdb_df.record_name == "ATOM  ")
        & (pdb_df.atom_name.isin([" CA "]))
        & (pdb_df.element_symbol.isin([" C"]))
    ]
    assert pdb_df.ca_atoms == ca_atoms, "C alpha atoms not correct."

    heavy_atoms = pdb_df[(~pdb_df.element_symbol.isin([" H", " D", " T"]))]
    assert pdb_df.heavy_atoms == heavy_atoms, "Heavy atoms not correct."

    hetero_atoms = pdb_df[(pdb_df.record_name == "HETATM")]
    assert pdb_df.hetero_atoms == hetero_atoms, "Hetero atoms not correct."

    residues = pdb_df[(pdb_df.record_name == "ATOM  ")]
    assert pdb_df.residues == residues, "Residue atoms not correct."

    water = pdb_df[(pdb_df.record_name == "HETATM") & (pdb_df.residue_name == "HOH ")]
    assert pdb_df.water == water, "Water atoms not correct."

    assert pdb_df.n_atoms == len(pdb_df.atoms), "Number of atoms not correct."

    assert pdb_df.n_chains == 2, "Number of chains not correct."

    assert pdb_df.n_segments == 1, "Number of segments not correct."

    assert pdb_df.n_models == 1, "Number of models not correct."
    assert np.allclose(
        pdb_df.center_of_geometry, [4.1473727, 15.076957, 10.4202175]
    ), "Center of geometry not correct."

    assert np.allclose(
        pdb_df.center_of_mass, [4.1451645, 15.094206, 10.414299]
    ), "Center of mass not correct."

    assert (
        abs(pdb_df.radius_of_gyration - 741.9655564097923) <= 0.0000001
    ), "Center of mass not correct."

    assert np.allclose(
        pdb_df.atom_numbers([1, 2, 3]).distance_matrix,
        [2.194553, 6.47993342, 2.34410508],
    ), "Distance matrix not correct."


def test_property_setter():
    """Testing setting properties."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)

    pdb_df.RESIDUE_CODES = {"Ala": "A"}
    assert pdb_df.RESIDUE_CODES == {"ALA ": "A"}, "RESIDUE CODES not set correctly."
    pdb_df.RESIDUE_CODES = RESIDUE_CODES

    pdb_df.ELEMENT_MASSES = {"Ca": 1000}
    assert pdb_df.ELEMENT_MASSES == {
        "CA": 1000
    }, "ELEMENT_MASSES CODES not set correctly."
    pdb_df.ELEMENT_MASSES = ELEMENT_MASSES

    pdb_df.hash_random_state = 42
    assert pdb_df.hash_random_state == 42, "hash_random_state not set correctly."

    pdb_df.use_squared_distance = False
    assert not pdb_df.use_squared_distance, "use_squared_distance not set correctly."

    pdb_df.use_square_form = True
    assert pdb_df.use_square_form, "use_square_form not set correctly."


def test_get_bonds_by_distance():
    """Testing the get_bonds_by_distance function."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    bonds = pdb_df.residue_numbers(259).get_bonds_by_distance()
    expected = {
        (1, 2): 1,
        (2, 3): 1,
        (2, 5): 1,
        (3, 4): 2,
        (2091, 2092): 1,
        (2092, 2093): 1,
        (2092, 2095): 1,
        (2093, 2094): 2,
    }
    assert bonds == expected, "Covalent-only bonds incorrect for 5K9I 259."
    bonds = pdb_df.residue_numbers(259).get_bonds_by_distance(need_non_covalent=True)
    expected = {
        (1, 2): 1,
        (1, 3): 0.5,
        (1, 4): 0.5,
        (1, 5): 0.5,
        (2, 3): 1,
        (2, 4): 0.5,
        (2, 5): 1,
        (3, 4): 2,
        (3, 5): 0.5,
        (4, 5): 0.5,
        (2091, 2092): 1,
        (2091, 2093): 0.5,
        (2091, 2094): 0.5,
        (2091, 2095): 0.5,
        (2092, 2093): 1,
        (2092, 2094): 0.5,
        (2092, 2095): 1,
        (2093, 2094): 2,
        (2093, 2095): 0.5,
        (2094, 2095): 0.5,
    }
    assert bonds == expected, "Covalent+non-covalent bonds incorrect for 5K9I 259."


def test_get_bonds_by_template():
    """Testing the get_bonds_by_template function."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    bonds = pdb_df.residue_numbers([259, 260]).get_bonds_by_template()
    assert bonds[(1, 2)] == "SING", "N-CA bond incorrect for 5K9I"
    assert bonds[(3, 6)] == "SING", "Peptide bond incorrect for 5K9I"
    bonds = pdb_df.chain_ids(["A"]).residue_numbers([295, 601]).get_bonds_by_template()
    assert bonds[(296, 4211)] == "SING", "O44-Lys bond incorrect for 5K9I"

    pdb = read_pdb(
        pdb_id="1eru",
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
        pdb_file_dir=Path(CFD, "test_files"),
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    bonds = pdb_df.residue_numbers([32, 35]).get_bonds_by_template()
    assert bonds[(239, 240)] == "SING", "N-CA bond incorrect for 1ERU"
    assert bonds[(244, 261)] == "SING", "Disulfide bond incorrect for 1ERU"
    os.remove(Path(CFD, "test_files", "1ERU.pdb"))
    shutil.rmtree("./template_files")


test_get_bonds_by_template()


def test_get_residue_list():
    """Testing the get_residue_list class method."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    residue_list = PDBDataFrame.get_residue_list(pdb_df)

    assert len(residue_list) == 520, "get_residue_list gets length wrong."

    assert residue_list[0] == (
        "A",
        "ALA ",
        259,
    ), "get_residue_list gets first residue wrong."


def test_get_masses():
    """Testing the get_masses method."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    masses = pdb_df.atoms.get_masses()

    assert masses.shape == (4347,), "get_masses gets length wrong."

    assert np.isclose(masses[0], 14.007), "get_masses gets first atom mass wrong."


def test_get_distance_matrix():
    """Testing the get_distance_matrix class method."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)

    # Defaults: use R2 and use condensed matrix
    condensed_dis_matrix = PDBDataFrame.get_distance_matrix(pdb_df)
    assert condensed_dis_matrix.shape == (
        9446031,
    ), "get_distance_matrix gets condensed matrix shape wrong for self distance matrix."
    message = "get_distance_matrix gets 1st element wrong "
    message += "for condensed matrix for self distance matrix."
    assert np.isclose(condensed_dis_matrix[0], 2.194553), message

    # Use R and square matrix
    dis_matrix = PDBDataFrame.get_distance_matrix(
        pdb_df, use_r2=False, square_form=True
    )
    assert dis_matrix.shape == (
        4347,
        4347,
    ), "get_distance_matrix gets square matrix shape wrong for self distance matrix."
    message = "get_distance_matrix gets element [0, 1] wrong for square matrix "
    message += "while using r for self distance matrix."
    assert np.isclose(dis_matrix[0, 1], 1.48140238), message

    # Use 1 point in other_data
    dis_matrix = PDBDataFrame.get_distance_matrix(pdb_df, other_data=(0, 0, 0))
    assert dis_matrix.shape == (
        4347,
        1,
    ), "get_distance_matrix gets matrix shape wrong for distance matrix to one point."
    message = "get_distance_matrix gets element [0, 0] wrong for distance matrix "
    message += "while using r for distance matrix to one point."
    assert np.isclose(dis_matrix[0], 753.49101945), message

    # Use multiple point in other_data
    dis_matrix = PDBDataFrame.get_distance_matrix(
        pdb_df,
        other_data=((0, 0, 0), (1, 1, 1)),
    )
    assert dis_matrix.shape == (
        4347,
        2,
    ), "get_distance_matrix gets matrix shape wrong for distance matrix to two points."
    message = "get_distance_matrix gets element [0, 0] wrong for distance matrix "
    message += "while using r for distance matrix to two points."
    assert np.isclose(dis_matrix[0, 0], 753.49101945), message

    # Use another PDBDataFrame
    dis_matrix = PDBDataFrame.get_distance_matrix(
        pdb_df,
        pdb_df.ca_atoms,
    )
    message = "get_distance_matrix gets matrix shape wrong for distance matrix "
    message += "to another PDBDataFrame."
    assert dis_matrix.shape == (4347, 522), message
    message = "get_distance_matrix gets element [0, 0] wrong for distance matrix "
    message += "while using r for distance matrix to another PDBDataFrame."
    assert np.isclose(dis_matrix[0, 0], 2.194553), message


def test_rmsd():
    """Testing the rmsd method."""
    file_path = [CFD, "test_files", "1G03.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    m1 = pdb_df.nmr_models(1).coords.values
    m2 = pdb_df.nmr_models(20).coords.values
    expected_raw = np.sqrt(np.sum((m1 - m2) ** 2) / len(m1))
    assert np.isclose(
        expected_raw, pdb_df.rmsd(align=False)[-1]
    ), "Unaligned RMSD between model 1 and 20 are different."

    m1 = m1 - np.mean(m1, axis=0)
    m2 = m2 - np.mean(m2, axis=0)
    rot = Rotation.align_vectors(m1, m2)[0]
    m2 = rot.apply(m2)
    expected_aligned = np.sqrt(np.sum((m1 - m2) ** 2) / len(m1))
    assert np.isclose(
        expected_aligned, pdb_df.rmsd(align=True)[-1]
    ), "Aligned RMSD between model 1 and 20 are different."


def test_filter_num_col():
    """Testing the _filter_num_col method."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)

    # One int
    sub_df = pdb_df._filter_num_col(1, "atom_number")
    expected = pdb_df[pdb_df.atom_number.isin([1])]
    assert sub_df == expected, "Filter by 1 integer failed with default keywords."

    # One float
    sub_df = pdb_df._filter_num_col(1, "x_coord")
    expected = pdb_df[np.abs(pdb_df.x_coord.values - 1) <= 0.01]
    assert sub_df == expected, "Filter by 1 float failed with default keywords."

    # A list of ints
    sub_df = pdb_df._filter_num_col([1, 2], "atom_number")
    expected = pdb_df[pdb_df.atom_number.isin([1, 2])]
    assert (
        sub_df == expected
    ), "Filter by a list of two integer failed with default keywords."

    # Mixed types in list
    except_message = "Only 'int' is allowed in 'value' if it is a list, "
    except_message += "but <class 'float'> was in [1, 2.0]."
    with pytest.raises(ValueError) as exception_info:
        sub_df = pdb_df._filter_num_col([1, 2.0], "atom_number")
    assert (
        str(exception_info.value) == except_message
    ), "Did not catch the error of using mixed types in list."

    # One str
    except_message = "Only 'int', 'float', or 'list[int]' are allowed in 'value' "
    except_message += "but <class 'str'> was put."
    with pytest.raises(ValueError) as exception_info:
        sub_df = pdb_df._filter_num_col("1", "atom_number")
    assert (
        str(exception_info.value) == except_message
    ), "Did not catch the error of using a str type."

    # Allowed columns
    except_message = "Only '['atom_number', 'residue_number', 'x_coord', 'y_coord', "
    except_message += "'z_coord', 'occupancy', 'b_factor', 'nmr_model']' are allowed "
    except_message += "in 'num_col_name' but atom_name was put."
    with pytest.raises(ValueError) as exception_info:
        sub_df = pdb_df._filter_num_col(1, "atom_name")
    assert (
        str(exception_info.value) == except_message
    ), "Did not catch the error of using a wrong column name."

    # Wrong relations
    except_message = "Only '['<', '<=', '=', '>=', '>']' are allowed "
    except_message += "in 'relation' but ? was put."
    with pytest.raises(ValueError) as exception_info:
        sub_df = pdb_df._filter_num_col(1, "atom_number", relation="?")
    assert (
        str(exception_info.value) == except_message
    ), "Did not catch the error of using a wrong 'relation'."

    # relation for list
    warning_message = "'relation' is ignored when a list is provided to 'value'."
    with pytest.warns(RuntimeWarning) as warning_info:
        sub_df = pdb_df._filter_num_col([1, 2], "atom_number", relation="?")
    assert (
        warning_info[0].message.args[0] == warning_message
    ), "Did not catch the warning of using 'relation' when a list is given."

    # relation and invert
    sub_df_1 = pdb_df._filter_num_col(2, "atom_number", relation="<=", invert=False)
    sub_df_2 = pdb_df._filter_num_col(2, "atom_number", relation=">", invert=True)
    assert (
        sub_df_1 == sub_df_2
    ), "<= and invert=False should have been the same as > and invert=True"


def test_atom_selections():
    """Testing atom selection functions."""
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)

    # record_name
    sub_df = pdb_df.record_names(["HETATM"])
    expected = pdb_df[pdb_df.record_name == "HETATM"]
    assert sub_df == expected, "record_names selection failed."

    # atom_number
    sub_df = pdb_df.atom_numbers([1, 2, 3])
    expected = pdb_df[pdb_df.atom_number.isin([1, 2, 3])]
    assert sub_df == expected, "atom_numbers selection failed."

    # atom_name
    sub_df = pdb_df.atom_names(["CA", "C"])
    expected = pdb_df[pdb_df.atom_name.isin([" CA ", " C  "])]
    assert sub_df == expected, "atom_names selection failed."

    # alt_loc
    sub_df = pdb_df.alt_locs(["B"])
    expected = pdb_df[pdb_df.alt_loc.isin(["B", " "])]
    assert (
        sub_df.alt_loc.unique() == expected.alt_loc.unique()
    ).all(), "alt_locs selection failed."

    # residue_name
    sub_df = pdb_df.residue_names(["ALA"])
    expected = pdb_df[pdb_df.residue_name.isin(["ALA "])]
    assert sub_df == expected, "residue_names selection failed."

    # chain_id
    sub_df = pdb_df.chain_ids(["A"])
    expected = pdb_df[pdb_df.chain_id.isin(["A"])]
    assert sub_df == expected, "chain_ids selection failed."

    # residue_number
    sub_df = pdb_df.residue_numbers([300, 301])
    expected = pdb_df[pdb_df.residue_number.isin([300, 301])]
    assert sub_df == expected, "residue_numbers selection failed."

    # insertion
    sub_df = pdb_df.insertions([" "])
    expected = pdb_df[pdb_df.insertion.isin([" "])]
    assert sub_df == expected, "insertions selection failed with insertion=' '."
    assert (
        len(pdb_df.insertions([""])) == 0
    ), "insertions selection failed with insertion=''."

    # x_coord
    sub_df = pdb_df.x_coords(1.0)
    expected = pdb_df[pdb_df.x_coord.values <= 1.0]
    assert sub_df == expected, "x_coords selection failed."

    # y_coord
    sub_df = pdb_df.y_coords(1.0, relation=">")
    expected = pdb_df[pdb_df.y_coord.values > 1.0]
    assert sub_df == expected, "y_coords selection failed."

    # z_coord
    sub_df = pdb_df.z_coords(1.0, relation="=", epsilon=1.0)
    expected = pdb_df[np.abs(pdb_df.z_coord.values - 1.0) <= 1.0]
    assert sub_df == expected, "z_coords selection failed."

    # occupancy
    sub_df = pdb_df.occupancies(1.0, relation="=", epsilon=0.0001)
    expected = pdb_df[np.abs(pdb_df.occupancy.values - 1.0) <= 0.0001]
    assert sub_df == expected, "occupancies selection failed."

    # b_factor
    sub_df = pdb_df.b_factors(50, relation=">=")
    expected = pdb_df[pdb_df.b_factor.values >= 50]
    assert sub_df == expected, "b_factors selection failed."

    # segment_id
    sub_df = pdb_df.segment_ids(["    "])
    expected = pdb_df
    assert sub_df == expected, "segment_ids selection failed."

    # element_symbol
    sub_df = pdb_df.element_symbols(["O", "S"])
    expected = pdb_df[pdb_df.element_symbol.isin([" O", " S"])]
    assert sub_df == expected, "element_symbols selection failed."

    # charge
    sub_df = pdb_df.charges(["  "])
    expected = pdb_df
    assert sub_df == expected, "charges selection failed."

    # distances
    first_atom = pdb_df.head(1)
    sub_df = pdb_df.distances(first_atom, cut_off=1.49, to="all")
    expected = pdb_df.atom_numbers([1, 2])
    assert (
        sub_df == expected
    ), "distances selection failed with a PDBDataFrame instance."

    sub_df = pdb_df.distances([-1.972, 3.107, -27.202], cut_off=1.49, to="all")
    expected = pdb_df.atom_numbers([1, 2])
    assert sub_df == expected, "distances selection failed with one point."

    points = [
        [-1.972, 3.107, -27.202],
        [-0.601, 3.663, -27.278],
    ]
    sub_df = pdb_df.distances(points, cut_off=1.49, to="any")
    expected = pdb_df.atom_numbers([1, 2])
    assert sub_df == expected, "distances selection failed with two points."

    sub_df_1 = pdb_df.distances(points, cut_off=0, to="any", invert=True)
    sub_df_2 = pdb_df.distances(points, cut_off=0, to="all", invert=False)
    assert (
        sub_df_2 == sub_df_1
    ), "to='any' and invert=True not the same as to='all' and invert=False."


def test_nmr_models():
    """Testing atom selection using nmr_model."""
    file_path = [CFD, "test_files", "1G03.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)

    # nmr_model
    sub_df = pdb_df.nmr_models([1, 2, 20])
    expected = pdb_df.atoms[pdb_df.atoms.nmr_model.isin([1, 2, 20])]
    assert sub_df == expected, "nmr_models selection failed."
    assert (
        sub_df.nmr_model.unique() == [1, 2, 20]
    ).all(), "Not all 1,2&20 models are selected."

    # Without nmr_model column
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)
    pdb = read_pdb(
        pdb_file=test_file,
        category_names=["_atom_site"],
        allow_chimera=True,
        need_ter_lines=True,
    )
    df = pdb["_atom_site"]
    pdb_df = PDBDataFrame(df)
    sub_df = pdb_df.nmr_models([1, 2, 20])
    expected = pdb_df
    assert sub_df == expected, "nmr_models selection failed when there is only 1 model."
