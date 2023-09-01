from __future__ import annotations

import os
import warnings

import pandas as pd  # type: ignore

IMPLEMENTED_MOL2_CATS = ["ATOM", "MOLECULE", "BOND"]

ATOM_COL_NAMES = (
    "atom_id",  # int
    "atom_name",  # str
    "x",  # float
    "y",  # float
    "z",  # float
    "atom_type",  # str
    "subst_id",  # int, optional
    "subst_name",  # str, optional
    "charge",  # float, optional
    "status_bit",  # str, optional
)

BOND_COL_NAMES = (
    "bond_id",  # int
    "origin_atom_id",  # int
    "target_atom_id",  # int
    "bond_type",  # str
    "status_bit",  # str, optional
)


def read_mol2(
    mol2_file: str | os.PathLike,
    category_names: list | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Read a mol2 file's categories into a dict of Pandas DataFrames.

    Args:
        mol2_file (str|os.PathLike): file name for a PDB file.
        category_names (list|None; defaults to None): a list of names for the categories as to the .mol2 file format.
            If None, ["ATOM", "MOLECULE", "BOND"] is used.

    Returns:
        dict[str, pd.DataFrame]: A dict of {category_name: pd.DataFrame of the info belongs to the category}
    """  # noqa
    data: dict[str, pd.DataFrame] = {}
    if category_names is None:
        category_names = ["ATOM", "MOLECULE", "BOND"]
    for category_name in category_names:
        if category_name not in IMPLEMENTED_MOL2_CATS:
            implemented = ", ".join(IMPLEMENTED_MOL2_CATS)
            raise NotImplementedError(
                f"""Only {implemented} categories are implemented for the MOL2 format.
                Create an issue at https://github.com/Ruibin-Liu/pdbx2df if
                you want the {category_name} category implemented.
                """
            )
        data[category_name] = pd.DataFrame()

    category_block_lines: dict[str, list] = {}
    with open(mol2_file, "r", encoding="utf-8") as mol_f:
        line = mol_f.readline()
        while line:
            if line.startswith("@<TRIPOS>"):
                category_name = line.strip()[9:]
                if category_name not in category_names:
                    line = mol_f.readline()
                    continue
                category_block_lines[category_name] = []
                line = mol_f.readline()
                while line and line != "\n":
                    category_block_lines[category_name].append(
                        tuple(line.strip().split())
                    )
                    line = mol_f.readline()
            line = mol_f.readline()

    for category_name in category_names:
        if category_name not in category_block_lines:
            warnings.warn(
                f"The required category {category_name} is not in the file {mol2_file}.",
                RuntimeWarning,
                stacklevel=2,
            )
        elif category_name == "ATOM":
            width = len(category_block_lines[category_name][0])
            data[category_name] = pd.DataFrame(
                category_block_lines[category_name],
                columns=ATOM_COL_NAMES[0:width],
            )
            data[category_name] = _set_atom_df_dtypes(data[category_name])
        elif category_name == "BOND":
            width = len(category_block_lines[category_name][0])
            data[category_name] = pd.DataFrame(
                category_block_lines[category_name],
                columns=BOND_COL_NAMES[0:width],
            )
            data[category_name] = _set_bond_df_dtypes(data[category_name])
        elif category_name == "MOLECULE":
            data[category_name] = _get_molecule_df(category_block_lines[category_name])

    return data


def _get_molecule_df(molecule_lines: list) -> pd.DataFrame:
    """Turn the 'MOLECULE' lines into a Pandas DataFrame.

    Args:
        molecule_lines (list): a list of tuples corresponding to each line's content.

    Returns:
        pd.DataFrame: the 'MOLECULE' category as a Pandas DataFrame
    """
    molecule_attrs: dict[str, list[str] | list[int]] = {}
    line_0 = {"mol_name": [" ".join(molecule_lines[0])]}
    line_1_names = ["num_atoms", "num_bonds", "num_subst", "num_feat", "num_sets"]
    line_1 = {
        name: [int(value)] for name, value in zip(line_1_names, molecule_lines[1])
    }
    line_2 = {"mol_type": [molecule_lines[2][0]]}
    line_3 = {"charge_type": [molecule_lines[3][0]]}
    molecule_attrs = {**line_0, **line_1, **line_2, **line_3}
    if len(molecule_lines) > 4:
        line_4 = {"status_bits": [molecule_lines[4][0]]}
        molecule_attrs = {**molecule_attrs, **line_4}
    if len(molecule_lines) > 5:
        line_5 = {"mol_comment": [molecule_lines[5][0]]}
        molecule_attrs = {**molecule_attrs, **line_5}

    return pd.DataFrame(molecule_attrs)


def _set_atom_df_dtypes(data_df: pd.DataFrame) -> pd.DataFrame:
    """Set the data types for the 'ATOM' category

    Args:
        data_df (pd.DataFrame): original Pandas DataFrame for the 'ATOM' category with all strings.

    Returns:
        pd.DataFrame: the 'ATOM' Pandas DataFrame dtypes corrected for 'atom_id', 'x', 'y', 'z',
            ['subst_id', ['charge']].
    """
    data_df[["atom_id", "x", "y", "z"]] = data_df[["atom_id", "x", "y", "z"]].astype(
        {"atom_id": "int32", "x": "float32", "y": "float32", "z": "float32"}
    )
    if "subst_id" in data_df.columns:
        data_df["subst_id"] = data_df["subst_id"].astype("int32")
    if "charge" in data_df.columns:
        data_df["charge"] = data_df["charge"].astype("float32")

    return data_df


def _set_bond_df_dtypes(data_df: pd.DataFrame) -> pd.DataFrame:
    """Set the data types for the 'BOND' category

    Args:
        data_df (pd.DataFrame): original Pandas DataFrame for the 'BOND' category with all strings.

    Returns:
        pd.DataFrame: dtypes corrected for 'bond_id', 'origin_atom_id', 'target_atom_id'.
    """
    data_df[["bond_id", "origin_atom_id", "target_atom_id"]] = data_df[
        ["bond_id", "origin_atom_id", "target_atom_id"]
    ].astype({"bond_id": "int32", "origin_atom_id": "int32", "target_atom_id": "int32"})

    return data_df
