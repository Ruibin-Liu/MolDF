# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Mol2 format reading.

Read a Tripos ``.mol2`` file into a dictionary of ``pandas DataFrames``. Different
categories like 'ATOM' and 'BOND' are read into different DataFrame objects.

"""
from __future__ import annotations

import os
import warnings
from collections import defaultdict

import pandas as pd  # type: ignore

IMPLEMENTED_MOL2_CATS = ["ATOM", "MOLECULE", "BOND", "HEADER"]
"""MOL2 categories that are currently implemented."""


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
"""MOL2 ``ATOM`` column names."""


BOND_COL_NAMES = (
    "bond_id",  # int
    "origin_atom_id",  # int
    "target_atom_id",  # int
    "bond_type",  # str
    "status_bit",  # str, optional
)
"""MOL2 ``BOND`` column names."""


def read_mol2(
    mol2_file: str | os.PathLike,
    category_names: list | None = None,
) -> dict[str, pd.DataFrame]:
    """Reads a ``.mol2`` file's categories into a ``dict`` of ``Pandas DataFrame`` s.

    Args:
        mol2_file(required): file name for a PDB file.
        category_names (optional): a list of categories as to the ``.mol2`` file format.
            If ``None``, [``'ATOM'``, ``'MOLECULE'``, ``'BOND'``, ``'HEADER'``] is used.
            Defaults to **None**.

    Returns:
        A dict of ``category_name`` as keys(s) and ``pd.DataFrame`` as values.

    Raises:
        NotImplementedError: if ``category_names`` not a subset of
            [``'ATOM'``, ``'MOLECULE'``, ``'BOND'``, ``'HEADER'``]
    """
    data: dict[str, pd.DataFrame] = {}
    if category_names is None:
        category_names = ["ATOM", "MOLECULE", "BOND", "HEADER"]
    for category_name in category_names:
        if category_name not in IMPLEMENTED_MOL2_CATS:
            implemented = ", ".join(IMPLEMENTED_MOL2_CATS)
            raise NotImplementedError(
                f"""Only {implemented} categories are implemented for the MOL2 format.
                Create an issue at https://github.com/Ruibin-Liu/moldf if
                you want the {category_name} category to be implemented.
                """
            )
        data[category_name] = pd.DataFrame()

    category_block_lines: dict[str, list] = defaultdict(list)
    with open(mol2_file, "r", encoding="utf-8") as mol_f:
        line = mol_f.readline()
        while line:
            if line.startswith("#"):
                category_name = "HEADER"
                line = line.lstrip("#").strip()
                if line:
                    category_block_lines[category_name].append(tuple(line.split(": ")))
            elif line.startswith("@<TRIPOS>"):
                category_name = line.strip()[9:]
                if category_name not in category_names:
                    line = mol_f.readline()
                    continue
                line = mol_f.readline()
                while line and line != "\n" and line[0] != "@":
                    category_block_lines[category_name].append(
                        tuple(line.strip().split())
                    )
                    line = mol_f.readline()
                if line and line[0] == "@":
                    continue
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
        elif category_name == "HEADER":
            data[category_name] = _get_header_df(category_block_lines[category_name])

    return data


def _get_header_df(header_lines: list[tuple]) -> pd.DataFrame:
    """Turns the ``HEADER`` lines into a ``Pandas DataFrame``.
    The ``HEADER`` lines are those starting with one or multiple ``#`` symbols.

    Args:
        header_lines (required): a list of tuples corresponding to
            each line's content. Tuples are generate by splitting the lines by ``:``.

    Returns:
        ``Pandas DataFrame`` of The ``HEADER`` category
    """
    header_attrs: dict[str, list[str]] = defaultdict(list)
    n_no_name = 0
    for header_line in header_lines:
        header_line = tuple([i.strip() for i in header_line])
        if len(header_line) == 1:
            header_attrs[f"info_{n_no_name}"] = [header_line[0]]
        elif len(header_line) == 2:
            header_attrs[header_line[0]] = [header_line[1]]
        else:
            message = f"The line {header_line} has > 2 items separated"
            message += " by ':'. moldf uses the first as column name"
            message += " and concatenate the rest as value."
            warnings.warn(
                message,
                RuntimeWarning,
                stacklevel=2,
            )
            value = ";".join(header_line[1:])
            header_attrs[header_line[0]] = [value]

    return pd.DataFrame(header_attrs)


def _get_molecule_df(molecule_lines: list[tuple]) -> pd.DataFrame:
    """Turns the ``MOLECULE`` lines into a ``Pandas DataFrame``.

    Args:
        molecule_lines (required): a list of tuples corresponding to
            each line's content.

    Returns:
        ``Pandas DataFrame`` of The ``MOLECULE`` category
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
    """Sets the data types for the ``ATOM`` category.

    Args:
        data_df (required): original ``Pandas DataFrame``
            for the ``ATOM`` category with all strings.

    Returns:
        ``Pandas DataFrame`` of The ``ATOM`` category
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
    """Sets the data types for the ``BOND`` category

    Args:
        data_df (required): original ``Pandas DataFrame``
            for the ``BOND`` category with all strings.

    Returns:
        ``Pandas DataFrame`` of The ``BOND`` category
    """
    data_df[["bond_id", "origin_atom_id", "target_atom_id"]] = data_df[
        ["bond_id", "origin_atom_id", "target_atom_id"]
    ].astype({"bond_id": "int32", "origin_atom_id": "int32", "target_atom_id": "int32"})

    return data_df
