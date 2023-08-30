from __future__ import annotations

import os
import warnings

import pandas as pd  # type: ignore

IMPLEMENTED_PDB_CATS = ["_atom_site"]
ATOM_SITE = ("ATOM", "HETATM", "TER")
NMR_MDL = ("NUMMDL", "MODEL", "ENDMDL")


def read_pdb(
    pdb_file: str | os.PathLike,
    category_names: list | None = None,
    allow_chimera: bool = True,
) -> dict:
    """
    Read a pdb file categories into Pandas DataFrame.

    Args:
        pdb_file (str|os.PathLike): file name for a PDB file.
        category_names (list|None; defaults to None): a list of names for the categories as to the mmCIF file format.
            If None, "_atom_site" is used.
            To be consistent with the PDBx file format, the following category names are used to refer
            to block(s) in a PDB file and only they are supported:
            1. _atom_site: 'ATOM', 'HETATM', and 'TER' lines
            2. TBD
        allow_chimera (bool; defaults to True): whether to allow Chimera-formatted PDB files.

    Returns:
        A dict of {category_name: pd.DataFrame of the info belongs to the category}
    """  # noqa
    data: dict[str, pd.DataFrame] = {}
    if not category_names:
        category_names = ["_atom_site"]
    for category_name in category_names:
        if category_name not in IMPLEMENTED_PDB_CATS:
            implemented = ", ".join(IMPLEMENTED_PDB_CATS)
            raise ValueError(f"Only {implemented} are implented for the PDB format.")
        data[category_name] = pd.DataFrame()

    num_nmr_models = 0
    n_model_lines = 0
    nmr_model = None
    n_endmdl_lines = 0
    atom_site_rows = []

    with open(pdb_file, "r", encoding="utf-8") as pf:
        for line in pf:
            if line.startswith(NMR_MDL):
                if line.startswith("NUMMDL"):
                    num_nmr_models = int(line[6:].strip())
                elif line.startswith("MODEL"):
                    n_model_lines += 1
                    nmr_model = int(line[6:].strip())
                elif line.startswith("ENDMDL"):
                    n_endmdl_lines += 1
            elif line.startswith(ATOM_SITE):
                atom_site_rows.append(
                    _split_atom_line(
                        line, nmr_model=nmr_model, allow_chimera=allow_chimera
                    )
                )
    if num_nmr_models != n_model_lines or num_nmr_models != n_endmdl_lines:
        warnings.warn(
            f"The NUMMDL says {num_nmr_models} NMR models, but only {n_model_lines} found.",
            RuntimeWarning,
            stacklevel=2,
        )
    column_names = [
        "record_name",
        "atom_number",
        "atom_name",
        "alt_loc",
        "residue_name",
        "chain_id",
        "residue_number",
        "insertion",
        "x_coord",
        "y_coord",
        "z_coord",
        "occupancy",
        "b_factor",
        "segment_id",
        "element_symbol",
        "charge",
        "nmr_model",
    ]
    data["_atom_site"] = pd.DataFrame(atom_site_rows, columns=column_names)
    if n_model_lines == 0:
        data["_atom_site"].drop(columns=["nmr_model"], inplace=True)

    return data


def _split_atom_line(
    line: str, nmr_model: int | None = None, allow_chimera: bool = True
) -> tuple:
    """Internal function to parse a single line belonging to 'ATOM', 'HETATM', and 'TER' lines

    Args:
        line (str): A 'ATOM', 'HETATM', and 'TER' line
        nmr_model (int|None; defaults to None): the NMR model number for the line
        allow_chimera (bool; defaults to True): try to parse as a Chimera-formatted PDB file.

    Returns:
        tuple: parsed values
    """
    atom_name = line[12:16].strip()
    alt_loc = line[16].strip()
    chain_id = line[21]
    residue_number = int(line[22:26].lstrip())
    insertion = line[26].strip()

    x_coord_s = line[30:38].lstrip()
    x_coord = float(x_coord_s) if x_coord_s else None
    y_coord_s = line[38:46].lstrip()
    y_coord = float(y_coord_s) if y_coord_s else None
    z_coord_s = line[46:54].lstrip()
    z_coord = float(z_coord_s) if z_coord_s else None

    occupancy_s = line[54:60].lstrip()
    occupancy = float(occupancy_s) if occupancy_s else None
    b_factor_s = line[60:66].lstrip()
    b_factor = float(b_factor_s) if b_factor_s else None

    segment_id = line[72:76].rstrip()
    element_symbol = line[76:78].lstrip()
    charge = line[78:80].strip()

    if allow_chimera:
        record_name = line[0:5].rstrip()
        if record_name == "HETAT":
            record_name = "HETATM"
        atom_number = int(line[5:11].lstrip("M").lstrip(" "))
        residue_name = line[17:21].strip()
    else:
        record_name = line[0:6].rstrip()
        atom_number = int(line[6:11].lstrip())
        residue_name = line[17:20].strip()

    return (
        record_name,
        atom_number,
        atom_name,
        alt_loc,
        residue_name,
        chain_id,
        residue_number,
        insertion,
        x_coord,
        y_coord,
        z_coord,
        occupancy,
        b_factor,
        segment_id,
        element_symbol,
        charge,
        nmr_model,
    )
