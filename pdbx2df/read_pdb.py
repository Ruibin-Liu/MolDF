from __future__ import annotations

import os
import warnings

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

IMPLEMENTED_PDB_CATS = ["_atom_site"]
ATOM_SITE = ("ATOM", "HETATM", "TER")
NMR_MDL = ("NUMMDL", "MODEL", "ENDMDL")


def read_pdb(
    pdb_file: str | os.PathLike,
    category_names: list | None = None,
    allow_chimera: bool = True,
    need_ter_lines: bool = True,
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
        need_ter_lines (bool; defaults to True): whether to read the TER lines into the DataFrame.

    Returns:
        A dict of {category_name: pd.DataFrame of the info belongs to the category}
    """  # noqa
    data: dict[str, pd.DataFrame] = {}
    if category_names is None:
        category_names = ["_atom_site"]
    for category_name in category_names:
        if category_name not in IMPLEMENTED_PDB_CATS:
            implemented = ", ".join(IMPLEMENTED_PDB_CATS)
            raise ValueError(f"Only {implemented} are implemented for the PDB format.")
        data[category_name] = pd.DataFrame()

    num_nmr_models = 0
    n_model_lines = 0
    nmr_model = -1  # -1 means there is no NMR model
    n_endmdl_lines = 0
    # atom_site_rows = []
    column_name_types = [
        ("record_name", "U6"),
        ("atom_number", "i4"),
        ("atom_name", "U4"),
        ("alt_loc", "U1"),
        ("residue_name", "U4"),
        ("chain_id", "U1"),
        ("residue_number", "i2"),
        ("insertion", "U1"),
        ("x_coord", "f4"),
        ("y_coord", "f4"),
        ("z_coord", "f4"),
        ("occupancy", "f4"),
        ("b_factor", "f4"),
        ("segment_id", "U4"),
        ("element_symbol", "U2"),
        ("charge", "U2"),
        ("nmr_model", "i2"),
    ]

    if not need_ter_lines:
        all_reads = ("ATOM", "HETATM") + NMR_MDL
    else:
        all_reads = ATOM_SITE + NMR_MDL  # type: ignore

    file_stat = os.stat(pdb_file)
    total_lines = int(file_stat.st_size / 81)
    # array = np.zeros(total_lines, column_name_types)
    n_record = 0
    n_lines = 0
    with open(pdb_file, "r", encoding="utf-8") as pf:
        for line in pf:
            if line.startswith(all_reads):
                if (not n_record) and line.startswith(("ATOM", "MODEL")):
                    array = np.zeros(total_lines - n_lines, column_name_types)
                if line.startswith(("ATOM", "HETATM")):
                    array[n_record] = _split_atom_line(
                        line,
                        nmr_model=nmr_model,
                        allow_chimera=allow_chimera,
                        is_ter_line=False,
                    )
                    n_record += 1
                elif need_ter_lines and line.startswith("TER"):
                    array[n_record] = _split_atom_line(
                        line,
                        nmr_model=nmr_model,
                        allow_chimera=allow_chimera,
                        is_ter_line=True,
                    )
                    n_record += 1
                elif line.startswith("MODEL"):
                    n_model_lines += 1
                    nmr_model = int(line[6:].strip())
                elif line.startswith("ENDMDL"):
                    n_endmdl_lines += 1
                    if num_nmr_models and n_endmdl_lines == num_nmr_models:
                        break
                else:  # line.startswith("NUMMDL")
                    num_nmr_models = int(line[6:].strip())
            if not n_record:
                n_lines += 1

    if num_nmr_models != n_model_lines:
        warnings.warn(
            f"The NUMMDL says {num_nmr_models} NMR models, but only {n_model_lines} found.",
            RuntimeWarning,
            stacklevel=2,
        )

    df = pd.DataFrame.from_records(array[:n_record])
    data["_atom_site"] = df

    if not n_model_lines:
        data["_atom_site"].drop(columns=["nmr_model"], inplace=True)

    return data


def _split_atom_line(
    line: str,
    nmr_model: int = -1,
    allow_chimera: bool = True,
    is_ter_line: bool = False,
) -> tuple:
    """Internal function to parse a single line belonging to 'ATOM', 'HETATM', or 'TER' lines

    Args:
        line (str): A 'ATOM', 'HETATM', or 'TER' line.
        nmr_model (int; defaults to -1): the NMR model number for the line; default -1 means not an NMR model
        allow_chimera (bool; defaults to True): try to parse as a Chimera-formatted PDB file.
        is_ter_line (bool; defaults to False): whether the line starts with 'TER'.

    Returns:
        tuple: parsed values
    """
    if not (is_ter_line or allow_chimera):
        return (
            line[0:6],
            line[6:11],
            line[12:16],
            line[16],
            line[17:20],
            line[21],
            line[22:26],
            line[26],
            line[30:38],
            line[38:46],
            line[46:54],
            line[54:60],
            line[60:66],
            line[72:76],
            line[76:78],
            line[78:80],
            nmr_model,
        )
    elif is_ter_line:
        if allow_chimera:
            record_name = line[0:5] + " "
            residue_name = line[17:21]
        else:
            record_name = line[0:6]
            residue_name = line[17:20]
        return (
            record_name,
            line[6:11],
            line[12:16],
            line[16],
            residue_name,
            line[21],
            line[22:26],
            line[26],
            None,
            None,
            None,
            None,
            None,
            line[72:76],
            line[76:78],
            line[78:80],
            nmr_model,
        )
    else:
        if line[0:5] == "HETAT":
            record_name = "HETATM"
            atom_number = line[6:11]
        else:
            record_name = line[0:5] + " "
            atom_number = line[5:11]
        return (
            record_name,
            atom_number,
            line[12:16],
            line[16],
            line[17:21],
            line[21],
            line[22:26],
            line[26],
            line[30:38],
            line[38:46],
            line[46:54],
            line[54:60],
            line[60:66],
            line[72:76],
            line[76:78],
            line[78:80],
            nmr_model,
        )
