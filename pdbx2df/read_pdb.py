from __future__ import annotations

import io
import os
import warnings

import pandas as pd  # type: ignore

IMPLEMENTED_PDB_CATS = ["_atom_site"]
ATOM_SITE = ["ATOM  ", "HETATM", "TER   "]


def read_pdb(pdb_file: str | os.PathLike, category_names: list | None = None) -> dict:
    """
    Read a pdb file categories into Pandas DataFrame.

    Args:
        pdb_file (str|os.PathLike): file name for a PDB file.
        category_names (list|None; defaults None): a list of names for the categories as to the mmCIF file format.
            If None, "_atom_site" is used.
            To be consistent with the PDBx file format, the following category names are used to refer
            to block(s) in a PDB file and only they are supported:
            1. _atom_site: 'ATOM' and 'HETATM'
            2. TBD

    Returns:
        A dict of {category_name: pd.DataFrame of the info belongs to the category}
    """
    data: dict[str, pd.DataFrame] = {}
    if not category_names:
        category_names = ["_atom_site"]
    for category_name in category_names:
        if category_name not in IMPLEMENTED_PDB_CATS:
            implemented = ", ".join(IMPLEMENTED_PDB_CATS)
            raise ValueError(f"Only {implemented} are implented for the PDB format.")
        data[category_name] = pd.DataFrame()

    atom_site_lines = ""
    with open(pdb_file, "r") as pf:
        for line in pf:
            if "_atom_site" in category_names and line[0:3] == "TER":
                # This is for pdbfixer fixed PDB files whose TER lines are non standard.
                # This condition section can be removed if the above problem is fixed in pdbfixer.
                line = line.rstrip()
                line_len = len(line)
                line = line + " " * (80 - line_len) + "\n"
            elif len(line) != 81:
                warnings.warn(
                    f"Line {line} has non-standard length {len(line) - 1}, not 80; skipped",
                    RuntimeWarning,
                    stacklevel=2,
                )
                continue
            if "_atom_site" in category_names and line[0:6] in ATOM_SITE:
                atom_site_lines += line
    if "_atom_site" in category_names:
        atom_site_buffer = io.StringIO()
        atom_site_buffer.writelines(atom_site_lines)
        atom_site_buffer.seek(0)
        col_widths = [6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2]
        col_names = [
            "record_name",
            "atom_number",
            "blank_1",
            "atom_name",
            "alt_loc",
            "residue_name",
            "blank_2",
            "chain_id",
            "residue_number",
            "insertion",
            "blank_3",
            "x_coord",
            "y_coord",
            "z_coord",
            "occupancy",
            "b_factor",
            "blank_4",
            "segment_id",
            "element_symbol",
            "charge",
        ]
        df_atom_site = pd.read_fwf(atom_site_buffer, widths=col_widths, names=col_names)
        atom_site_buffer.close()
        col_names = [col_name for col_name in col_names if "blank" not in col_name]
        df_atom_site = df_atom_site[col_names]

        str_names = [
            "atom_name",
            "alt_loc",
            "residue_name",
            "chain_id",
            "insertion",
            "segment_id",
            "element_symbol",
        ]
        df_atom_site[str_names] = df_atom_site[str_names].fillna("")

        data["_atom_site"] = df_atom_site

    return data
