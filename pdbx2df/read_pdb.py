from __future__ import annotations

import io
import os
import warnings

import pandas as pd  # type: ignore

IMPLEMENTED_PDB_CATS = ["_atom_site"]
ATOM_SITE = ["ATOM  ", "HETATM", "TER   "]
NMR_MDL = ["NUMMDL", "MODEL ", "ENDMDL"]


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

    atom_site_lines = ""
    n_nmr_modes = 0
    n_model_lines = 0
    n_endmdl_lines = 0
    nmr_model = None
    atom_site_line_nmr_model = []
    with open(pdb_file, "r", encoding="utf-8") as pf:
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
            if "_atom_site" in category_names and line[0:6] in NMR_MDL:
                if line[0:6] == "MODEL ":
                    n_model_lines += 1
                    nmr_model = int(line[6:].strip())
                elif line[0:6] == "ENDMDL":
                    n_endmdl_lines += 1
                else:
                    n_nmr_modes = int(line[6:].strip())
            if "_atom_site" in category_names and line[0:6] in ATOM_SITE:
                atom_site_lines += line
                if nmr_model is not None:
                    atom_site_line_nmr_model.append(nmr_model)
        if n_model_lines != n_endmdl_lines:
            raise ValueError(
                f"NMR records MODEL ({n_model_lines} lines) and ENDMDL ({n_endmdl_lines} lines) not matched."
            )
        if n_nmr_modes != n_model_lines:
            warnings.warn(
                f"The NUMMDL says {n_nmr_modes} NMR models, but only {n_model_lines} found.",
                RuntimeWarning,
                stacklevel=2,
            )
    if "_atom_site" in category_names:
        atom_site_buffer = io.StringIO()
        atom_site_buffer.writelines(atom_site_lines)
        atom_site_buffer.seek(0)
        if allow_chimera:
            col_widths = [5, 6, 1, 4, 1, 4, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2]
            col_names = [
                "record_name",
                "atom_number",
                "blank_1",
                "atom_name",
                "alt_loc",
                "residue_name",
                # "blank_2", # removed to be compatible with Chimera-formatted PDBs
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
        else:
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
        assert sum(col_widths) == 80
        assert len(col_widths) == len(col_names)
        df_atom_site = pd.read_fwf(atom_site_buffer, widths=col_widths, names=col_names)
        atom_site_buffer.close()

        col_names = [col_name for col_name in col_names if "blank" not in col_name]
        df_atom_site = df_atom_site[col_names]
        if allow_chimera:
            df_atom_site = _fix_chimera(df_atom_site)

        if atom_site_line_nmr_model:
            df_atom_site["nmr_model"] = atom_site_line_nmr_model

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


def _fix_chimera(df_atom_site: pd.DataFrame) -> pd.DataFrame:
    """Internal function to fix the 'record_name' and 'atom_number' columns in processing Chimera formatted PDBs.

    Args:
        df_atom_site (pd.DataFrame): original PDB dataframe.

    Returns:
        pd.DataFrame: updated PDB dataframe
    """  # noqa
    original_record_names = list(df_atom_site.record_name)
    record_names = [
        record_name if record_name in ["ATOM", "TER"] else record_name + "M"
        for record_name in original_record_names
    ]
    df_atom_site["record_name"] = record_names

    if df_atom_site["atom_number"].dtypes != "int64":
        original_atom_numbers = list(df_atom_site.atom_number)
        atom_numbers = [
            int(atom_number.lstrip("M").strip())
            for atom_number in original_atom_numbers
        ]
        df_atom_site["atom_number"] = atom_numbers

    return df_atom_site
