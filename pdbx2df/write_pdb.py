from __future__ import annotations

import os
from datetime import date

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

from .version import __version__ as pdbx2df_version

IMPLEMENTED_PDB_CATS = ["_atom_site"]


def write_pdb(
    pdb: dict[str, pd.DataFrame], file_name: str | os.PathLike | None = None
) -> None:
    """
    Write a dict of Pandas DataFrames into a PDB file.

    Args:
        pdb (dict[str, pd.DataFrame]): a dict of Pandas DataFrames to write.
        file_name (str|os.PathLike|None; defaults to None): file name to write a PDB file.
            If None, "pdbx2df_output.pdb" will be used as the file name.

    Returns:
        None
    """  # noqa
    if not file_name:
        file_name = "pdbx2df_output.pdb"

    if not isinstance(pdb, dict):
        raise TypeError(f"pdb has to be a dict but {type(pdb)} is providied.")

    implemented = ", ".join(IMPLEMENTED_PDB_CATS)
    for key in pdb.keys():
        if key not in IMPLEMENTED_PDB_CATS:
            raise ValueError(f"Only {implemented} are implented for the PDB format.")
    if "_atom_site" in pdb.keys():
        df = pdb["_atom_site"]
        df.fillna("", inplace=True)
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
        df["blank_1"] = ""
        df["blank_2"] = ""
        df["blank_3"] = ""
        df["blank_4"] = ""
        df = df[col_names]
        formats = "%-6s%+5s%+1s%4s%+1s%+3s%+1s%+1s%+4s%+1s%-3s%+8s%+8s%+8s%+6s%+6s%+6s%-4s%+2s%+2s"
        today = date.today().strftime("%Y-%m-%d")
        header = "{:>6s} {:>3s} {:>1s} {:>1s} {:>9s} {:>49s}".format(
            "REMARK", "", "1", "", "", f"CREATED WITH PDBX2DF {pdbx2df_version} {today}"
        )
        np.savetxt(
            file_name,
            df.values,
            fmt=formats,
            header=header,
            comments="",
            footer="END" + " " * 77,
        )
