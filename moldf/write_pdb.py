# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""PDB format writing.

Write a dict of ``Pandas DataFrame`` back to a PDB file.

Currently, only the ``_atom_site`` category can be written back.

"""
from __future__ import annotations

import os
import warnings
from datetime import date

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

from .version import __version__ as moldf_version

IMPLEMENTED_PDB_CATS = ["_atom_site"]
"""PDB categories that are currently implemented."""


def write_pdb(
    pdb: dict[str, pd.DataFrame],
    file_name: str | os.PathLike | None = None,
    allow_chimera: bool = False,
) -> None:
    """Write a dict of ``Pandas DataFrame`` s into a PDB file.

    Args:
        pdb (required): a ``dict`` of ``Pandas DataFrame`` s to write.
        file_name (optional): file name to write a PDB file. If ``None``,
            ``moldf_output.pdb`` will be used as the file name.
            Defaults to **None**.
        allow_chimera (optional): whether to allow writing to Chimera-formatted PDB
            files. Defaults to **False**.

    Raises:
        TypeError: if ``pdb`` is not a valid dict of ``DataFrame``.
        ValueError: if the ``pdb`` contains other than supported categories.
    """
    if not file_name:
        file_name = "moldf_output.pdb"

    if not isinstance(pdb, dict):
        raise TypeError(f"pdb has to be a dict but {type(pdb)} is provided.")

    implemented = ", ".join(IMPLEMENTED_PDB_CATS)
    for key in pdb.keys():
        if key not in IMPLEMENTED_PDB_CATS:
            raise ValueError(f"Only {implemented} are implemented for the PDB format.")

    with open(file_name, "w", encoding="utf-8") as out_file:
        today = date.today().strftime("%Y-%m-%d")
        padding = " "
        tag = f"CREATED WITH moldf v{moldf_version} {today}  "
        header = f"REMARK   1{tag:>70}\n"
        out_file.write(header)

        if "_atom_site" in pdb.keys():
            df_atom_site = pdb["_atom_site"]
            str_names = [
                "atom_name",
                "alt_loc",
                "residue_name",
                "chain_id",
                "insertion",
                "segment_id",
                "element_symbol",
                "charge",
            ]
            df_atom_site[str_names] = df_atom_site[str_names].fillna("")
            n_nmr_models = 0
            if "nmr_model" in df_atom_site.columns:
                n_nmr_models = len(df_atom_site.nmr_model.unique())
                out_file.write(f"NUMMDL{n_nmr_models:>6}{padding:>68}\n")
            prev_nmr_model = None
            n_nmr_written = 0
            for _, row in df_atom_site.iterrows():
                record_name = row[
                    "record_name"
                ]  # 1-6 (normal) or 1-5 (Chimera);char;left
                atom_number = row[
                    "atom_number"
                ]  # 7-11 (normal) or 6-11 (Chimera);char;left
                # 12 (blank)
                atom_name = row[
                    "atom_name"
                ]  # 13-16 (2/1-l elements left-13/14 unless == 4);char
                alt_loc = row["alt_loc"]  # 17;char
                residue_name = row[
                    "residue_name"
                ]  # 18-20 (normal) or 18-21 (Chimera);char;right
                # 21 (blank) or NA (Chimera)
                chain_id = row["chain_id"]  # 22;char
                residue_number = row["residue_number"]  # 23-26;int;right
                insertion = row["insertion"]  # 27;char
                # 28-30 (blank)
                x_coord = row["x_coord"]  # 31-38 (8.3);right
                y_coord = float(row["y_coord"])  # 39-46 (8.3);right
                z_coord = float(row["z_coord"])  # 47-54 (8.3);right
                occupancy = float(row["occupancy"])  # 55-60 (6.2);right
                b_factor = float(row["b_factor"])  # 61-66 (6.2);right
                # 67-72 (blank)
                segment_id = row["segment_id"]  # 73-76;char;left
                element_symbol = row["element_symbol"]  # 77-78;char;right
                charge = row["charge"]  # 79-80;char(?);right?
                if n_nmr_models > 0:
                    nmr_model = row["nmr_model"]
                    if prev_nmr_model is None:
                        out_file.write(f"MODEL {nmr_model:>6}{padding:>68}\n")
                    elif nmr_model != prev_nmr_model:
                        out_file.write(f"ENDMDL{padding:>74}\n")
                        n_nmr_written += 1
                        out_file.write(f"MODEL {nmr_model:>6}{padding:>68}\n")
                    prev_nmr_model = nmr_model
                if len(element_symbol) == 1 and len(atom_name) < 4:
                    atom_name = " " + atom_name
                if allow_chimera:
                    if len(record_name) == 6:
                        warnings.warn(
                            f"Record name {record_name} length was 6 and is truncated to {record_name[:-1]}",  # noqa
                            RuntimeWarning,
                            stacklevel=2,
                        )
                        record_name = record_name[:-1]
                    atom_site_line = f"{record_name:<5s}{atom_number:>6d} {atom_name:<4s}{alt_loc:<1s}"  # noqa
                    if len(residue_name) < 4:
                        residue_name = residue_name + " "
                    atom_site_line += f"{residue_name:>4s}{chain_id:<1s}{residue_number:>4d}{insertion:<1s}   "  # noqa
                else:
                    atom_site_line = f"{record_name:<6s}{atom_number:>5d} {atom_name:<4s}{alt_loc:<1s}"  # noqa
                    atom_site_line += f"{residue_name:>3s} {chain_id:<1s}{residue_number:>4d}{insertion:<1s}   "  # noqa
                if np.isnan(x_coord):
                    x_coord_e, y_coord_e, z_coord_e, occupancy_e, b_factor_e = (
                        " ",
                        " ",
                        " ",
                        " ",
                        " ",
                    )
                    atom_site_line += f"{x_coord_e:>8s}{y_coord_e:>8s}{z_coord_e:>8s}{occupancy_e:>6s}"  # noqa
                    atom_site_line += f"{b_factor_e:>6s}      {segment_id:<4s}{element_symbol:>2s}{charge:>2s}\n"  # noqa
                else:
                    atom_site_line += (
                        f"{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}{occupancy:6.2f}"
                    )
                    atom_site_line += f"{b_factor:6.2f}      {segment_id:<4s}{element_symbol:>2s}{charge:>2s}\n"  # noqa
                out_file.write(atom_site_line)
            if n_nmr_models > 0:
                out_file.write(f"ENDMDL{padding:>74}\n")
