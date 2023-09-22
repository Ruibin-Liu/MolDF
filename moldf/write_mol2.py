# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""MOL2 format writing.

Write a dict of ``Pandas DataFrame`` back to a MOL2 file.

Currently, only the ``MOLECULE``, ``ATOM``, and ``BOND`` categories can be written back.

"""
from __future__ import annotations

import os
from datetime import date

import pandas as pd  # type: ignore

from .version import __version__ as moldf_version

IMPLEMENTED_MOL2_CATS = ["MOLECULE", "ATOM", "BOND", "HEADER"]
"""MOL2 categories that are currently implemented."""


def write_mol2(
    mol2: dict[str, pd.DataFrame],
    file_name: str | os.PathLike | None = None,
) -> None:
    """Write a dict of ``Pandas DataFrame`` s into a MOL2 file.
    See https://is.muni.cz/th/fzk5s/dp_jakub_Vana.pdf p19 for column definitions.

    Args:
        mol2 (required): a ``dict`` of ``Pandas DataFrame`` s to write.
        file_name (optional): file name to write a MOL2 file. If ``None``,
            ``moldf_output.mol2`` will be used as the file name.
            Defaults to **None**.

    Raises:
        TypeError: if ``mol2`` is not a valid dict of ``DataFrame``.
        ValueError: if the ``mol2`` contains other than supported categories.
    """
    if not file_name:
        file_name = "moldf_output.mol2"

    if not isinstance(mol2, dict):
        raise TypeError(f"'mol2' has to be a dict but {type(mol2)} is provided.")

    implemented = ", ".join(IMPLEMENTED_MOL2_CATS)
    for key in mol2:
        if key not in IMPLEMENTED_MOL2_CATS:
            raise ValueError(f"Only {implemented} are implemented for the MOL2 format.")

        if not isinstance(mol2[key], pd.DataFrame):
            raise TypeError(f"{mol2[key]} is not a Pandas DataFrame object.")

    with open(file_name, "w", encoding="utf-8") as out_file:
        out_file.write("###\n")
        today = date.today().strftime("%Y-%m-%d")
        out_file.write(f"### Created by moldf v{moldf_version} {today}\n")

        if "HEADER" in mol2:
            df_header = mol2["HEADER"]
            for col_name in df_header.columns:
                header_line = df_header[col_name].to_list()[0]
                if col_name.startswith("info_"):
                    out_file.write(f"### Original header: {header_line}\n")
                else:
                    out_file.write(f"### {col_name}: {header_line}\n")

        out_file.write("###\n\n")
        out_file.write("@<TRIPOS>MOLECULE\n")
        if "MOLECULE" in mol2:
            df_molecule = mol2["MOLECULE"]
            mol_name = df_molecule.mol_name.to_list()[0]
            num_atoms = df_molecule.num_atoms.to_list()[0]
            num_bonds = df_molecule.num_bonds.to_list()[0]
            num_subst = df_molecule.num_subst.to_list()[0]
            num_feat = df_molecule.num_feat.to_list()[0]
            num_sets = df_molecule.num_sets.to_list()[0]
            mol_type = df_molecule.mol_type.to_list()[0]
            charge_type = df_molecule.charge_type.to_list()[0]

            status_bits, mol_comment = "", ""
            if status_bits in df_molecule.columns:
                status_bits = df_molecule.status_bits.to_list()[0]
                if mol_comment in df_molecule.columns:
                    mol_comment = df_molecule.mol_comment.to_list()[0]
            out_file.write(f"{mol_name}\n")
            out_file.write(
                f" {num_atoms} {num_bonds} {num_subst} {num_feat} {num_sets}\n"
            )
            out_file.write(f"{mol_type}\n")
            out_file.write(f"{charge_type}\n")
            if status_bits:
                out_file.write(f"{status_bits}\n")
                if mol_comment:
                    out_file.write(f"{mol_comment}\n")

        if "ATOM" in mol2:
            out_file.write("\n@<TRIPOS>ATOM\n")
            for _, row in mol2["ATOM"].iterrows():
                atom_id = row["atom_id"]
                atom_name = row["atom_name"]
                x = row["x"]
                y = row["y"]
                z = row["z"]
                atom_type = row["atom_type"]
                try:
                    subst_id = row["subst_id"]
                except KeyError:
                    subst_id = ""
                try:
                    subst_name = row["subst_name"]
                except KeyError:
                    subst_name = ""
                try:
                    charge = row["charge"]
                except KeyError:
                    charge = ""
                try:
                    status_bit = row["status_bit"]
                except KeyError:
                    status_bit = ""

                line = f"{atom_id:>6d} {atom_name:<8s} {x:>10.4f} {y:>10.4f} {z:>10.4f}"
                line += f" {atom_type:<9s} {subst_id:<2d} {subst_name:<7s} "
                line += f"{charge:>10.4f} {status_bit}"
                line = line.rstrip()
                out_file.write(line + "\n")

        if "BOND" in mol2:
            out_file.write("@<TRIPOS>BOND\n")
            for _, row in mol2["BOND"].iterrows():
                bond_id = row["bond_id"]
                origin_atom_id = row["origin_atom_id"]
                target_atom_id = row["target_atom_id"]
                bond_type = row["bond_type"]
                try:
                    status_bit = row["status_bit"]
                except KeyError:
                    status_bit = ""
                line = f"{bond_id:>6d} {origin_atom_id:>6d} {target_atom_id:>6d} "
                line += f"    {bond_type:1s} {status_bit:s}"
                line = line.rstrip()
                out_file.write(line + "\n")
