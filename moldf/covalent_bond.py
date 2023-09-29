# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Gets covalent bonds by covalent radii cutoff or ligand templates."""
from __future__ import annotations

import os
import warnings
from itertools import product
from pathlib import Path

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

from moldf.read_pdbx import read_pdbx


def get_covalent_bond_cutoffs(
    element_symbols: list | set, single_radii_set: str | None = None
) -> tuple[dict, dict, dict]:
    """Gets the DataFrame for additive covalent radii of relevant elements.

    Args:
        element_symbols (required): a list or set of element symbols whose covalent
            radii cutoffs are returned.
        single_radii_set (optional): radii sets to use. If ``None``, ``single_C``
            is used as to Cordero (PMID 18478144). Another option is ``single_PA``
            which refers to Pyykkö's studies (PMID 19058281;19856342;15832398,
            and doi:10.1103/PhysRevB.85.024115). Defaults to **None**.

    Raises:
        ValueError: if the ``single_radii_set`` is not valid.

    Returns:
        a tuple of the three dictionaries of the distance cutoffs for the pairs of
            the queried element symbols: (single_bonds, double_bonds, triple_bonds).
    """
    if single_radii_set is None:
        single_radii_set = "single_C"
    if single_radii_set not in ["single_C", "single_PA"]:
        message = "The 'single_radii_set' has to be one of "
        message += "'single_C' and 'single_PA', but "
        message += f"{single_radii_set} is provided."
        raise ValueError(message)

    this_dir = os.path.dirname(__file__)
    cov_radii_file = Path(this_dir).absolute() / "covalent_bonds" / "covalent_radii.csv"
    covalent_radii = pd.read_csv(cov_radii_file)
    covalent_radii.replace("-", np.nan, inplace=True)
    covalent_radii[
        ["atomic_number", single_radii_set, "double", "triple"]
    ] = covalent_radii[["atomic_number", single_radii_set, "double", "triple"]].astype(
        float,
    )

    element_set = set(element_symbols)
    if "D" in element_set or "T" in element_set:
        element_set.update({"H"})
    covalent_radii = covalent_radii[covalent_radii.element_symbol.isin(element_set)]

    single_bonds = {}
    for first, second in product(element_set, repeat=2):
        first_radius = covalent_radii[covalent_radii.element_symbol == first][
            single_radii_set
        ].to_list()[0]
        second_radius = covalent_radii[covalent_radii.element_symbol == second][
            single_radii_set
        ].to_list()[0]
        first = first.rjust(2)
        second = second.rjust(2)
        single_bonds[(first, second)] = (
            (first_radius + second_radius) / 100 + 0.2
        ) ** 2
        single_bonds[(second, first)] = single_bonds[(first, second)]
        # TODO: make it possible to use 'single_PA' as radii criteria option.
        # Is it possible to fit a function purely based on closest 6 distances?
    double_bonds = {}
    for first, second in product(element_set, repeat=2):
        first_radius = covalent_radii[
            covalent_radii.element_symbol.str.upper() == first
        ]["double"].to_list()[0]
        second_radius = covalent_radii[
            covalent_radii.element_symbol.str.upper() == second
        ]["double"].to_list()[0]
        first = first.rjust(2)
        second = second.rjust(2)
        double_bonds[(first, second)] = (
            (first_radius + second_radius) / 100 + 0.05
        ) ** 2
        double_bonds[(second, first)] = double_bonds[(first, second)]
    triple_bonds = {}
    for first, second in product(element_set, repeat=2):
        first_radius = covalent_radii[
            covalent_radii.element_symbol.str.upper() == first
        ]["triple"].to_list()[0]
        second_radius = covalent_radii[
            covalent_radii.element_symbol.str.upper() == second
        ]["triple"].to_list()[0]
        first = first.rjust(2)
        second = second.rjust(2)
        triple_bonds[(first, second)] = (
            (first_radius + second_radius) / 100 + 0.05
        ) ** 2
        triple_bonds[(second, first)] = triple_bonds[(first, second)]

    return single_bonds, double_bonds, triple_bonds


def get_residue_template(
    residue_name: str,
    parent_name: str | None = None,
    residue_template_file: str | os.PathLike | None = None,
    save_template_file: bool = True,
    template_file_dir: str | os.PathLike | None = None,
) -> dict:
    """Gets the intra- covalent bonds of a residue or ligand.

    Args:
        residue_name (required): residue or ligand name to query.
        parent_name (optional): parent name of the residue or ligand. This is useful if
            the ``residue_name`` in the PDB file is used by RCSB for a different residue
            or ligand. If ``None``, it is the same as ``residue_name``.
            Defaults to **None**.
        residue_template_file (optional): residue or ligand template file name/path.
            If ``None``, this function queries RCSB by the ``parent_name``. If it's not
            a path, this function looks at the current working directory first and then
            ``template_file_dir`` if it is not ``None``.
            Defaults to **None**.
        save_template_file (optional): whether to save the downloaded (from RCSB)
            template file. Defaults to **True**.
        template_file_dir (optional): directory to save fetched template files.
            If ``None`` but ``save_template_file`` is ``True``, './template_files'
            directory is used. Defaults to **None**.

    Returns:
        a dictionary of bonds for the residue. The dictionary key format is
            (atom_id_1: str, atom_id_2: str), and the dictionary value format is
            (bond_order: str, is_aromatic: bool, stereo_flag: str)

    Raises:
        ValueError: if a ``.cif`` file for the residue/ligand cannot be downloaded from
            RCSB for any reason if ``residue_template_file`` is ``None``. Check the url
            provided in the error message to download manually if possible.
        FileNotFoundError: if ``residue_template_file`` cannot be found if not ``None``.

    Warnings:
        RuntimeWarning: if the template file contains other residue/ligands.
    """
    if parent_name is None:
        parent_name = residue_name
    parent_name = parent_name.upper()

    if residue_template_file is None:
        if template_file_dir is None:
            template_file_dir = "./template_files"
            template_file_dir = Path(template_file_dir)
            if not template_file_dir.exists():
                template_file_dir.mkdir(parents=True, exist_ok=True)
        residue_template = read_pdbx(
            pdb_id=parent_name,
            save_pdbx_file=save_template_file,
            pdbx_file_dir=template_file_dir,
            category_names=["_chem_comp_bond"],
        )["_chem_comp_bond"]
        residue_template_file = Path(template_file_dir) / f"{parent_name}.cif"
    else:
        residue_template_file = Path(residue_template_file)
        if not residue_template_file.exists():
            if template_file_dir is not None:
                template_file_dir = Path(template_file_dir)
                residue_template_file = Path(template_file_dir, residue_template_file)
            if not residue_template_file.exists():
                raise FileNotFoundError(f"File {residue_template_file} not found.")
        residue_template = read_pdbx(
            pdbx_file=residue_template_file,
            category_names=["_chem_comp_bond"],
        )["_chem_comp_bond"]

    if residue_template.comp_id.unique() != [parent_name]:
        message = f"The {residue_template_file} contains other compounds, "
        message += f"but only {parent_name} is used."
        warnings.warn(message, RuntimeWarning, stacklevel=2)
        residue_template = residue_template[residue_template.comp_id == parent_name]

    results: dict = {}
    atom_id_1 = residue_template.atom_id_1
    atom_id_2 = residue_template.atom_id_2
    bond_order = residue_template.value_order
    aromatic_flag = residue_template.pdbx_aromatic_flag
    stereo_flag = residue_template.pdbx_stereo_config
    for id1, id2, bo, af, sf in zip(
        atom_id_1, atom_id_2, bond_order, aromatic_flag, stereo_flag
    ):
        is_aromatic = True
        if af == "Y":
            is_aromatic = False
        results[(id1, id2)] = (
            bo.upper(),  # bond_order
            is_aromatic,
            sf,  # stereo_flag
        )
        results[(id2, id1)] = results[(id1, id2)]

    return results
