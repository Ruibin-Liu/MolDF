# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""PDBx/mmCIF format reading.

Reads a PDBx file into a dict of ``Pandas DataFrame`` s.
The dict keys are read from the mmCIF category names automatically.
For each category, the column names of the corresponding ``DataFrame``
are read form the category attributes automatically.

For example:
::

  _audit_conform.dict_name       mmcif_pdbx.dic
  _audit_conform.dict_version    5.355
  _audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic

In this category, the category name is ``_audit_conform``, so the returned dict
key is ``_audit_conform``. The category attributes are ``dict_name``,
``dict_version``, and ``dict_location``, so the returned dict value as a
``DataFrame`` has the exact column names.
"""
from __future__ import annotations

import io
import os
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd  # type: ignore

from .split_line import split_line

AF2_MODEL = 4
"""For AlphaFold structures, the version to use."""


def read_pdbx(
    pdbx_file: str | os.PathLike[str] | None = None,
    pdb_id: str | None = None,
    save_pdbx_file: bool = True,
    pdbx_file_dir: str | os.PathLike | None = None,
    category_names: list | None = None,
) -> dict[str, pd.DataFrame]:
    """Reads a ``.cif`` file's categories into a ``dict`` of ``Pandas DataFrame`` s.

    Args:
        pdb_id (optional): PDB/Uniprot ID. Required if ``pdbx_file`` is ``None``.
            Defaults to **None**.
        pdbx_file (optional): file name for a PDBx/mmCIF file. Used over `pdb_id`.
            Defaults to **None**.
        category_names (optional): a list of categories in the mmCIF file format.
            If ``None``, ``all`` is used and all categories will be processed.
            Defaults to **None**.
        save_pdbx_file (optional): whether to save the fetched PDBx file from RCSB
            to ``pdbx_file_dir``. Defaults to **False**.
        pdbx_file_dir (optional): directory to save fetched PDBx files. If ``None`` but
            ``save_pdbx_file`` is ``True``, './PDBx_files' is used.
            Defaults to **None**.

    Returns:
        A dict of ``Pandas DataFrame`` s corresponding to required categories.

    Raises:
        ValueError: if none of ``pdb_id`` or ``pdbx_file`` is provided, or if ``pdb_id``
            is given but cannot the PDB file cannot be downloaded from RCSB, or the
            PDB file is corrupted like no end-line symbol, or some content is irregular.
        FileNotFoundError: if ``pdbx_file`` cannot be found.
    """
    data: dict[str, pd.DataFrame] = {}
    if pdb_id is None and pdbx_file is None:
        raise ValueError("At least one of pdb_id and pdbx_file has to be given.")
    elif pdbx_file is None:
        pdb_id = str(pdb_id).upper()
        if len(pdb_id) == 4:
            pdbx_file_url = f"https://files.rcsb.org/view/{pdb_id.upper()}.cif"
        elif len(pdb_id) < 4:
            pdbx_file_url = f"https://files.rcsb.org/ligands/view/{pdb_id.upper()}.cif"
        else:
            pdbx_file_url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id.upper()}-F1-model_v{AF2_MODEL}.cif"
        try:
            with urllib.request.urlopen(pdbx_file_url) as response:
                raw_data = response.read()
            text = raw_data.decode("utf-8")
            pdbx_file_handle: io.TextIOWrapper | io.StringIO = io.StringIO(text)
            if save_pdbx_file:
                if pdbx_file_dir is None:
                    pdbx_file_dir = "./PDBx_files"
                pdbx_file_dir = Path(pdbx_file_dir)
                if not pdbx_file_dir.exists():
                    pdbx_file_dir.mkdir(parents=True, exist_ok=True)
                file_path = Path(pdbx_file_dir, f"{pdb_id}.cif")
                with open(file_path, "w", encoding="utf-8") as p_file:
                    p_file.write(text)
        except urllib.error.HTTPError as http_error:
            raise ValueError(
                f"Cannot download PDBx file from url {pdbx_file_url}."
            ) from http_error
    else:
        pdbx_file = Path(pdbx_file)
        if not pdbx_file.exists():
            raise FileNotFoundError(f"File {pdbx_file} not found.")
        pdbx_file_handle = open(pdbx_file, "r", encoding="utf-8")
    if not category_names:
        category_names = ["all"]
    category_name = ""
    with pdbx_file_handle:
        line = pdbx_file_handle.readline()  # title line
        _ = line.strip().split("_")[1]
        processing_category = False
        while line:
            if line[0] == "#":
                processing_category = True
                loop_category = False
                category_lines = ""
                category_cols = []
                category_dict = defaultdict(list)
                line = pdbx_file_handle.readline()  # category first line
                if not line:
                    break
                if line == "loop_\n":
                    loop_category = True
                    line = pdbx_file_handle.readline()
                in_multiple_line = (
                    False  # lines quoted by ';', they need to be specially treated
                )
                while line.strip() not in ["#", "#   #", "##"]:
                    if line[0] == "_" and not in_multiple_line:
                        category_name, category_nth_col = (
                            line.strip().split(" ")[0].split(".")
                        )
                        category_cols.append(category_nth_col)
                    if (
                        category_name in category_names or "all" in category_names
                    ):  # in a required category
                        if line[0] == ";":
                            in_multiple_line = True  # a ';' quoted record begins
                            line = line.replace(
                                '"', "'"
                            )  # to avoid conflict with the next action
                            line = '"' + line[1:]
                            # these ; symbols are not contents but a pair of container
                            # symbols for a single record
                            # but other ';' are parts of a single record, so we use "
                            # to make it easier to parse
                            if len(line.strip()) == 1:
                                in_multiple_line = False  # a ';' quoted record ends
                        elif in_multiple_line:
                            line = line.replace(
                                '"', "'"
                            )  # to avoid conflict with the above action
                        category_lines += line
                    line = pdbx_file_handle.readline()
                    if not line:
                        raise ValueError(
                            f"{pdbx_file} not normally ended. Download it again?"
                        )
                if "all" not in category_names and category_name not in category_names:
                    continue
                category_lines = category_lines.replace(
                    "\n", " "
                )  # cannot rely on newline to parse anyway.
                records = split_line(category_lines)
                if not loop_category:
                    if len(records) != 2 * len(category_cols):
                        raise ValueError(
                            f"{pdbx_file} category {category_name} has irregular contents: {records}"  # noqa
                        )
                    for i, rec in enumerate(records):
                        if i % 2 == 0:
                            col_name = rec.split(".")[1]
                            if col_name in [
                                "pdbx_seq_one_letter_code",
                                "pdbx_seq_one_letter_code_can",
                            ]:
                                records[i + 1] = records[i + 1].replace(" ", "")
                            category_dict[rec.split(".")[1]] = [records[i + 1]]
                else:
                    if (
                        records[len(category_cols)][0] == "_"
                        or records[len(category_cols) - 1][0] != "_"
                        or not len(records) % len(category_cols) == 0
                    ):
                        raise ValueError(
                            f"""{pdbx_file} category {category_name} has irregular contents:
                            {len(records)} for {len(category_cols)}"""  # noqa
                        )
                    record_data = records[len(category_cols) :]
                    for i, i_data in enumerate(record_data):
                        col = i % len(category_cols)
                        if category_cols[col] in [
                            "pdbx_seq_one_letter_code",
                            "pdbx_seq_one_letter_code_can",
                        ]:
                            i_data = i_data.replace(
                                " ", ""
                            )  # delete the gaps in the one_letter_code
                        category_dict[category_cols[col]].append(i_data)
                data[category_name] = pd.DataFrame(category_dict)
                for c_name in category_names:
                    if c_name not in data:
                        break
                        # at least one required_category is not processed,
                        # continue the main while loop
                else:
                    break
                    # all required_categories are processed,
                    # break the main while loop

            if not processing_category:
                line = pdbx_file_handle.readline()
    if "_atom_site" in data:
        col_dtypes = {
            "id": "int",
            "label_seq_id": "int",
            "Cartn_x": "float",
            "Cartn_y": "float",
            "Cartn_z": "float",
            "occupancy": "float",
            "B_iso_or_equiv": "float",
            "auth_seq_id": "int",
            "pdbx_PDB_model_num": "int",
            "pdbx_sifts_xref_db_num": "int",
        }
        col_dtypes = {
            col_name: data_type
            for col_name, data_type in col_dtypes.items()
            if col_name in data["_atom_site"]
        }

        data["_atom_site"] = data["_atom_site"].astype(col_dtypes)
    return data
