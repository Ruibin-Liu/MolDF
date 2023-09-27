# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""PDB format reading.

Reads a PDB file, including Chimera compatible ones, into a dict of
``Pandas DataFrame`` s.

Atom coordinates and sequences are currently supported by the following category names:

    ``_atom_site``:  ``ATOM``, ``HETATM``, ``TER``, ``NUMMDL``, ``MODEL``, and
    ``ENDMDL`` lines.

    ``_seq_res``: ``SEQRES`` lines.

"""
from __future__ import annotations

import io
import os
import urllib.request
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

from .constants import AMINO_ACIDS

IMPLEMENTED_PDB_CATS = ["_atom_site", "_seq_res"]
"""PDB categories that are currently implemented."""

ATOM_SITE = ("ATOM", "HETATM", "TER")
"""``_atom_site`` primary lines."""

NMR_MDL = ("NUMMDL", "MODEL", "ENDMDL")
"""``_atom_site`` additional lines."""

AF2_MODEL = 4
"""For AlphaFold structures, the version to use."""


def read_pdb(
    pdb_file: str | os.PathLike[str] | None = None,
    pdb_id: str | None = None,
    category_names: list | None = None,
    save_pdb_file: bool = True,
    pdb_file_dir: str | os.PathLike | None = None,
    allow_chimera: bool = True,
    need_ter_lines: bool = True,
) -> dict[str, pd.DataFrame]:
    """Reads a ``.pdb`` file's categories into a ``dict`` of ``Pandas DataFrame`` s.

    Args:
        pdb_id (optional): PDB/Uniprot ID. Required if ``pdb_file`` is ``None``.
            Defaults to **None**.
        pdb_file (optional): file name for a PDB file. Used over `pdb_id`.
            Defaults to **None**.
        category_names (optional): a list of categories similar to the mmCIF format.
            If ``None``, ``_atom_site`` is used.
            To be consistent with the PDBx format, the following category names are
            used to refer to block(s) in a PDB file and only they are supported:

                1. ``_atom_site``: ``ATOM``, ``HETATM``, and ``TER`` lines and possible
                ``NUMMDL``, ``MODEL``, and ``ENDMDL`` lines.

                2. ``_seq_res``: ``SEQRES`` lines.

            Defaults to **None**.
        save_pdb_file (optional): whether to save the fetched PDB file from RCSB
            to ``pdb_file_dir``. Defaults to **True**.
        pdb_file_dir (optional): directory to save fetched PDB files. If ``None`` but
            ``save_pdb_file`` is ``True``, './PDB_files' is used.
            Defaults to **None**.
        allow_chimera (optional): whether to allow Chimera-formatted PDB files.
            Defaults to **True**.
        need_ter_lines (optional): whether to read the ``TER`` lines into the
            ``DataFrame``. Defaults to **True**.

    Returns:
        A dict of ``Pandas DataFrame`` s corresponding to required categories.

    Raises:
        ValueError: if none of ``pdb_id`` or ``pdbx_file`` is provided, or if ``pdb_id``
            is given but cannot the PDB file cannot be downloaded from RCSB,
        NotImplementedError: if `category_names` not a subset of allowed names.
        FileNotFoundError: if ``pdbx_file`` cannot be found.
    """
    data: dict[str, pd.DataFrame] = {}
    if pdb_id is None and pdb_file is None:
        raise ValueError("At least one of pdb_id and pdb_file has to be given.")
    elif pdb_file is None:
        pdb_id = str(pdb_id).upper()
        if len(pdb_id) == 4:
            pdb_file_url = f"https://files.rcsb.org/view/{pdb_id.upper()}.pdb"
        else:
            pdb_file_url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id.upper()}-F1-model_v{AF2_MODEL}.pdb"
        try:
            with urllib.request.urlopen(pdb_file_url) as response:
                raw_data = response.read()
            text = raw_data.decode("utf-8")
            pdb_file_handle: io.TextIOWrapper | io.StringIO = io.StringIO(text)
            if save_pdb_file:
                if pdb_file_dir is None:
                    pdb_file_dir = "./PDB_files"
                pdb_file_dir = Path(pdb_file_dir)
                if not pdb_file_dir.exists():
                    pdb_file_dir.mkdir(parents=True, exist_ok=True)
                file_path = Path(pdb_file_dir, f"{pdb_id}.pdb")
                with open(file_path, "w", encoding="utf-8") as p_file:
                    p_file.write(text)
        except urllib.error.HTTPError as http_error:
            raise ValueError(
                f"Cannot download PDB file from url {pdb_file_url}."
            ) from http_error
    else:
        pdb_file = Path(pdb_file)
        if not pdb_file.exists():
            raise FileNotFoundError(f"File {pdb_file} not found.")
        pdb_file_handle = open(pdb_file, "r", encoding="utf-8")
    if category_names is None:
        category_names = ["_atom_site"]
    for category_name in category_names:
        if category_name not in IMPLEMENTED_PDB_CATS:
            implemented = ", ".join(IMPLEMENTED_PDB_CATS)
            raise NotImplementedError(
                f"Only {implemented} are implemented for the PDB format."
            )
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
        all_atom_site_reads = ("ATOM", "HETATM") + NMR_MDL
    else:
        all_atom_site_reads = ATOM_SITE + NMR_MDL  # type: ignore

    if pdb_file is not None:  # This check is just for mypy
        file_stat = os.stat(pdb_file)
        total_lines = int(file_stat.st_size / 81)
    else:
        total_lines = 100000
    n_record = 0
    n_lines_till_atom_lines = 0
    with pdb_file_handle:
        line = pdb_file_handle.readline()
        for line in pdb_file_handle:
            if line.startswith(all_atom_site_reads) and "_atom_site" in category_names:
                if (not n_record) and line.startswith(("ATOM", "HETATM", "MODEL")):
                    array = np.zeros(
                        total_lines - n_lines_till_atom_lines, column_name_types
                    )
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
            elif line.startswith("SEQRES") and "_seq_res" in category_names:
                chain_residues: dict[str, list[str]] = defaultdict(list)
                chain_lengths: dict[str, int] = defaultdict(int)
                while line.startswith("SEQRES"):
                    items: list[str] = line.strip().split()
                    chain_id: str = items[2]
                    chain_length: int = int(items[3])
                    chain_lengths[chain_id] = chain_length
                    chain_residues[chain_id].extend(items[4:])
                    if not n_record and "_atom_site" in category_names:
                        n_lines_till_atom_lines += 1
                    line = pdb_file_handle.readline()

            if not n_record and "_atom_site" in category_names:
                n_lines_till_atom_lines += 1

    if (
        "_atom_site" in category_names
        and num_nmr_models != n_model_lines
        and num_nmr_models != 0
    ):
        warnings.warn(
            f"The NUMMDL says {num_nmr_models} NMR models, but {n_model_lines} found.",
            RuntimeWarning,
            stacklevel=2,
        )

    if "_atom_site" in category_names:
        df_atom_site = pd.DataFrame.from_records(array[:n_record])
        data["_atom_site"] = df_atom_site
        if not n_model_lines:
            data["_atom_site"].drop(columns=["nmr_model"], inplace=True)

    if "_seq_res" in category_names:
        chain_sequences: dict[str, str] = {}
        for chain_id, residues in chain_residues.items():
            chain_sequences[chain_id] = "".join([AMINO_ACIDS[i] for i in residues])
            if len(residues) != chain_lengths[chain_id]:
                warnings.warn(
                    f"Cols 14-17 value {chain_lengths[chain_id]} contrasts real length {len(residues)}.",  # noqa
                    RuntimeWarning,
                    stacklevel=2,
                )
        chain_ids = chain_sequences.keys()
        chain_num_residues = [len(chain_sequences[chain_id]) for chain_id in chain_ids]
        chain_seqs = [chain_sequences[chain_id] for chain_id in chain_ids]
        data["_seq_res"] = pd.DataFrame(
            {
                "chain_id": chain_ids,
                "chain_sequence": chain_seqs,
                "chain_length": chain_num_residues,
            }
        )

    return data


def _split_atom_line(
    line: str,
    nmr_model: int = -1,
    allow_chimera: bool = True,
    is_ter_line: bool = False,
) -> tuple:
    """Internal function to parse a single line belonging to ``ATOM``, ``HETATM``,
    or ``TER`` lines.

    Args:
        line (required): A ``ATOM``, ``HETATM``, or ``TER`` line.
        nmr_model (optional): the NMR model number for the line; ``-1`` means
            not an NMR model. Defaults to **-1**.
        allow_chimera (optional): try to parse as a Chimera-formatted PDB file.
            Defaults to **True**
        is_ter_line (optional): whether the line starts with ``TER``.
            Defaults to **False**.

    Returns:
        parsed values as a tuple.
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
