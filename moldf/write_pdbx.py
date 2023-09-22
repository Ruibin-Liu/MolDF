# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""PDBx/mmCIF format writing.

Write a dict of ``Pandas DataFrame`` back to a PDBx file.

"""
from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import pandas as pd  # type: ignore


def write_pdbx(
    pdbx: dict[str, pd.DataFrame], file_name: str | os.PathLike | None = None
) -> None:
    """Writes a ``dict`` of ``Pandas DataFrame`` s into a PDBx file.

    Args:
        pdbx (required): a ``dict`` of ``Pandas DataFrame`` s to write.
        file_name (optional): file name to write a PDBx file. If ``None``,
            ``moldf_output.cif`` will be used as the file name.
            Defaults to **None**.

    Raises:
        TypeError: if ``pdbx`` is not a valid dict of ``DataFrame``.
    """
    if not file_name:
        file_name = "moldf_output.cif"

    if not isinstance(pdbx, dict):
        raise TypeError(f"pdbx has to be a dict but {type(pdbx)} is provided.")
    multi_record: dict[str, int] = defaultdict(bool)
    max_tag_length: dict[str, int] = defaultdict(int)
    for category_name, records in pdbx.items():
        if not isinstance(records, pd.DataFrame):
            raise TypeError(
                f"pdbx values have to be Pandas DataFrames but {category_name} is a {type(records)}."  # noqa
            )
        if len(records) > 1:
            multi_record[category_name] = True
        for col in records.columns:
            tag_length = len(category_name) + 1 + len(col)
            max_tag_length[category_name] = max(
                max_tag_length[category_name], tag_length
            )
    with open(file_name, "w", encoding="utf-8") as f:
        # write header
        target_name = Path(file_name).name
        if ".cif" == target_name[-4:]:
            f.write(f"data_{target_name[:-4]}\n")
        else:
            f.write(f"data_{target_name}\n")
        # write each category
        for category_name, records in pdbx.items():
            # categories that only have a record
            f.write("#\n")
            if not multi_record[category_name]:
                for col in records.columns:
                    tag = f"{category_name}.{col}"
                    f.write(f"{tag:{max_tag_length[category_name]+3}}")
                    content = records[col].iloc[0]
                    if '"' in content and "'" in content:
                        raise ValueError(
                            f"'{content}' cannot be written into a PDBx file."
                        )
                    elif "'" in content:
                        content = f'"{content}"'
                    elif '"' in content:
                        content = f"'{content}'"
                    elif " " in content:
                        content = f"'{content}'"

                    content_length = len(content)
                    if tag_length + content_length > 130:
                        content = content.strip('"').strip("'")
                        f.write("\n;")
                        if category_name == "_struct_ref":
                            for i in range(0, content_length // 80):
                                f.write(f"{content[80*i:80*(i+1)]}\n")
                        else:
                            f.write(f"{content}\n")
                        f.write(";\n")
                    else:
                        f.write(f"{content}\n")

            # categories that have multiple records
            else:
                max_col_length = defaultdict(int)
                for col in records.columns:
                    if records[col].dtype == "int":
                        max_col_length[col] = len(str(max(records[col])))
                    elif col == "occupancy":
                        max_col_length[col] = 4
                    elif records[col].dtype == "float":
                        max_int_width = max(
                            len(str(int(max(records[col])))),
                            len(str(int(min(records[col])))),
                        )
                        max_col_length[col] = max_int_width + 4
                        if col == "B_iso_or_equiv":
                            max_col_length[col] = max_int_width + 3
                    else:
                        max_col_length[col] = max(records[col].str.len())
                        if records[col].str.contains(" ").any():
                            max_col_length[col] = max_col_length[col] + 1

                    f.write(f"{category_name}.{col}\n")
                for _, record in records.iterrows():
                    for col in records.columns:
                        content = record[col]
                        pad_length = max_col_length[col]
                        if isinstance(content, str):
                            if '"' in content and "'" in content:
                                raise ValueError(
                                    f"'{content}' cannot be written into a PDBx file."
                                )
                            elif "'" in content:
                                content = f'"{content}"'
                            elif '"' in content:
                                content = f"'{content}'"
                            elif " " in content:
                                content = f"'{content}'"

                            f.write(f"{content:<{pad_length+1}}")
                        elif isinstance(content, int):
                            f.write(f"{content:<{pad_length+1}}")
                        elif isinstance(content, float) and col in [
                            "Cartn_x",
                            "Cartn_y",
                            "Cartn_z",
                        ]:
                            f.write(f"{content:<{pad_length+1}.3f}")
                        elif isinstance(content, float) and col in [
                            "occupancy",
                            "B_iso_or_equiv",
                        ]:
                            f.write(f"{content:<{pad_length+1}.2f}")
                    f.write("\n")

        f.write("#")
