from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd  # type: ignore


def write_pdbx(pdbx: dict[str, pd.DataFrame], file_name: str) -> None:
    """
    Write a dict of Pandas DataFrames into a PDBx file.

    Args:
        pdbx (dict[str, pd.DataFrame]): a dict of Pandas DataFrames to write.
        file_name (str): file name to write the PDBx file.

    Returns:
        None
    """
    # Validate
    if not (type(pdbx) is dict):
        raise TypeError(f"pdbx has to be a dict but {type(pdbx)} is providied.")
    multi_record: dict[str, int] = defaultdict(bool)
    max_tag_length: dict[str, int] = defaultdict(int)
    for category_name, records in pdbx.items():
        if not (type(records) is pd.DataFrame):
            raise TypeError(
                f"pdbx values have to be Pandas Dataframes but {category_name} is a {type(records)}."
            )
        if len(records) > 1:
            multi_record[category_name] = True
        for col in records.columns:
            tag_length = len(category_name) + 1 + len(col)
            max_tag_length[category_name] = max(
                max_tag_length[category_name], tag_length
            )
    with open(file_name, "w") as f:
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
                    max_col_length[col] = max(records[col].str.len())
                    if records[col].str.contains(" ").any():
                        max_col_length[col] = max_col_length[col] + 2
                    f.write(f"{category_name}.{col}\n")
                for _, record in records.iterrows():
                    for col in records.columns:
                        content = record[col]
                        pad_length = max_col_length[col]
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
                    f.write("\n")

        f.write("#")
