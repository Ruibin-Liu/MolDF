from __future__ import annotations

import os
from collections import defaultdict

import pandas as pd  # type: ignore

from .split_line import split_line


def read_pdbx(pdbx_file: str | os.PathLike, category_names: list | None = None) -> dict:
    """
    Read a pdbx file categories into Pandas DataFrame.

    Args:
        pdbx_file (str|Pathlike): file name for a PDBx file.
        category_names (list|None; defaults to None): a list of names for the categories in the mmCIF file format.
            If None, "all" is used and all categories will be processed.

    Returns:
        A dict of {category_name: pd.DataFrame of the info belongs to the category}
    """  # noqa
    data: dict[str, pd.DataFrame] = {}
    if not category_names:
        category_names = ["all"]
    category_name = ""
    with open(pdbx_file, "r") as pf:
        line = pf.readline()  # title line
        _ = line.strip().split("_")[1]
        processing_category = False
        while line:
            if line[0] == "#":
                processing_category = True
                loop_category = False
                category_lines = ""
                category_cols = []
                category_dict = defaultdict(list)
                line = pf.readline()  # category first line
                if not line:
                    break
                if line == "loop_\n":
                    loop_category = True
                    line = pf.readline()
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
                            line = (
                                '"' + line[1:]
                            )  # these ; symbols are not contents but a pair of container symbols for a single record
                            # but other ';' are parts of a single record, so we use " to make it easier to parse
                            if len(line.strip()) == 1:
                                in_multiple_line = False  # a ';' quoted record ends
                        elif in_multiple_line:
                            line = line.replace(
                                '"', "'"
                            )  # to avoid conflict with the above action
                        category_lines += line
                    line = pf.readline()
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
                            f"{pdbx_file} category {category_name} has irregular contents: {records}"
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
                            {len(records)} for {len(category_cols)}"""
                        )
                    record_data = records[len(category_cols) :]
                    for i, rd in enumerate(record_data):
                        col = i % len(category_cols)
                        if category_cols[col] in [
                            "pdbx_seq_one_letter_code",
                            "pdbx_seq_one_letter_code_can",
                        ]:
                            rd = rd.replace(
                                " ", ""
                            )  # delete the gaps in the one_letter_code
                        category_dict[category_cols[col]].append(rd)
                data[category_name] = pd.DataFrame(category_dict)
                for bn in category_names:
                    if bn not in data.keys():
                        break  # at least one required_category is not processed, continue the main while loop
                else:
                    break  # all required_categorys are processed, break the main while loop

            if not processing_category:
                line = pf.readline()

    return data
