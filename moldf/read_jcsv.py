# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""JCSV format reading.

Reads a JCSV file into a dict of ``Pandas DataFrame`` s.
It is not limited to any molecular format.

"""
from __future__ import annotations

import csv
import os
import warnings
from collections import defaultdict

import pandas as pd  # type: ignore


def read_jcsv(
    jcsv_file: str | os.PathLike,
    category_names: list | None = None,
) -> dict[str, pd.DataFrame]:
    """Reads a JCSV file by name.

    Currently no molecular file repository has JCSV files so we can only read from a
    file name/path.

    Args:
        jcsv_file (required): JCSV file name/path.
        category_names (optional): a list of category names. If ``None``, all categories
            are read. Defaults to **None**.

    Returns:
        a dict of Pandas DataFrames for each category.

    Raises:
        TypeError: if ``category_names`` is not a list of strings.
        ValueError: if any of the ``category_names`` has double quotes or
            if the number of items in any line does not match the number of
            column names in the same category.
    """
    read_all = False
    if category_names is not None:
        if not isinstance(category_names, list):
            raise TypeError(f"{category_names} is not a list")
        for cat in category_names:
            if not isinstance(cat, str):
                raise TypeError(f"{cat} is not a str")
            elif '"' in cat:
                raise ValueError(f"{cat} has double quotes.")
    else:
        read_all = True

    results: dict[str, pd.DataFrame] = {}
    meta_found: bool | int = False
    with open(jcsv_file, "r") as jf:
        jf_reader = csv.reader(jf, delimiter=",", quotechar='"')
        for i, row in enumerate(jf_reader):
            if i == 0 and row[0] == "#jcsv_meta":
                meta_found = i + 1
                n_lines = _count_n_lines(jcsv_file)
            elif not meta_found:
                results = _read_jcsv_by_line(jcsv_file, category_names=category_names)
                break
            elif meta_found and i == meta_found:
                meta_col_names = row
                col_data: dict[str, list] = defaultdict(list)
            elif row[0][0] == "#":
                break
            elif meta_found:
                if len(row) != len(meta_col_names):
                    message = "Meta data has unmatched number"
                    message += f" of items in row '{row}' with the column"
                    message += f" names: {meta_col_names}"
                    raise ValueError(message)
                value: str | int = ""
                for col_name, value in zip(meta_col_names, row):
                    if col_name == "start_line_index":
                        value = int(value)
                    col_data[col_name].append(value)
    if meta_found:
        meta = list(zip(col_data["category"], col_data["start_line_index"]))
        start_line_index: int | str = 0
        for i, (category_name, start_line_index) in enumerate(meta[1:]):
            if read_all or (
                isinstance(category_names, list) and category_name in category_names
            ):
                start_line_index = int(start_line_index)
                skip_rows = [j for j in range(start_line_index)]
                if i < len(meta) - 2:
                    next_start = int(meta[i + 2][1])
                    ending_rows = [j for j in range(n_lines) if j > (next_start - 2)]
                    skip_rows.extend(ending_rows)
                results[category_name] = pd.read_csv(
                    jcsv_file, sep=",", quotechar='"', skiprows=skip_rows
                )

    return results


def _read_jcsv_by_line(
    jcsv_file: str | os.PathLike,
    category_names: list | None = None,
) -> dict[str, pd.DataFrame]:
    """Reads JCSV file line by line when the file has no meta data to select blocks.

    Args:
        jcsv_file (required): JCSV file name/path.
        category_names (optional): a list of category names. If ``None``, all categories
            are read. Defaults to **None**. It is passed by the ``read_jcsv`` caller, so
            it is not sanitized here.

    Returns:
        a dict of Pandas DataFrames for each category.

    Raises:
        ValueError: if the number of items in any line does not match the number of
            column names in the same category.
    """
    results: dict[str, pd.DataFrame] = {}

    read_all = False
    if category_names is None:
        read_all = True
    with open(jcsv_file, "r", encoding="utf-8") as jf:
        jf_reader = csv.reader(jf, delimiter=",", quotechar='"')
        req_cat_name_found: int | bool = False
        category_name = ""
        cols_data: dict[str, list] = defaultdict(list)
        if isinstance(category_names, list):
            category_names = list(category_names)
        for i, row in enumerate(jf_reader):
            if row[0][0] == "#":
                if read_all or (
                    isinstance(category_names, list) and row[0][1:] in category_names
                ):
                    if category_name:
                        results[category_name] = pd.DataFrame(cols_data)

                    category_name = row[0][1:]
                    cols_data = defaultdict(list)
                    req_cat_name_found = i + 1
                else:
                    req_cat_name_found = False
            elif req_cat_name_found and i == req_cat_name_found:
                col_names = row
            elif req_cat_name_found:
                if len(row) != len(col_names):
                    message = f"Category {category_name} has unmatched number"
                    message += f" of items in row '{row}' with the column"
                    message += f" names: {col_names}"
                    raise ValueError(message)
                for col_name, value in zip(col_names, row):
                    cols_data[col_name].append(value)

        if category_name:
            results[category_name] = pd.DataFrame(cols_data)

    if category_names is not None:
        for category_name in category_names:
            if category_name not in results:
                warnings.warn(
                    "Category {category_name} not in {jcsv_file}, not read",
                    RuntimeWarning,
                    stacklevel=2,
                )

    return results


def _count_n_lines(file_name: str | os.PathLike):
    """Gets the number of lines in a file.
    From https://stackoverflow.com/a/68385697/10094189

    Args:
        file_name (required): file name or path.
    """

    def _make_gen(reader):
        while True:
            b = reader(2**16)
            if not b:
                break
            yield b

    with open(file_name, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))

    return count
