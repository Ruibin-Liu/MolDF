# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Write any dict of Pandas DataFrame to JCSV."""
from __future__ import annotations

import os

import pandas as pd  # type: ignore


def write_jcsv(
    data: dict[str, pd.DataFrame],
    file_name: str | os.PathLike | None = None,
    **kwargs,
) -> None:
    """Write a dict of ``Pandas DataFrame`` s into a JCSV file.
    See https://github.com/Ruibin-Liu/JCSV for definitions.

    Args:
        data (required): a ``dict`` of ``Pandas DataFrame`` s to write.
        file_name (optional): file name to write a JCSV file. If ``None``,
            ``moldf_output.jcsv`` will be used as the file name.
            Defaults to **None**.
        **kwargs: keyword arguments for ``pd.DataFrame.to_csv``. Invalid ones are ignored.
            Check https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html

    Raises:
        TypeError: if ``data`` is not a valid dict of ``DataFrame``.
        NotImplementedError: if ``"`` in any column name.
    """
    if not file_name:
        file_name = "moldf_output.jcsv"

    if not isinstance(data, dict):
        raise TypeError(f"'data' has to be a dict but {type(data)} is provided.")
    for key in data:
        if not isinstance(data[key], pd.DataFrame):
            raise TypeError(f"{data[key]} is not a Pandas DataFrame object.")
        if '"' in key:
            raise NotImplementedError(f'Column name {key} has " symbol, not supported.')

    exclusion = ["self", "path_or_buf", "df", "formatter"]
    to_csv_keys = [
        key for key in pd.DataFrame.to_csv.__code__.co_varnames if key not in exclusion
    ]
    kwargs = {k: v for k, v in kwargs.items() if k in to_csv_keys}

    with open(file_name, "w", encoding="utf-8") as out_file:
        for key, df in data.items():
            out_file.write(f"#{key}\n")
            out_file.write(df.to_csv(**kwargs))
