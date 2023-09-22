# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Split a line in mmCIF files."""
from __future__ import annotations


def split_line(line: str, delimiter: str | None = None) -> list:
    """Splits a string line into tokens separated by ``delimiter`` s, assuming all
    ``'`` and ``"`` in the start character or following a `delimiter` are paired to
    quote a token.

    Args:
        line (required): line as a string
        delimiter (optional): ``delimiter`` to split the line.
            If ``None``, ``' '`` (one space) is used.
            Defaults to **None**.

    Returns:
        A list of tokens
    """
    if not delimiter:
        delimiter = " "
    words = []
    # without quotes, using shlex
    if '"' not in line and "'" not in line:
        if delimiter == " ":
            return line.split()
        return line.split(delimiter)

    # with quotes
    single_start = False
    double_start = False
    tmp: list[str] = []
    for i, char in enumerate(line):
        # quoted by single quotes ''
        if i == 0 and char == "'":  # line starting with '
            single_start = True
        elif (
            char == "'"
            and line[i - 1] == delimiter
            and not double_start
            and not single_start
        ):  # a new part quoted with '
            single_start = True
        elif char == "'" and single_start:  # a part quoted with ' ended
            single_start = False  # reset
        # quoted by double quotes ""
        elif i == 0 and char == '"':  # line starting with "
            double_start = True
        elif (
            char == '"'
            and line[i - 1] == delimiter
            and not double_start
            and not single_start
        ):  # a new part quoted with "
            double_start = True
        elif char == '"' and double_start:  # a part quoted with " ended
            double_start = False  # reset
        elif (
            char in [delimiter, "\n", "\r"] and not single_start and not double_start
        ):  # a part not quoted ended
            if tmp:
                words.append("".join(tmp))
                tmp = []
        else:  # Other characters including space in quotes
            tmp.append(char)
        if (
            tmp and i == len(line) - 1
        ):  # in case no '\r', '\n', or delimiter is at the end
            words.append("".join(tmp))
    if single_start or double_start:
        raise ValueError("Bad line: quotes not paired!")

    return words
