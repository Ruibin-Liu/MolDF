from __future__ import annotations


def split_line(line: str, delimeter: str | None = None) -> list:
    """
    Split a string line into tokens separated by delimeters, assuming all ' and " in the start character
    or following a delimeter are paired to quote a token.

    Args:
        line (str): line as a string
        delimeter (str|None; defaults to None): delimeter to split the line; if None, delimeter == ' '.

    Returns:
        A list of tokens: words
    """  # noqa
    if not delimeter:
        delimeter = " "
    words = []
    # wihtout quotes, using shlex
    if '"' not in line and "'" not in line:
        if delimeter == " ":
            return line.split()
        return line.split(delimeter)

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
            and line[i - 1] == delimeter
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
            and line[i - 1] == delimeter
            and not double_start
            and not single_start
        ):  # a new part quoted with "
            double_start = True
        elif char == '"' and double_start:  # a part quoted with " ended
            double_start = False  # reset
        elif (
            char in [delimeter, "\n", "\r"] and not single_start and not double_start
        ):  # a part not quoted ended
            if tmp:
                words.append("".join(tmp))
                tmp = []
        else:  # Other characters including space in quotes
            tmp.append(char)
        if (
            tmp and i == len(line) - 1
        ):  # in case no '\r', '\n', or delimeter is at the end
            words.append("".join(tmp))
    if single_start or double_start:
        raise ValueError("Bad line: quotes not paired!")

    return words
