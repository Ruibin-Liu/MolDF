from collections import defaultdict

import pandas as pd

from .split_line import split_line


def read_pdbx(pdbx_file: str, block_names: list = ["_pdbx_poly_seq_scheme"]) -> dict:
    """
    Read a pdbx file block into Pandas DataFrame.

    Args:
        pdbx_file: file path for a PDBx file.
        block_names: a list of names for the blocks (between two '#' lines) in a PDBx file that need to be read.

    Returns:
        A dict of {block_name: pd.DataFrame of the info belongs to the block}
    """
    data = {}
    with open(pdbx_file, "r") as pf:
        line = pf.readline()  # title line
        _ = line.strip().split("_")[1]
        processing_block = False
        while line:
            if line[0] == "#":
                processing_block = True
                loop_block = False
                block_lines = ""
                block_cols = []
                block_dict = defaultdict(list)
                line = pf.readline()  # block first line
                if not line:
                    break
                if line == "loop_\n":
                    loop_block = True
                    line = pf.readline()
                in_multiple_line = (
                    False  # lines quoted by ';', they need to be specially treated
                )
                while line[0] != "#":
                    if line[0] == "_":
                        block_name, block_nth_col = (
                            line.strip().split(" ")[0].split(".")
                        )
                        block_cols.append(block_nth_col)
                    if (
                        block_name in block_names or "all" in block_names
                    ):  # in a required block
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
                        block_lines += line
                    line = pf.readline()
                    if not line:
                        raise ValueError(
                            f"{pdbx_file} not normally ended. Download it again?"
                        )
                if block_name not in block_names:
                    continue
                block_lines = block_lines.replace(
                    "\n", " "
                )  # cannot rely on newline to parse anyway.
                records = split_line(block_lines)
                if not loop_block:
                    if len(records) != 2 * len(block_cols):
                        raise ValueError(
                            f"{pdbx_file} block {block_name} has irregular contents: {records}"
                        )
                    for i, rec in enumerate(records):
                        if i % 2 == 0:
                            col_name = rec.split(".")[1]
                            if col_name in [
                                "pdbx_seq_one_letter_code",
                                "pdbx_seq_one_letter_code_can",
                            ]:
                                records[i + 1] = records[i + 1].replace(" ", "")
                            block_dict[rec.split(".")[1]] = [records[i + 1]]
                else:
                    if (
                        records[len(block_cols)][0] == "_"
                        or records[len(block_cols) - 1][0] != "_"
                        or not len(records) % len(block_cols) == 0
                    ):
                        raise ValueError(
                            f"""{pdbx_file} block {block_name} has irregular contents:
                            {len(records)} for {len(block_cols)}"""
                        )
                    record_data = records[len(block_cols) :]
                    for i, rd in enumerate(record_data):
                        col = i % len(block_cols)
                        if block_cols[col] in [
                            "pdbx_seq_one_letter_code",
                            "pdbx_seq_one_letter_code_can",
                        ]:
                            rd = rd.replace(
                                " ", ""
                            )  # delete the gaps in the one_letter_code
                        block_dict[block_cols[col]].append(rd)
                data[block_name] = pd.DataFrame(block_dict)
                for bn in block_names:
                    if bn not in data.keys():
                        break  # at least one required_block is not processed, continue the main while loop
                else:
                    break  # all required_blocks are processed, break the main while loop

            if not processing_block:
                line = pf.readline()

    return data
