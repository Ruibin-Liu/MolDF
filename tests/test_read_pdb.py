import os
import sys

import pandas as pd  # type: ignore

from pdbx2df.read_pdb import read_pdb

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_read_pdb():
    """
    Test read_pdb function
    """
    file_path = [CFD, "test_files", "5K9I.pdb"]
    test_file = f"{os.sep}".join(file_path)

    # Basic one category read
    pdb = read_pdb(test_file, category_names=["_atom_site"])

    compare_to = [CFD, "test_files", "5K9I.csv"]
    compare_to = f"{os.sep}".join(compare_to)
    df_expected = pd.read_csv(compare_to)
    str_names = [
        "atom_name",
        "alt_loc",
        "residue_name",
        "chain_id",
        "insertion",
        "segment_id",
        "element_symbol",
    ]
    df_expected[str_names] = df_expected[str_names].fillna("")
    pd.testing.assert_frame_equal(pdb["_atom_site"], df_expected)
