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
    expected = {"_atom_site": pd.read_pickle("test_files/5K9I.pkl")}
    pd.testing.assert_frame_equal(pdb["_atom_site"], expected["_atom_site"])
