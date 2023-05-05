import os
import sys

import pandas as pd  # type: ignore
from pdbx2df.read_pdbx import read_pdbx
from pdbx2df.write_pdbx import write_pdbx

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_write_pdbx():
    """
    Test write_pdbx function
    """
    file_path = [CFD, "test_files", "1VII.cif"]
    pdbx = read_pdbx(file_path, category_names=['all'])
    
    write_to = [CFD, "test_files", "Test.cif"]
    write_pdbx(pdbx, file_name=write_to)

