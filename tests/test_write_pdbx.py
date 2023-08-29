import filecmp
import os
import sys

from pdbx2df.read_pdbx import read_pdbx
from pdbx2df.write_pdbx import write_pdbx

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_write_pdbx():
    """
    Test write_pdbx function
    """
    test_file = [CFD, "test_files", "1VII.cif"]
    test_file = f"{os.sep}".join(test_file)
    pdbx = read_pdbx(test_file, category_names=["all"])

    compare_to = [CFD, "test_files", "1VII_pdbx2df.cif"]
    compare_to = f"{os.sep}".join(compare_to)

    write_to = [CFD, "test_files", "1VII_test.cif"]
    write_to = f"{os.sep}".join(write_to)

    write_pdbx(pdbx, file_name=write_to)
    assert filecmp.cmp(compare_to, write_to)
    os.remove(write_to)
