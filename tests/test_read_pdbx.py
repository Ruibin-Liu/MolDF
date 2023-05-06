import os
import sys

import pandas as pd  # type: ignore

from pdbx2df.read_pdbx import read_pdbx

sys.path.append("..")
CFD = os.path.dirname(__file__)


def test_read_pdbx():
    """
    Test read_pdbx function
    """
    file_path = [CFD, "test_files", "1VII.cif"]
    test_file = f"{os.sep}".join(file_path)

    # Basic one category read
    pdbx = read_pdbx(test_file, category_names=["_entry"])
    expected = {"_entry": pd.DataFrame({"id": ["1VII"]})}
    pd.testing.assert_frame_equal(pdbx["_entry"], expected["_entry"])

    # Multiple categories read
    pdbx = read_pdbx(test_file, category_names=["_entry", "_entity_name_com"])
    expected = {
        "_entry": pd.DataFrame({"id": ["1VII"]}),
        "_entity_name_com": pd.DataFrame(
            {"entity_id": ["1"], "name": ["HP-36, R42-76"]}
        ),
    }
    pd.testing.assert_frame_equal(pdbx["_entry"], expected["_entry"])
    pd.testing.assert_frame_equal(
        pdbx["_entity_name_com"], expected["_entity_name_com"]
    )

    # category with ';'
    pdbx = read_pdbx(test_file, category_names=["_pdbx_nmr_refine"])
    expected = {
        "_pdbx_nmr_refine": pd.DataFrame(
            {
                "entry_id": ["1VII"],
                "method": ["DISTANCE GEOMETRY SIMULATED ANNEALING"],
                "details": [
                    "THE X-PLOR (R6)1/6 NOE POTENTIAL WAS USED FOR NOE'S INVOLVING NON-STEREOSPECIFICALLY ASSIGNED METHYL, METHYLENE, AND AROMATIC PROTONS. NO ATTRACTIVE POTENTIALS WERE USED IN CALCULATING THE STRUCTURES. THE VAN DER WAALS CUTOFF USED FOR THE X-PLOR REPEL FUNCTION WAS 0.75 ANGSTROMS. AFTER DISTANCE GEOMETRY AND REGULARIZATION, EACH STRUCTURE WAS SUBJECTED TO ONE ROUND OF SIMULATED ANNEALING FROM 2000K TO 100K OVER 2000 STEPS. THIS IS THE AVERAGE OF 29 STRUCTURES MINIMIZED USING ONLY REPULSIVE POTENTIALS. "  # noqa
                ],
                "software_ordinal": ["1"],
            }
        )
    }
    pd.testing.assert_frame_equal(
        pdbx["_pdbx_nmr_refine"], expected["_pdbx_nmr_refine"]
    )
