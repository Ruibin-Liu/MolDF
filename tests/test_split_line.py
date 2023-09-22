# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
"""Tests for splitting lines in PDBx/mmCIF files."""
import sys

import pytest  # type: ignore

from moldf.split_line import split_line

sys.path.append("..")


def test_split_line():
    """
    Test split_line function
    """
    # Test line start with '  and " quoted by ''
    line = """'a "T"' in Test
    """
    expected = ['a "T"', "in", "Test"]
    assert split_line(line) == expected

    # Test line start with " and ' quoted by ""
    line = """"a 't'" in Test
    """
    expected = ["a 't'", "in", "Test"]
    assert split_line(line) == expected

    # Test '' and "" both used as quotes
    line = """
    DA  'DNA linking'       y "2'-DEOXYADENOSINE-5'-MONOPHOSPHATE" ? 'C10 H14 N5 O6 P' 331.222
    """  # noqa
    expected = [
        "DA",
        "DNA linking",
        "y",
        "2'-DEOXYADENOSINE-5'-MONOPHOSPHATE",
        "?",
        "C10 H14 N5 O6 P",
        "331.222",
    ]
    assert split_line(line) == expected

    # Test ending with '\n'
    line = "Test\n return symbol in middle"
    expected = ["Test", "return", "symbol", "in", "middle"]
    assert split_line(line) == expected

    # Test ending with '\t'
    line = "This is a test\t"
    expected = ["This", "is", "a", "test"]
    assert split_line(line) == expected

    # Test ending with multiple delimiter
    line = "This is a test    "
    expected = ["This", "is", "a", "test"]
    assert split_line(line) == expected

    # Test no ending symbol
    line = "This is a test"
    expected = ["This", "is", "a", "test"]
    assert split_line(line) == expected

    # Test not not paired ' or "
    line = 'This "test fails'
    with pytest.raises(ValueError) as e:
        assert split_line(line)
    assert str(e.value) == "Bad line: quotes not paired!"
    line = "This 'test fails"
    with pytest.raises(ValueError) as e:
        assert split_line(line)
    assert str(e.value) == "Bad line: quotes not paired!"

    # Test , as delimiter
    line = 'This,is,  a,"test"'
    expected = ["This", "is", "  a", "test"]
    assert split_line(line, delimiter=",") == expected
