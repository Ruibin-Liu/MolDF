"""Common IO utilities."""

from io import TextIOWrapper
from pathlib import Path 

def open_file(
    pdb_file: Path,
) -> TextIOWrapper:
    """Open a file for reading, handling compression if necessary.

    Args:
        pdb_file: Path to the file to open.

    Returns:
        File handle for the opened file.
    """
    
    if not pdb_file.exists():
        raise FileNotFoundError(f"File {pdb_file} not found.")
    
    # Check compression 
    if pdb_file.suffix == ".gz":
        import gzip
        pdb_file_handle = gzip.open(pdb_file, "rt", encoding="utf-8")
    elif pdb_file.suffix == ".bz2":
        import bz2
        pdb_file_handle = bz2.open(pdb_file, "rt", encoding="utf-8")
    elif pdb_file.suffix == ".xz":
        import lzma
        pdb_file_handle = lzma.open(pdb_file, "rt", encoding="utf-8")
    elif pdb_file.suffix == ".zip":
        import zipfile
        with zipfile.ZipFile(pdb_file, "r") as zip_file:
            pdb_file_handle = TextIOWrapper(zip_file.open(zip_file.namelist()[0], "r"))
    else:
        pdb_file_handle = open(pdb_file, "rt", encoding="utf-8")
    
    return pdb_file_handle
