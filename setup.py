import os

__package_name__ = "moldf"
__author__ = "Ruibin Liu"

from setuptools import find_packages, setup  # type: ignore

if os.path.exists("README.md"):
    long_description = open("README.md").read()
else:
    long_description = "MolDF - Super lightweight and fast mmCIF/PDB/MOL2 file parser into Pandas DataFrames and backwards writer."  # noqa

with open("requirements.txt") as f:
    REQUIREMENTS = f.read().strip().split("\n")

ver = {}  # type: ignore
with open("moldf/version.py", "r") as vf:
    exec(vf.read(), ver)

setup(
    name=__package_name__,
    version=ver["__version__"],
    author=__author__,
    author_email="ruibinliuphd@gmail.com",
    description="Super lightweight and fast mmCIF/PDB/MOL2 file parser into Pandas DataFrames and backwards writer.",  # noqa
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Ruibin-Liu/MolDF",
    project_urls={
        "Bug Tracker": "https://github.com/Ruibin-Liu/MolDF/issues",
        "Documentation": "https://moldf.readthedocs.io/en/stable/",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    package_data={"moldf": ["covalent_bonds/*"]},
    install_requires=REQUIREMENTS,
    python_requires=">=3.7",
)
