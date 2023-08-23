import os

__package_name__ = "pdbx2df"
__author__ = "Ruibin Liu"

from setuptools import find_packages, setup  # type: ignore

if os.path.exists("README.md"):
    long_description = open("README.md").read()
else:
    long_description = (
        "pdbx2df - A python package to parse PDBx file into Pandas DataFrames."
    )

with open("requirements.txt") as f:
    REQUIREMENTS = f.read().strip().split("\n")

ver = {}  # type: ignore
with open("pdbx2df/version.py", "r") as vf:
    exec(vf.read(), ver)

setup(
    name=__package_name__,
    version=ver["__version__"],
    author=__author__,
    author_email="ruibinliuphd@gmail.com",
    description="A python package to parse PDBx file into Pandas DataFrames.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Ruibin-Liu/pdbx2df",
    project_urls={"Bug Tracker": "https://github.com/Ruibin-Liu/pdbx2df/issues"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    install_requires=REQUIREMENTS,
    python_requires=">=3.7",
)
