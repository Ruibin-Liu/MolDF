# MolDF
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/MolDF
""" ``PDBDataFrame`` as a subclass of ``Pandas DataFrame``.

Several features are added to make PDB data more accessible and selectable:

1. Properties like ``sequences``, ``heavy_atoms``, ``backbone``, and ``water`` are
directly accessed by ``.`` operation.

2. Atom selection by using methods whose names are just the column names plus ``s``
(plural form). For example, selecting atoms by names is simply
``df.atom_names([names])`` where ``atom_name`` is the column name
and ``atom_names`` is the selection function. Each selection returns a
``PDBDataFrame`` object as well, which means we can chain selections one by one
like ``df.atom_names([names]).residue_numbers([numbers])``.

3. Distance matrix as a ``@property`` and ``@classmethod``.
"""
from __future__ import annotations

import functools
import warnings
from collections import defaultdict
from collections.abc import Iterable
from itertools import combinations

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from scipy.spatial.distance import cdist, pdist, squareform  # type: ignore
from scipy.spatial.transform import Rotation  # type: ignore
from typing_extensions import Self

from .constants import AMINO_ACIDS, ELEMENT_MASSES
from .covalent_bond import get_covalent_bond_cutoffs, get_residue_template

RESIDUE_CODES = AMINO_ACIDS
"""dict[str, str], turn 3-, 2-, and 1-letter residue codes to 1-letter codes."""


PDBX_COLS = {
    "record_name": "group_PDB",
    "atom_number": "id",
    "atom_name": "label_atom_id",
    "alt_loc": "label_alt_id",
    "residue_name": "label_comp_id",
    "chain_id": "label_asym_id",
    "residue_number": "label_seq_id",
    "insertion": "pdbx_PDB_ins_code",
    "x_coord": "Cartn_x",
    "y_coord": "Cartn_y",
    "z_coord": "Cartn_z",
    "occupancy": "occupancy",
    "b_factor": "B_iso_or_equiv",
    "segment_id": "label_entity_id",
    "element_symbol": "type_symbol",
    "charge": "pdbx_formal_charge",
    "nmr_model": "pdbx_PDB_model_num",
}
"""dict[str, str], PDB and mmCIF column name dictionary."""


class PDBDataFrame(pd.DataFrame):
    """Pandas DataFrame with extended attributes and methods for PDB data.

    It enables Pythonic atom selection methods and convenient ``.`` accessing to common
    PDB structure properties.

    Args:
        *args: all ``pd.DataFrame`` positional arguments. For example, the
            ``_atom_site`` dataframe returned by reading a PDB file.
        pdb_format (optional): PDB format in the underlying provided data.
            If ``None``, ``PDB`` is assumed. Defaults to **None**.
        use_squared_distance (optional): whether to use squared distance
            when calculating distance matrix. Defaults to **True**.
        use_square_form (optional): whether to use a square matrix
            for the distance matrix. Defaults to **False**.
        **kwargs: all ``pd.DataFrame`` acceptable keyword arguments.

    Returns:
        A ``PDBDataFrame`` instance.

    Example
    -------
    >>> from moldf import read_pdb, PDBDataFrame
    >>> pdb = read_pdb(pdb_id='1vii')
    >>> pdb_df = pdb['_atom_site']
    >>> pdb_df = PDBDataFrame(pdb_df)

    Warnings
    --------
    This subclass uses a custom ``__hash__`` function for caching some calculations.
    And thus a custom ``__eq__`` function is also implemented. For other typical
    ``DataFrame`` operations, use those ``.all()``, ``.any()``, ``.bool()`` functions
    to do comparison.

    """

    _metadata = [
        "_use_squared_distance",
        "_use_square_form",
        "_is_chimera",
        "_RESIDUE_CODES",
        "_ELEMENT_MASSES",
        "_pdb_format",
    ]

    def __init__(
        self,
        *args,
        pdb_format: str | None = None,
        use_squared_distance: bool = True,
        use_square_form: bool = False,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        self._pdb_format: str | None = pdb_format
        if self._pdb_format is None:
            self._pdb_format = "PDB"
        if self._pdb_format.lower() in ["mmcif", "pdbx"]:
            self._pdbx_to_pdb()
        self._use_squared_distance: bool = use_squared_distance
        self._use_square_form: bool = use_square_form
        self._hash_random_state: int = 0
        self._is_chimera = False
        self._RESIDUE_CODES: dict[str, str] = {}
        self._ELEMENT_MASSES: dict[str, float] = {}
        self._ter_line_removed: bool = False
        self._atoms: Self | None = None

    @property
    def _constructor(self):
        return PDBDataFrame

    def __hash__(self) -> int:
        """Uses head X coords to hash; for distance matrix calculation cache."""
        try:
            sample_atom_numbers = self.sample(
                5, random_state=self.hash_random_state, replace=True
            )["atom_number"]
        except ValueError:
            sample_atom_numbers = []
        return hash(tuple(self[self.atom_number.isin(sample_atom_numbers)].x_coord))

    def __eq__(self, other) -> bool:
        """Uses head X coords to compare; for distance matrix calculation cache."""
        return self.__hash__() == other.__hash__()

    def _pdbx_to_pdb(self, keep_original: bool = False):
        """Converts PDBx '_atom_site' DataFrame to PDB format.

        Args:
            keep_original (optional): whether to keep the original columns in the PDBx
                '_atom_site' DataFrame. Defaults to **False**.
        """
        pdbx_cols = {k: v for k, v in PDBX_COLS.items() if v in self.columns}
        for pdb_name, pdbx_name in pdbx_cols.items():
            self[pdb_name] = self[pdbx_name]
        if not keep_original:
            drop_columns = [col for col in self.columns if col not in pdbx_cols.keys()]
            self.drop(columns=drop_columns, inplace=True)
        self._pdb_format = "PDBx"

    @property
    def pdb_format(self) -> str:
        """
        The format of the current PDBDataFrame.
        """
        return self._pdb_format  # type: ignore

    @property
    def RESIDUE_CODES(self) -> dict[str, str]:
        """
        A dict of ``residue_name`` as keys and ``residue_code`` as values, where
        ``residue_code`` is a 1-character code used in sequences. **Settable**.
        """
        if not self._RESIDUE_CODES:
            res_name_width = 3
            if self.is_chimera:
                res_name_width = 4
            self._RESIDUE_CODES = {
                res.upper().ljust(res_name_width): code
                for res, code in RESIDUE_CODES.items()
            }

        return self._RESIDUE_CODES

    @RESIDUE_CODES.setter
    def RESIDUE_CODES(self, residue_codes: dict[str, str]) -> None:
        res_name_width = 3
        if self.is_chimera:
            res_name_width = 4

        self._RESIDUE_CODES = {
            res.upper().ljust(res_name_width): code
            for res, code in residue_codes.items()
        }

    @property
    def ELEMENT_MASSES(self) -> dict[str, float]:
        """
        A dict of ``element_symbol`` as keys and ``element_mass`` as values, where
        ``element_mass`` is taken from NIST. **Settable**.
        """
        if not self._ELEMENT_MASSES:
            self._ELEMENT_MASSES = {
                e.upper().rjust(2): mass for e, mass in ELEMENT_MASSES.items()
            }
        return self._ELEMENT_MASSES

    @ELEMENT_MASSES.setter
    def ELEMENT_MASSES(self, element_masses: dict[str, float]):
        self._ELEMENT_MASSES = {
            e.upper().rjust(2): mass for e, mass in element_masses.items()
        }

    @property
    def is_chimera(self) -> bool:
        """
        Whether the original read-in PDB was Chimera compatible format.
            The main effect is the ``residue_name`` str width is 4 in the
            Chimera compatible format instead of 3 as in the standard PDB
            format.  **Not settable**.
        """
        residue_name_set = self.atoms.residue_name.unique()
        try:
            residue_name_set = self.atoms.residue_name.unique()
            if len(residue_name_set) > 0 and len(residue_name_set[0]) == 4:
                self._is_chimera = True
        except AttributeError:
            pass  # 'residue_name' not in self.coords

        return self._is_chimera

    @property
    def hash_random_state(self) -> int:
        """The ``random_state`` used in the ``__hash__`` function. **Settable**."""
        return self._hash_random_state

    @hash_random_state.setter
    def hash_random_state(self, random_state: int) -> None:
        self._hash_random_state = random_state

    @property
    def use_squared_distance(self) -> bool:
        """
        Whether R or R^2 is used in distance matrix calculations.
            Using R^2 saves computation time. **Settable**.
        """
        return self._use_squared_distance

    @use_squared_distance.setter
    def use_squared_distance(self, use_r2: bool) -> None:
        self._use_squared_distance = use_r2

    @property
    def use_square_form(self) -> bool:
        """
        Whether the distance matrix will be in a square form.
            Using square form consumes less memory. **Settable**.
        """
        return self._use_square_form

    @use_square_form.setter
    def use_square_form(self, square_form: bool) -> None:
        self._use_square_form = square_form

    @property
    def atoms(self) -> Self:
        """Gets atoms in the ``ATOM`` and ``HETATM`` entries.
        In other words, removing 'TER' lines.

        Returns:
            sub ``PDBDataFrame``.
        """
        if self._atoms is None:
            if self._ter_line_removed:
                self._atoms = self
            else:
                self._atoms = self.record_names(["ATOM  ", "HETATM"])
                self._ter_line_removed = True
        return self._atoms

    @property
    def coords(self) -> Self:
        """
        Gets the ``x_coord``, ``y_coord``, and ``z_coord`` columns only.
            Use ``pdb_df.coords.values`` to get the underlying
            Numpy array of the coordinates.
        """
        return self.atoms[["x_coord", "y_coord", "z_coord"]]

    @property
    def element_set(self) -> set:
        """Gets the set of element symbols."""
        elements = self.atoms.element_symbol.unique()
        return set([e.strip().upper() for e in elements])

    @property
    @functools.lru_cache()
    def bonds(self) -> dict:
        """Gets the list of bonds. Each bond is represented as a pair of
            ``atom_number`` values.

        Raises:
            ValueError: if the list of ``atom_number`` is not a set.
        """
        return self.get_bonds_by_template()

    def get_bonds_by_distance(
        self,
        single_radii_set: str | None = None,
        need_non_covalent: bool = False,
        non_covalent_cutoff: float = 4.5,
    ) -> dict:
        """Gets all the bonds purely by covalent radii constraints.

        Args:
            single_radii_set (optional): radii sets to use. If ``None``, ``single_C``
                is used as to Cordero (PMID 18478144). Another option is ``single_PA``
                which refers to Pyykkö's studies (PMID 19058281;19856342;15832398,
                and doi:10.1103/PhysRevB.85.024115). Defaults to **None**.
            need_non_covalent (optional): whether non-covalent 'bonding' should be
                included. Defaults to **False**.
            non_covalent_cutoff (optional): distance cutoff for non-covalent 'bonding'.

        Raises:
            ValueError: if the list of ``atom_number`` is not unique or if the
                ``single_radii_set`` is not valid.

        Returns:
            a dictionary of bonds with tuple of ``atom_number`` as keys and bond types
                as values.
        """
        atoms = self.atoms
        if len(atoms) != len(atoms.atom_number.unique()):
            raise ValueError("The 'atom_number' list has repetitive values.")
        if single_radii_set is None:
            single_radii_set = "single_C"
        if single_radii_set not in ["single_C", "single_PA"]:
            message = "The 'single_radii_set' has to be one of "
            message += "'single_C' and 'single_PA', but "
            message += f"{single_radii_set} is provided."
            raise ValueError(message)

        results: dict = {}

        element_set = self.element_set
        if "D" in element_set or "T" in element_set:
            element_set.update({"H"})
        single_bonds, double_bonds, triple_bonds = get_covalent_bond_cutoffs(
            element_set
        )

        # No-bond, triple-bond, and default single bond
        self.use_square_form = False
        non_covalent_cutoff = non_covalent_cutoff**2
        dis_matrix = self.atoms.distance_matrix
        element_matrix = combinations(self.atoms.element_symbol.to_list(), 2)
        atom_number_matrix = combinations(self.atoms.atom_number.to_list(), 2)
        for i, (elements, atom_numbers) in enumerate(
            zip(element_matrix, atom_number_matrix)
        ):
            bond_type: int | float = 0
            dis = dis_matrix[i]
            first_element, second_element = elements
            if dis <= triple_bonds[(first_element, second_element)]:
                bond_type = 3
            elif dis <= double_bonds[(first_element, second_element)]:
                bond_type = 2
            elif dis <= single_bonds[(first_element, second_element)]:
                bond_type = 1
            elif need_non_covalent and dis < non_covalent_cutoff:
                bond_type = 0.5
            atom_number_1, atom_number_2 = atom_numbers
            if bond_type > 0:
                results[(atom_number_1, atom_number_2)] = bond_type

        return results

    def get_bonds_by_template(self) -> dict:
        """Gets covalent bonds based on residue/ligand templates."""
        residue_name_sets = self.residue_name.unique()

        intra_bonds_dict = {}
        for residue_name in residue_name_sets:
            intra_bonds_dict[residue_name] = get_residue_template(
                residue_name=residue_name.strip()
            )

        bonds: dict = {}
        # intro bonds
        all_residues = PDBDataFrame.get_residue_list(self.atoms, include_heteros=True)
        for chain_id, residue_name, residue_number in all_residues:
            residue = (
                self.atoms.chain_ids([chain_id])
                .residue_numbers(residue_number)
                .residue_names([residue_name])
            )
            name_matrix = combinations(residue.atom_name, 2)
            number_matrix = combinations(residue.atom_number, 2)
            for names, numbers in zip(name_matrix, number_matrix):
                first_name, second_name = names
                first_name = first_name.strip()
                second_name = second_name.strip()
                first_number, second_number = numbers
                bond_type = intra_bonds_dict[residue_name].get(
                    (first_name, second_name)
                )
                if bond_type:
                    bonds[(first_number, second_number)] = bond_type[0]
        # inter bonds for peptide bonds
        carb_atoms = self.atoms.atom_names(["C"])
        nitro_atoms = self.atoms.atom_names(["N"])
        if len(carb_atoms) != len(nitro_atoms):
            cn_matrix = PDBDataFrame.get_distance_matrix(carb_atoms, nitro_atoms)
            for carb_index, carb_number in enumerate(carb_atoms.atom_number):
                for nitro_index, nitro_number in enumerate(nitro_atoms.atom_number):
                    if cn_matrix[carb_index, nitro_index] < 2.7889:
                        bonds[(carb_number, nitro_number)] = "SING"
        else:
            for i in range(len(carb_atoms) - 1):
                carb_atom = carb_atoms.iloc[i]
                carb_xyz = (carb_atom.x_coord, carb_atom.y_coord, carb_atom.z_coord)
                carb_number = carb_atom.atom_number
                nitro_atom = nitro_atoms.iloc[i + 1]
                nitro_xyz = (nitro_atom.x_coord, nitro_atom.y_coord, nitro_atom.z_coord)
                nitro_number = nitro_atom.atom_number
                dis = (
                    (carb_xyz[0] - nitro_xyz[0]) ** 2
                    + (carb_xyz[1] - nitro_xyz[1]) ** 2
                    + (carb_xyz[2] - nitro_xyz[2]) ** 2
                )
                if dis < 2.7889:
                    bonds[(carb_number, nitro_number)] = "SING"
        # inter bonds for disulfide bonds
        sulfur_atoms = self.residues.residue_names(["CYS"]).atom_names(["SG"])
        s_dis_matrix = sulfur_atoms.distance_matrix
        s_atom_matrix = combinations(sulfur_atoms.atom_number, 2)
        for i, numbers in enumerate(s_atom_matrix):
            if s_dis_matrix[i] < 9.0:  # 3.0 A is used as cutoff
                bonds[(numbers[0], numbers[1])] = "SING"
        # inter bonds for other pairs
        dis_matrix = self.atoms.distance_matrix
        res_name_matrix = combinations(self.atoms.residue_name, 2)
        atom_number_matrix = combinations(self.atoms.atom_number, 2)
        element_symbol_matrix = combinations(self.atoms.element_symbol, 2)
        element_set = self.element_set
        if "D" in element_set or "T" in element_set:
            element_set.update({"H"})
        single_bonds, double_bonds, triple_bonds = get_covalent_bond_cutoffs(
            element_set
        )
        for i, (names, numbers, symbols) in enumerate(
            zip(res_name_matrix, atom_number_matrix, element_symbol_matrix)
        ):
            first_res = names[0].strip()
            second_res = names[1].strip()
            first_ele = symbols[0]
            second_ele = symbols[1]
            if first_res == "HOH" or second_res == "HOH":
                continue
            if not (first_res in RESIDUE_CODES and second_res in RESIDUE_CODES):
                if dis_matrix[i] <= triple_bonds[(first_ele, second_ele)]:
                    bonds[numbers] = "TRIP"
                elif dis_matrix[i] <= double_bonds[(first_ele, second_ele)]:
                    bonds[numbers] = "DOUB"
                elif dis_matrix[i] <= single_bonds[(first_ele, second_ele)]:
                    bonds[numbers] = "SING"
        return bonds

    @property
    @functools.lru_cache()
    def sequences(self) -> dict[str, str]:
        """
        Gets the sequences for each chain as a dict of ``chain_id``
            as key(s) and ``chain_sequence`` as value(s).
        """
        chain_sequence: dict[str, str] = defaultdict(str)
        for resi_info in self.residue_list:
            chain_sequence[resi_info[0]] += self.RESIDUE_CODES[resi_info[1]]
        return chain_sequence

    @property
    def residue_list(self) -> list[tuple]:
        """
        Gets all residues as a list of tuple
            (``chain_id``, ``residue_name``, ``residue_number``)
        """
        return PDBDataFrame.get_residue_list(self)

    @property
    def backbone(self) -> Self:
        """Gets backbone or N+CA+C+O atoms.

        Returns:
            sub ``PDBDataFrame``
        """
        return self.residues.atom_names(
            ["N", "CA", "C", "O"],
            suppress_warning=True,
        ).element_symbols(["C", "N", "O"])

    @property
    def side_chain(self) -> Self:
        """Gets side chain or NOT N+CA+C+O atoms.

        Returns:
            sub ``PDBDataFrame``.
        """
        return self.residues.atom_names(
            ["N", "CA", "C", "O"],
            invert=True,
            suppress_warning=True,
        ).element_symbols(["C", "N", "O", "S"])

    @property
    def ca_atoms(self) -> Self:
        """Gets the alpha carbon (CA) atoms.

        Returns:
            sub ``PDBDataFrame``.
        """
        return self.residues.atom_names(
            ["CA"],
            suppress_warning=True,
        ).element_symbols(["C"])

    @property
    def heavy_atoms(self) -> Self:
        """Gets the heavy or NOT hydrogen atoms.

        Returns:
            sub ``PDBDataFrame``.
        """
        return self.element_symbols(["H", "D", "T"], invert=True)

    @property
    def hetero_atoms(self) -> Self:
        """Gets the hetero (``HETATM``) atoms.

        Returns:
            sub ``PDBDataFrame``.
        """
        return self.record_names(["HETATM"])

    @property
    def residues(self) -> Self:
        """Gets the residue (``ATOM``) atoms.

        Returns:
            sub ``PDBDataFrame``.
        """
        return self.record_names(["ATOM  "])

    @property
    def water(self) -> Self:
        """Gets all water atoms.

        Returns:
            sub ``PDBDataFrame``.
        """
        return self.atoms.residue_names(["HOH"])

    @property
    def n_atoms(self) -> int:
        """Gets the number of atoms."""
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        """Gets the number of residues."""
        return len(self.residues)

    @property
    def n_chains(self) -> int:
        """Gets the number of chains."""
        ter_lines = self[self.record_name == "TER   "]
        ter_residues = PDBDataFrame.get_residue_list(ter_lines)

        oxt_lines = self[self.atom_name == " OXT"]
        oxt_residues = PDBDataFrame.get_residue_list(oxt_lines)

        n_chain_ids = len(self.chain_id.unique())

        ter_oxt_residues = ter_residues
        for oxt_residue in oxt_residues:
            if oxt_residue not in ter_oxt_residues:
                ter_oxt_residues.append(oxt_residue)

        return max(len(ter_oxt_residues), n_chain_ids)

    @property
    def n_segments(self) -> int:
        """Gets the number of segments."""
        return len(self.atoms.segment_id.unique())

    @property
    def n_models(self) -> int:
        """Gets the number of models."""
        if "nmr_model" in self.columns:
            return len(self.nmr_model.unique())
        return 1

    @property
    @functools.lru_cache()
    def center_of_geometry(self) -> np.ndarray:
        """Gets the center of geometry as a ``(3, )`` ``np.ndarray``."""
        return np.mean(self.coords.values, axis=0)

    @property
    @functools.lru_cache()
    def center_of_mass(self) -> np.ndarray:
        """Gets the center of mass as a ``(3, )`` ``np.ndarray``."""
        masses = self.atoms.get_masses()
        masses = masses / masses.sum()
        return np.sum(self.coords.values * masses[:, None], axis=0)

    @property
    @functools.lru_cache()
    def radius_of_gyration(self) -> float:
        """Gets the radius of gyration"""
        com = self.center_of_mass
        com_t = (com[0], com[1], com[2])  # type: ignore
        dist_to_com = PDBDataFrame.get_distance_matrix(self.atoms, com_t, use_r2=True)
        masses = self.atoms.get_masses()
        masses = masses / masses.sum()
        return np.sum(dist_to_com * masses[:, None], axis=0)[0]

    @functools.lru_cache()
    def get_masses(self) -> np.ndarray:
        """Gets the masses for all atoms in the current dataframe."""
        masses = np.zeros(len(self.atoms), dtype="float32")
        for i, element in enumerate(self.atoms.element_symbol):
            masses[i] = self.ELEMENT_MASSES[element]
        return masses

    @property
    def distance_matrix(self) -> np.ndarray:
        """Gets the distance matrix."""
        return PDBDataFrame.get_distance_matrix(
            self.atoms,
            use_r2=self.use_squared_distance,
            square_form=self.use_square_form,
        )

    def rmsd(
        self,
        other: Self | np.ndarray | None = None,
        align: bool = True,
        weights: list | None = None,
        selection: Self | list | None = None,
    ) -> list | float:
        """Calculates RMSD 1) among sets of coordinates in one ``PDBDataFrame`` with
        multiple ``nmr_model`` s or 2) two sets of coordinates in two ``PDBDataFrames``.

        Args:
            other (optional): the other ``PDBDataFrame`` or (N, 3) ``numpy.ndarray``
                to calculate RMSD against. If ``None``, ``self`` should contain at least
                two sets of coordinates (``nmr_model`` has >= 2 unique values).
                Defaults to **None**.
            align (optional): whether to align the structures before calculating RMSD.
                If ``False``, the ``weights`` and ``selection`` keywords are ignored.
                Defaults to **True**.
            weights (optional): a list of weights for all the atoms in ``selection`` to
                do structure alignment. If ``None``, all coordinates in the
                ``selection`` or ``self`` have the same weights. Defaults to **None**.
            selection (optional): a list of ``atom_number`` s in ``self`` or a
                PDBDataFrame after the filtering methods.
                If ``None``, all coordinates in ``self`` are
                used for structure alignment. Defaults to **None**.

        Returns:
            RMSD or a list of RMSD's.

        Raises:
            ValueError: if dimensionalities mismatch among ``self``, ``other``,
                ``weights``, and ``selection`` if they are not ``None``;
                or ``atom_number`` s in ``self`` are not unique.
            TypeError: if ``other``, ``weights``, and ``selection`` have unsupported
                types.
        """
        result: list = []

        other_coords_list: list = []
        if "nmr_model" in self.columns and len(self.nmr_model.unique()) >= 2:
            all_models = self.nmr_model.unique()
            first_model = self.nmr_models(int(all_models[0]))
            if len(first_model.atom_number.unique()) != len(first_model):
                raise ValueError("'atom_number's in 'self' are not unique.")
            n_atoms = len(first_model)
            if other is None:
                for other_index in all_models[1:]:
                    other_model = self.nmr_models(int(other_index))
                    other_coords_list.append(other_model.coords.values)
        else:
            first_model = self
            if len(first_model.atom_number.unique()) != len(first_model):
                raise ValueError("'atom_number's in 'self' are not unique.")
            n_atoms = len(first_model)
            if other is None:
                message = "'self' has only one set of coordinates and 'other' "
                message += "is not provided."
                raise ValueError(message)
                # message = "'other' is not a PDBDataFrame or np.ndarray instance "
                # message += f"but a {type(other)}."
                # raise ValueError(message)
            elif isinstance(other, type(self)):
                if (self.atom_number != other.atom_number).any():
                    message = "'self' and 'other' have mismatched 'atom_number' values."
                    raise ValueError(message)
                other_coords_list.append(other.coords.values)
            elif isinstance(other, np.ndarray):
                if other.shape != (n_atoms, 3):
                    message = f"'other' shape is {other.shape} but expected to be"
                    message += f" ({n_atoms}, 3)."
                    raise ValueError(message)
                other_coords_list.append(other)
            else:
                raise TypeError(f"Unsupported type {type(other)} for 'other'.")

        first_coords = first_model.coords.values
        align_weights = np.ones(n_atoms)
        if align:
            if selection is not None:
                align_weights = np.zeros(n_atoms)
                if not isinstance(selection, (type(self), list)):
                    message = f"'selection' type is {type(selection)}, but list"
                    message += " or PDBDataFrame is expected."
                    raise TypeError(message)
                n_selection = len(selection)
                if len(selection) > n_atoms:
                    message = f"'selection' length is {len(selection)}, "
                    message += f"larger than the number of atoms {n_atoms}."
                    raise ValueError(message)
                if isinstance(selection, type(self)):
                    selection = self.atom_number.to_list()
                mask = first_model.atom_number.isin(selection)
                align_weights[mask] = 1.0

                if weights is not None:
                    if not isinstance(weights, list):
                        message = f"'weights' type is {type(weights)}, but list"
                        message += " is expected."
                        raise TypeError(message)
                    if len(weights) != n_selection:
                        message = f"'weights' length is {len(weights)}, "
                        message += "not equal to the selection length"
                        message += f" {n_selection}."
                        raise ValueError(message)
                    align_weights[mask] = np.array(weights)
            elif weights is not None:
                if not isinstance(weights, list):
                    message = f"'weights' type is {type(weights)}, but list"
                    message += " is expected."
                    raise TypeError(message)
                if len(weights) != n_atoms:
                    message = f"'weights' length is {len(weights)}, "
                    message += "not equal to the number of atoms"
                    message += f" {n_atoms} without 'selection'."
                    raise ValueError(message)
                align_weights = np.array(weights)

        first_centered = False
        for other_coords in other_coords_list:
            if align:
                if not first_centered:
                    first_coords = first_coords - first_model.center_of_geometry
                    first_centered = True
                other_cog = np.mean(other_coords, axis=0)
                other_coords = other_coords - other_cog
                rot = Rotation.align_vectors(
                    first_coords, other_coords, weights=align_weights
                )[0]
                other_coords = rot.apply(other_coords)

            rms = np.sqrt(np.sum((first_coords - other_coords) ** 2) / n_atoms)
            result.append(rms)

        if len(result) == 1:
            return result[0]

        return result

    def record_names(self, names: list[str], invert: bool = False) -> Self:
        """Filter by ``record_name``.

        Args:
            names (required): a list of ``record_name`` s.
            invert (optional): whether to invert the selection. Defaults to **False**.
        Returns:
            sub ``PDBDataFrame``
        """
        names = [name.strip().upper() for name in names]
        if self.pdb_format.upper() == "PDB":  # type: ignore
            names = [name.ljust(6) for name in names]
        if invert:
            return self[~self.record_name.isin(names)]
        return self[self.record_name.isin(names)]

    def atom_numbers(
        self,
        numbers: list[int] | int,
        relation: str | None = None,
        invert: bool = False,
    ) -> Self:
        """Filter by ``atom_number``.

        Args:
            numbers (required): one or a list of ``atom_number`` s.
            relation (optional): ``atom_number`` relationship to ``numbers``.
                If ``numbers`` is an integer, it has to be one of ``<``, ``<=``, ``=``,
                ``>=``, and ``>``. If ``None``, ``<=`` is used. Ignored if a list of
                integers are provided to ``numbers``. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            numbers, "atom_number", relation=relation, invert=invert
        )

    def atom_names(
        self,
        names: list[str],
        names_2c: list[str] | None = None,
        invert: bool = False,
        suppress_warning: bool = False,
    ) -> Self:
        """Filter by ``atom_name``.

        Args:
            names (required): a list of ``atom_name`` s whose ``element_symbols`` have
                only one character. Atoms in common residues and ligands should be
                provide here like ``C, H, O, N, S, P, F``.
            names_2c (optional): a list of ``atom_name`` s whose ``element_symbols``
                have two characters like ion (``FE``) and chloride (``CL``).
                Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            suppress_warning: whether to suppress the warning message about
                possible conflicts between ``names`` and ``names_2c``.
                Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        atom_name_strings = [name.strip().upper() for name in names]
        if self.pdb_format.upper() == "PDB":  # type: ignore
            atom_name_strings = []
            for name in names:
                if len(name) == 4:
                    atom_name_strings.append(name)
                else:
                    atom_name_strings.append(f" {name}".ljust(4))
                if (
                    len(name) == 2
                    and f"{name[0]}{name[1].lower()}" in self.ELEMENT_MASSES
                ):
                    if f" {name[0]}" not in self.ELEMENT_MASSES:
                        atom_name_strings.append(f"{name}".ljust(4))
                        # eg 'MG' where ' M' is not a legal ele in self.ELEMENT_MASSES.
                        continue
                    if suppress_warning:
                        continue
                    message = f"Atom name {name} is an atom of element {name[0]} "
                    message += f"but not element {name[0]}{name[1].lower()}."
                    message += "If you want the latter, put it in the 'names_2c' list."
                    warnings.warn(
                        message,
                        RuntimeWarning,
                        stacklevel=2,
                    )
            if names_2c is not None:
                for name in names_2c:
                    atom_name_strings.append(f"{name}".ljust(4))

        if invert:
            return self[~self.atom_name.isin(atom_name_strings)]
        return self[self.atom_name.isin(atom_name_strings)]

    def alt_locs(self, locs: list[str], invert: bool = False) -> Self:
        """Filter by ``alt_loc``.

        Args:
            locs (required): a list of ``alt_loc`` s.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub PDBDataFrame
        """
        if invert:
            return self[~self.alt_loc.isin(locs)]
        return self[self.alt_loc.isin(list(locs) + [" "])]

    def residue_names(self, names: list[str], invert: bool = False) -> Self:
        """Filter by ``residue_names``.

        Args:
            names (required): a list of ``residue_name`` s
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        names = [name.strip().upper() for name in names]
        if self.pdb_format.upper() == "PDB" and self.is_chimera:  # type: ignore
            names = [(name + " ").upper().rjust(4) for name in names]
        if invert:
            return self[~self.residue_name.isin(names)]
        return self[self.residue_name.isin(names)]

    def chain_ids(self, ids: list[str], invert: bool = False) -> Self:
        """Filter by ``chain_id``.

        Args:
            ids (required): a list of ``chain_id`` s.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        if invert:
            return self[~self.chain_id.isin(ids)]
        return self[self.chain_id.isin(ids)]

    def residue_numbers(
        self,
        numbers: list[int] | int,
        relation: str | None = None,
        invert: bool = False,
    ) -> Self:
        """Filter by ``residue_number``.

        Args:
            numbers (required): one or a list of ``residue_number`` s.
            relation (optional): ``residue_number`` relationship to ``numbers``.
                If ``numbers`` is an integer, it has to be one of ``<``, ``<=``, ``=``,
                ``>=``, and ``>``. If ``None``, '<=' is used. Ignored if a list of
                integers are provided to ``numbers``. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            numbers, "residue_number", relation=relation, invert=invert
        )

    def insertions(self, codes: list[str], invert: bool = False) -> Self:
        """Filter by ``insertion``.

        Args:
            codes (required): a list of ``insertion`` codes.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        if invert:
            return self[~self.insertion.isin(codes)]
        return self[self.insertion.isin(codes)]

    def x_coords(
        self,
        value: float,
        relation: str | None = None,
        invert: bool = False,
        epsilon: float = 0.01,
    ) -> Self:
        """Filter by ``x_coord``.

        Args:
            value (required): value to select ``x_coord`` s.
            relation (optional): ``x_coord`` relationship to ``value``.
                It has to be one of ``'<'``, ``'<='``, ``'='``, ``'>='``, and ``'>'``.
                If ``None``, ``'<='`` is used. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            epsilon (optional): atoms ``abs(x_coord - value)`` <= ``epsilon``
                are selected when ``invert`` = ``False``. Defaults to **0.01**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            value, "x_coord", relation=relation, invert=invert, epsilon=epsilon
        )

    def y_coords(
        self,
        value: float,
        relation: str | None = None,
        invert: bool = False,
        epsilon: float = 0.01,
    ) -> Self:
        """Filter by ``y_coord``.

        Args:
            value (required): value to select ``y_coord`` s.
            relation (optional): ``y_coord`` relationship to ``value``.
                It has to be one of ``'<'``, ``'<='``, ``'='``, ``'>='``, and ``'>'``.
                If ``None``, ``'<='`` is used. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            epsilon (optional): atoms ``abs(y_coord - value)`` <= ``epsilon``
                are selected when ``invert`` = ``False``. Defaults to **0.01**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            value, "y_coord", relation=relation, invert=invert, epsilon=epsilon
        )

    def z_coords(
        self,
        value: float,
        relation: str | None = None,
        invert: bool = False,
        epsilon: float = 0.01,
    ) -> Self:
        """Filter by ``z_coord``.

        Args:
            value (required): value to select ``z_coord`` s.

            relation (optional): ``z_coord`` relationship to ``value``.
                It has to be one of ``'<'``, ``'<='``, ``'='``, ``'>='``, and ``'>'``.
                If ``None``, ``'<='`` is used. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            epsilon (optional): atoms ``abs(z_coord - value)`` <= ``epsilon``
                are selected when ``invert`` = ``False``. Defaults to **0.01**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            value, "z_coord", relation=relation, invert=invert, epsilon=epsilon
        )

    def occupancies(
        self,
        value: float,
        relation: str | None = None,
        invert: bool = False,
        epsilon: float = 0.01,
    ) -> Self:
        """Filter by ``occupancy``.

        Args:
            value (required): value to select ``occupancy`` s.
            relation (optional): ``occupancy`` relationship to ``value``.
                It has to be one of ``'<'``, ``'<='``, ``'='``, ``'>='``, and ``'>'``.
                If ``None``, ``'<='`` is used. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            epsilon (optional): atoms ``abs(occupancy - value)`` <= ``epsilon``
                are selected when ``invert`` = ``False``. Defaults to **0.01**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            value, "occupancy", relation=relation, invert=invert, epsilon=epsilon
        )

    def b_factors(
        self,
        value: float,
        relation: str | None = None,
        invert: bool = False,
        epsilon: float = 0.01,
    ) -> Self:
        """Filter by ``b_factor``.

        Args:
            value (required): value to select ``b_factor`` s.
            relation (optional): ``b_factor`` relationship to ``value``.
                It has to be one of ``'<'``, ``'<='``, ``'='``, ``'>='``, and ``'>'``.
                If ``None``, ``'<='`` is used. Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            epsilon (optional): atoms ``abs(b_factor - value)`` <= ``epsilon``
                are selected when ``invert`` = ``False``. Defaults to **0.01**.

        Returns:
            sub ``PDBDataFrame``
        """
        return self._filter_num_col(
            value, "b_factor", relation=relation, invert=invert, epsilon=epsilon
        )

    def segment_ids(self, ids: list[str], invert: bool = False) -> Self:
        """Filter by ``segment_id``.

        Args:
            ids (required): a list of ``segment_id`` s.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        ids = [i.strip() for i in ids]
        if self.pdb_format.upper() == "PDB":  # type: ignore
            ids = [i.ljust(4) for i in ids]
        if invert:
            return self[~self.segment_id.isin(ids)]
        return self[self.segment_id.isin(ids)]

    def element_symbols(self, symbols: list[str], invert: bool = False) -> Self:
        """Filter by ``element_symbol``.

        Args:
            symbols (required): a list of ``element_symbol`` s.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        symbols = [symbol.strip().upper() for symbol in symbols]
        if self.pdb_format.upper() == "PDB":  # type: ignore
            symbols = [symbol.rjust(2) for symbol in symbols]
        if invert:
            return self[~self.element_symbol.isin(symbols)]
        return self[self.element_symbol.isin(symbols)]

    def charges(self, charges: list[str], invert: bool = False) -> Self:
        """Filter by ``charge``.

        Args:
            charges (required): a list of ``charge`` s.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``

        Notes: ``charge`` is ``2-char`` string in the PDB specifications.
        """
        charges = [charge.strip() for charge in charges]
        if self.pdb_format.upper() == "PDB":  # type: ignore
            charges = [charge.rjust(2) for charge in charges]
        if invert:
            return self[~self.charge.isin(charges)]
        return self[self.charge.isin(charges)]

    def nmr_models(
        self, models: list[int] | int, relation: str | None = None, invert: bool = False
    ) -> Self:
        """Filter by ``nmr_model``.

        Args:
            models (required): one or a list of ``nmr_model`` ids.
            relation (optional): ``nmr_model`` relationship to ``models``.
                If ``models`` is an integer, it has to be one of ``'<'``, ``'<='``,
                ``'='``, ``'>='``, and ``'>'``. If ``None``, ``'<='`` is used.
                Ignored if a list of integers are provided to ``models``.
                Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        if "nmr_model" in self.columns:
            return self._filter_num_col(
                models, "nmr_model", relation=relation, invert=invert
            )
        return self

    def distances(
        self,
        other: np.ndarray | Self | Iterable,
        cut_off: float = np.inf,
        to: str | None = None,
        invert: bool = False,
    ) -> Self:
        """Filter by ``distance`` to a reference point or group of atoms.

        Args:
            other: the other group's coordinate(s).
            cut_off: the distance cutoff to filter.
            to: if ``other`` is a group atoms, using which method to determine
                whether the atoms meet the cut_off distance. If None, ``COM`` or
                center of mass is used if ``other`` is ``PDBDataFrame``, and ``COG``
                or center of geometry is used if ``other`` is ``np.ndarray`` or
                ``Iterable``. The following are allowed:

                    ``com``, ``center of mass``, ``center_of_mass``:
                        use the center of mass for the ``other``.
                    ``cog``, ``center of geometry``, ``center_of_geometry``:
                        use the center of geometry for the ``other``.
                    ``all``:
                        whether all the pair-distances meet the ``cut_off``
                        distance criteria.
                    ``any``:
                        whether any of the pair-distances meets the ``cut_off``
                        distance criteria.
            invert: whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
        if not (isinstance(other, type(self)) or isinstance(other, Iterable)):
            message = "Only 'PDBDataFrame', 'np.ndarray', and 'Iterable' types "
            message += f"are supported, not {type(other)} for 'other'."
            raise TypeError(message)

        if isinstance(other, Iterable) and not isinstance(other, type(self)):
            other = np.asanyarray(other, dtype="float32")
            if other.shape != (3,):
                if len(other.shape) != 2 or other.shape[1] != 3:
                    message = "An 'Iterable' input should have a (N, 3) or (3,) shape, "
                    message += f"not '{other.shape}' in 'other'."
                    raise TypeError(message)
            else:
                other = other.reshape(1, 3)

        cut_off = cut_off**2

        allowed_tos = [
            "com",
            "center of mass",
            "center_of_mass",
            "cog",
            "center of geometry",
            "center_of_geometry",
            "any",
            "all",
        ]
        if to is None or to.lower() not in allowed_tos:
            if to is not None:
                message = "Only center of mass (COM), center of geometry (COG) "
                message += f"'all', or 'any' is supported. '{to}' is reset to "
                if isinstance(other, type(self)):
                    message += "'COM' for 'to'."
                else:
                    message += "'COG' for 'to'."
                warnings.warn(
                    message,
                    RuntimeWarning,
                    stacklevel=2,
                )
            to = "COG"
            if isinstance(other, type(self)):
                to = "COM"
        elif to.lower() in [
            "com",
            "center of mass",
            "center_of_mass",
        ] and not isinstance(other, type(self)):
            message = "Only center of geometry (COG), 'all', or 'any' is supported "
            message += f"if 'other' is 'Iterable'. '{to}' is reset to 'COG' for 'to'."
            to = "COG"
            warnings.warn(
                message,
                RuntimeWarning,
                stacklevel=2,
            )

        other_data: np.ndarray | tuple | Self = np.asanyarray(other)
        if isinstance(other, type(self)):
            if to.lower() in ["com", "center of mass", "center_of_mass"]:
                to = "COM"
                other_data = np.asanyarray(other.center_of_mass)
                other_data = (other_data[0], other_data[1], other_data[2])
            elif to.lower() in ["cog", "center of geometry", "center_of_geometry"]:
                to = "COG"
                other_data = np.asanyarray(other.center_of_geometry)
                other_data = (other_data[0], other_data[1], other_data[2])
            else:  # 'any' or 'all'
                to = to.upper()
                other_data = other

        else:
            if to.lower() in ["cog", "center of geometry", "center_of_geometry"]:
                to = "COG"
                other_data = np.mean(other, axis=0)
                other_data = (other_data[0], other_data[1], other_data[2])
            else:  # 'any' or 'all'
                to = to.upper()
                other_data = tuple(
                    [
                        (other[i, 0], other[i, 1], other[i, 2])
                        for i in range(other.shape[0])
                    ]
                )

        distance_matrix = PDBDataFrame.get_distance_matrix(
            self,
            other_data=other_data,
        )
        mask = distance_matrix <= cut_off
        if to in ["COG", "COM"]:
            if invert:
                return self.atoms[~mask[:, 0]]
            return self.atoms[mask[:, 0]]

        if (to == "ALL" and not invert) or (to == "ANY" and invert):
            return self.atoms[mask.all(axis=1)]
        return self.atoms[mask.any(axis=1)]

    def _filter_num_col(
        self,
        value: int | float | list[int],
        num_col_name: str,
        relation: str | None = None,
        invert: bool = False,
        epsilon: float | int = 0.01,
        suppress_warning: bool = False,
    ) -> Self:
        """Generic function to do filter by a numerical column.

        Args:
            value (required): value(s) to select by the column given by the
                ``num_col_name`` input.
            num_col_name (required): one of ``atom_number``, ``residue_number``,
                ``x_coord``, ``y_coord``, or ``z_coord``, ``occupancy``, ``b_factor``,
                and  ``nmr_model``. Note: the ``charge`` column is not numerical by
                ``PDB`` format.
            relation (optional): ``x/y/z_coord`` relationship to ``value``.
                It has to be one of ``'<'``, ``'<='``, ``'='``, ``'>='``, and ``'>'``.
                If ``None``, ``'<='`` is used.
                Ignored if a list of integers are provided to ``value``.
                Defaults to **None**.
            invert (optional): whether to invert the selection. Defaults to **False**.
            epsilon (optional): atoms ``abs``(``num_col_value`` - ``value``) <=
                ``epsilon`` are selected when ``invert`` = ``False`` and
                ``relation`` = ``'='``. Ignored if a list of integers are provided to
                ``value``. Defaults to **0.01**.
            suppress_warning (optional): whether to suppress warnings.
                Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``

        Raises:
            ValueError: if xyz not in [``atom_number``, ``residue_number``,
                ``x_coord``, ``y_coord``, or ``z_coord``,
                ``occupancy``, ``b_factor``, ``nmr_model``] or ``relation`` not in
                [``'<'``, ``'<='``, ``'='``, ``'>='``, ``'>'``] when
                selecting on float cols.
        """
        allowed_num_col_names = [
            "atom_number",
            "residue_number",
            "x_coord",
            "y_coord",
            "z_coord",
            "occupancy",
            "b_factor",
            "nmr_model",
        ]
        if num_col_name not in allowed_num_col_names:
            message = f"Only '{allowed_num_col_names}' are allowed in 'num_col_name' "
            message += f"but {num_col_name} was put."
            raise ValueError(message)

        allowed_relations = ["<", "<=", "=", ">=", ">"]
        if relation is None and not isinstance(value, list):
            if isinstance(value, float):
                relation = "<="
            elif isinstance(value, int):
                relation = "="
            else:
                message = "Only 'int', 'float', or 'list[int]' are allowed in 'value' "
                message += f"but {type(value)} was put."
                raise ValueError(message)
        elif not isinstance(value, list) and relation not in allowed_relations:
            message = f"Only '{allowed_relations}' are allowed in 'relation' "
            message += f"but {relation} was put."
            raise ValueError(message)
        elif isinstance(value, list):
            for v in value:
                if not isinstance(v, int):
                    message = "Only 'int' is allowed in 'value' if it is a list, "
                    message += f"but {type(v)} was in {value}."
                    raise ValueError(message)
            if relation is not None:
                relation = None
                if not suppress_warning:
                    message = "'relation' is ignored "
                    message += "when a list is provided to 'value'."
                    warnings.warn(
                        message,
                        RuntimeWarning,
                        stacklevel=2,
                    )

        if relation == "<":
            if invert:
                return self.atoms[self.atoms[num_col_name].values >= value]
            return self.atoms[self.atoms[num_col_name].values < value]
        elif relation == "<=":
            if invert:
                return self.atoms[self.atoms[num_col_name].values > value]
            return self.atoms[self.atoms[num_col_name].values <= value]
        elif relation == "=":
            if invert:
                return self.atoms[
                    np.abs(self.atoms[num_col_name].values - value) >= epsilon
                ]
            return self.atoms[np.abs(self.atoms[num_col_name].values - value) < epsilon]
        elif relation == ">=":
            if invert:
                return self.atoms[self.atoms[num_col_name].values < value]
            return self.atoms[self.atoms[num_col_name].values >= value]
        elif relation == ">":
            if invert:
                return self.atoms[self.atoms[num_col_name].values <= value]
            return self.atoms[self.atoms[num_col_name].values > value]

        # relation is None -> a list of numbers
        if invert:
            return self.atoms[~np.isin(self.atoms[num_col_name].values, value)]
        return self.atoms[np.isin(self.atoms[num_col_name].values, value)]

    @classmethod
    def get_residue_list(
        cls, pdb_df: Self, include_heteros: bool = False
    ) -> list[tuple]:
        """Gets the list of residues given a ``PDBDataFrame`` object.

        Args:
            pdb_df (required): a ``PDBDataFrame`` object.
            include_heteros (optional): whether to include hetero ligands.
                Defaults to **False**.

        Returns:
            a list of residues as (``chain_id``, ``residue_name``, ``residue_number``).
        """
        all_residues: list[tuple] = []
        for chain, residue_name, residue_number in zip(
            pdb_df["chain_id"], pdb_df["residue_name"], pdb_df["residue_number"]
        ):
            if not include_heteros and residue_name not in pdb_df.RESIDUE_CODES:
                continue
            residue = (chain, residue_name, residue_number)
            if len(all_residues) == 0:
                all_residues.append(residue)
            elif residue != all_residues[-1]:
                all_residues.append(residue)
        return all_residues

    @classmethod
    @functools.lru_cache()
    def get_distance_matrix(
        cls,
        pdb_df: Self,
        other_data: Self | tuple | None = None,
        use_r2: bool = True,
        square_form: bool = False,
    ) -> np.ndarray:
        """Calculates the distance matrix given a ``PDBDataFrame`` object and
        (optional) reference data.

        Args:
            pdb_df (required): a ``PDBDataFrame`` object.
            other_data (optional): the coordinates of to calculate the distances
                against. Defaults to **None**.
            use_r2 (optional): whether to use r^2 or r for distance matrix.
                Defaults to **True**.
            square_form (optional): whether to output a square form of the
                density matrix. If two ``PDBDataFrame``s are different or
                ``other_data`` is not a ``PDBDataFrame``, ``square_form`` is ignored.
                Defaults to **False**.


        Returns:
            distance matrix (squared or condensed form)

        Raises:
            ValueError: if ``other_data`` is not of ``PDBDataFrame|tuple|None`` type or
                wrong shape if it is a ``tuple``.
        """
        cols = ["x_coord", "y_coord", "z_coord"]

        if other_data is None or (
            isinstance(other_data, cls) and pdb_df.atoms == other_data.atoms
        ):
            if not square_form:
                if use_r2:
                    return pdist(pdb_df.atoms[cols].values, "sqeuclidean")
                return pdist(pdb_df.atoms[cols].values)
            if use_r2:
                return squareform(pdist(pdb_df.atoms[cols].values, "sqeuclidean"))
            return squareform(pdist(pdb_df.atoms[cols].values))

        elif isinstance(other_data, cls):
            if use_r2:
                return cdist(
                    pdb_df.atoms[cols].values,
                    other_data.atoms[cols].values,
                    "sqeuclidean",
                )
            return cdist(
                pdb_df.atoms[cols].values, other_data.atoms[cols].values, "euclidean"
            )

        elif isinstance(other_data, tuple):
            other_array = np.asanyarray(other_data)
            if other_array.shape != (3,):
                if not (len(other_array.shape) == 2 and other_array.shape[1] == 3):
                    message = "'other_data' expects a shape of (N, 3) or (3,) "
                    message += f"if given a tuple, but {other_array.shape} was given."
                    raise ValueError(message)
            else:
                other_array = other_array.reshape(1, 3)

            if use_r2:
                return cdist(pdb_df.atoms[cols].values, other_array, "sqeuclidean")
            return cdist(pdb_df.atoms[cols].values, other_array, "euclidean")
        raise ValueError("'other_data' type has to be tuple or PDBDataFrame or None.")
