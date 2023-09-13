# pdbx2df
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/pdbx2df
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

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from scipy.spatial.distance import cdist, pdist, squareform  # type: ignore
from typing_extensions import Self

from .constants import AMINO_ACIDS, ELEMENT_MASSES

"""dict[str, str], turn 3-, 2-, and 1-letter residue codes to 1-letter codes."""
RESIDUE_CODES = AMINO_ACIDS


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
    >>> from pdbx2df import read_pdb, PDBDataFrame
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
        self.pdb_format: str | None = pdb_format
        if self.pdb_format is None:
            self.pdb_format = "PDB"
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
        try:
            if len(self.head(1).residue_name[0]):
                self._is_chimera = True
        except AttributeError:
            pass  # 'residue_name' not in self.coords
        except IndexError:
            pass  # empty

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
                self._atoms = self[self.record_name.isin(["ATOM  ", "HETATM"])]
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
        return np.mean(self.coords.values * masses[:, None], axis=0)

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

    def record_names(self, names: list[str], invert: bool = False) -> Self:
        """Filter by ``record_name``.

        Args:
            names (required): a list of ``record_name`` s.

            invert (optional): whether to invert the selection. Defaults to **False**.

        Returns:
            sub ``PDBDataFrame``
        """
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
        atom_name_strings: list[str] = []
        for name in names:
            if len(name) == 4:
                atom_name_strings.append(name)
            else:
                atom_name_strings.append(f" {name}".ljust(4))
            if len(name) == 2 and f"{name[0]}{name[1].lower()}" in self.ELEMENT_MASSES:
                if f" {name[0]}" not in self.ELEMENT_MASSES:
                    atom_name_strings.append(f"{name}".ljust(4))
                    # eg 'MG' where ' M' is not a legal element in self.ELEMENT_MASSES.
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
        if invert:
            return self[~self.residue_name.str.strip().isin(names)]
        return self[self.residue_name.str.strip().isin(names)]

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
        symbols = [symbol.upper().rjust(2) for symbol in symbols]
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
    def get_residue_list(cls, pdb_df: Self) -> list[tuple]:
        """Gets the list of residues given a ``PDBDataFrame`` object.

        Args:
            pdb_df (required): a ``PDBDataFrame`` object.

        Returns:
            a list of residues as (``chain_id``, ``residue_name``, ``residue_number``).
        """
        all_residues: list[tuple] = []
        for chain, residue_name, residue_number in zip(
            pdb_df["chain_id"], pdb_df["residue_name"], pdb_df["residue_number"]
        ):
            if residue_name not in pdb_df.RESIDUE_CODES:
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
