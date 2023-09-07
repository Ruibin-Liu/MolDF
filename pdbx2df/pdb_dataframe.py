# pdbx2df
# Author: Ruibin Liu <ruibinliuphd@gmail.com>
# License: MIT
# Code Repository: https://github.com/Ruibin-Liu/pdbx2df
"""PdbDataFrame as a subclass of Pandas DataFrame"""
from __future__ import annotations

import functools
from functools import _lru_cache_wrapper

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from scipy.spatial.distance import pdist, squareform  # type: ignore

from .constants import AMINO_ACIDS, ELEMENT_MASSES


class PdbDataFrame(pd.DataFrame):
    """Pandas DataFrame with extended attributes and methods for PDB data."""

    _metadata = [
        # "_n_atoms",
        # "_n_residues",
        # "_n_chains",
        # "_n_models",
        # "_n_segments",
        # "_residues",
        # "_COM",
        # "_COG",
        # "_sequence",
        "_use_squared_distance",
    ]

    def __init__(
        self,
        *args,
        pdb_format: str | None = None,
        use_squared_distance: bool = True,
        **kwargs,
    ):
        """Init the PdbDataFrame class like a Pandas DataFrame with some extra keywords.

        Args:
            pdb_format (str|None; defaults to None): PDB format in the underlying provided data.
                If None, "PDB" is assumed.
            use_squared_distance (bool; defaults to True): whether to use squared distance
                when calculating distance matrix.
        """
        super().__init__(*args, **kwargs)
        self.pdb_format: str | None = pdb_format
        if self.pdb_format is None:
            self.pdb_format = "PDB"
        self._use_squared_distance: bool = use_squared_distance

    @property
    def _constructor(self):
        return PdbDataFrame

    def __hash__(self):
        """Use head X coords to hash; for distance matrix calculation cache."""
        return hash(tuple(self.head().x_coord))

    def __eq__(self, other):
        """Use head X coords to compare; for distance matrix calculation cache."""
        return hash(tuple(self.head().x_coord)) == hash(tuple(other.head().x_coord))

    @property
    def atoms(self):
        """Returns 'ATOM' and 'HETATM' entries.

        Returns:
            PdbDataFrame: sub DataFrame.
        """
        return self[self.record_name.str.strip().isin(["ATOM", "HETATM"])]

    @property
    def backbone(self):
        """Returns backbone or N+CA+C+O atoms.

        Returns:
            PdbDataFrame: sub DataFrame.
        """
        return self.atom_names(["N", "CA", "C", "O"])

    @property
    def side_chain(self):
        """Returns side chain or NOT N+CA+C+O atoms.

        Returns:
            PdbDataFrame: sub DataFrame.
        """
        return self.atom_names(["N", "CA", "C", "O"], invert=True)

    @property
    def ca_atoms(self):
        """Returns CA atoms.

        Returns:
            PdbDataFrame: sub DataFrame.
        """
        return self.atom_names(["CA"])

    @property
    def heavy_atoms(self):
        """Returns heavy or NOT hydrogen atoms.

        Returns:
            PdbDataFrame: sub DataFrame.
        """
        return self.element_symbols(["H", "D", "T"], invert=True)

    @property
    def n_atoms(self) -> int:
        """Returns number of atoms."""
        return len(self[self.record_name.str.strip().isin(["ATOM", "HETATM"])])

    @property
    def n_residues(self) -> int:
        """Returns number of residues."""
        return len(self.residues)

    @property
    def n_chains(self) -> int:
        """Returns number of chains."""
        ter_lines = self[self.record_name.str.strip() == "TER"]
        ter_residues = PdbDataFrame.get_residues(ter_lines)

        oxt_lines = self[self.atom_name.str.strip() == "OXT"]
        oxt_residues = PdbDataFrame.get_residues(oxt_lines)

        n_chain_ids = len(self.chain_id.unique())

        ter_oxt_residues = ter_residues
        for oxt_residue in oxt_residues:
            if oxt_residue not in ter_oxt_residues:
                ter_oxt_residues.append(oxt_residue)

        return max(len(ter_oxt_residues), n_chain_ids)

    @property
    def n_segments(self) -> int:
        """Returns number of segments."""
        return len(self.atoms.segment_id.unique())

    @property
    def n_models(self) -> int:
        """Returns number of models."""
        if "nmr_model" in self.columns:
            return len(self.nmr_model.unique())
        return 1

    @property
    def residues(self) -> list[tuple]:
        """Returns all residues as a tuple of (chain_id, residue_name, residue_number)"""
        return PdbDataFrame.get_residues(self)

    @property
    @functools.lru_cache()
    def sequence(self) -> str:
        """Returns the sequence."""
        return "".join([AMINO_ACIDS[i[1]] for i in self.residues])

    @property
    @functools.lru_cache()
    def center_of_geometry(self) -> tuple[float, float, float]:
        """Returns the center of geometry."""
        return (
            np.mean(self.atoms.x_coord),
            np.mean(self.atoms.y_coord),
            np.mean(self.atoms.z_coord),
        )

    @property
    @functools.lru_cache()
    def center_of_mass(self) -> tuple[float, float, float]:
        """Returns the center of mass."""
        masses = self.atoms.get_masses()
        total_mass = sum(masses)
        return (
            sum(np.array(self.atoms.x_coord) * np.array(masses)) / total_mass,
            sum(np.array(self.atoms.y_coord) * np.array(masses)) / total_mass,
            sum(np.array(self.atoms.z_coord) * np.array(masses)) / total_mass,
        )

    def get_masses(self) -> list:
        """Returns the masses for all atoms."""
        masses = [ELEMENT_MASSES[atom.strip()] for atom in self.element_symbol]
        return masses

    @property
    def use_squared_distance(self) -> bool:
        """Returns whether R or R2 is used in distance matrix calculations."""
        return self._use_squared_distance

    @property
    def distance_matrix(self) -> np.ndarray:
        """Returns the distance matrix."""
        return PdbDataFrame.get_distance_matrix(
            self.atoms,
            use_r2=self.use_squared_distance,
        )  # type: ignore

    def record_names(self, names: list[str], invert: bool = False):
        """Filter by 'record_name'.

        Args:
            names (list[str]): a list of 'record_name's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.record_name.str.strip().isin(names)]
        return self[self.record_name.str.strip().isin(names)]

    def atom_numbers(self, numbers: list[int], invert: bool = False):
        """Filter by 'atom_number'.

        Args:
            numbers (list[int]): a list of 'atom_number's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.atom_number.isin(numbers)]
        return self[self.atom_number.isin(numbers)]

    def atom_names(self, names: list[str], invert: bool = False):
        """Filter by 'atom_name'.

        Args:
            names (list[str]): a list of 'atom_name's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.atom_name.str.strip().isin(names)]
        return self[self.atom_name.str.strip().isin(names)]

    def alt_locs(self, locs: list[str], invert: bool = False):
        """Filter by 'atom_number'.

        Args:
            locs (list[str]): a list of 'alt_loc's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.alt_loc.str.strip().isin(locs)]
        return self[self.alt_loc.str.strip().isin(list(locs) + [""])]

    def residue_names(self, names: list[str], invert: bool = False):
        """Filter by 'residue_names'.

        Args:
            names (list[str]): a list of 'residue_name's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.residue_name.str.strip().isin(names)]
        return self[self.residue_name.str.strip().isin(names)]

    def chain_ids(self, ids: list[str], invert: bool = False):
        """Filter by 'chain_id'.

        Args:
            ids (list[str]): a list of 'chain_id's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.chain_id.isin(ids)]
        return self[self.chain_id.isin(ids)]

    def residue_numbers(self, numbers: list[int], invert: bool = False):
        """Filter by 'residue_number'.

        Args:
            numbers (list[int]): a list of 'residue_number's
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.residue_number.isin(numbers)]
        return self[self.residue_number.isin(numbers)]

    def insertions(self, codes: list[str], invert: bool = False):
        """Filter by 'insertion'.

        Args:
            codes (list[str]): a list of 'insertion' codes
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.insertion.isin(codes)]
        return self[self.insertion.isin(codes)]

    def x_coords(self, value: float, invert: bool = False, epsilon: float = 0.01):
        """Filter by 'x_coord'.

        Args:
            value (float): value to select 'x_coord's.
            invert (bool; defaults to False): whether to invert the selection.
            epsilon (float; defaults to 0.001): atoms |x_coord - value| <= epsilon are selected when invert=False.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[abs(self.x_coord - value) > epsilon]
        return self[abs(self.x_coord - value) <= epsilon]

    def y_coords(self, value: float, invert: bool = False, epsilon: float = 0.01):
        """Filter by 'y_coord'.

        Args:
            value (float): value to select 'y_coord's.
            invert (bool; defaults to False): whether to invert the selection.
            epsilon (float; defaults to 0.001): atoms |y_coord - value| <= epsilon are selected when invert=False.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[abs(self.y_coord - value) > epsilon]
        return self[abs(self.y_coord - value) <= epsilon]

    def z_coords(self, value: float, invert: bool = False, epsilon: float = 0.01):
        """Filter by 'z_coord'.

        Args:
            value (float): value to select 'z_coord's.
            invert (bool; defaults to False): whether to invert the selection.
            epsilon (float; defaults to 0.001): atoms |z_coord - value| <= epsilon are selected when invert=False.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[abs(self.y_coord - value) > epsilon]
        return self[abs(self.y_coord - value) <= epsilon]

    def occupancies(self, value: float, invert: bool = False, epsilon: float = 0.01):
        """Filter by 'occupancy'.

        Args:
            value (float): value to select 'occupancy's.
            invert (bool; defaults to False): whether to invert the selection.
            epsilon (float; defaults to 0.001): atoms |occupancy - value| <= epsilon are selected when invert=False.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[abs(self.occupancy - value) > epsilon]
        return self[abs(self.occupancy - value) <= epsilon]

    def b_factors(self, value: float, invert: bool = False, epsilon: float = 0.01):
        """Filter by 'b_factor'.

        Args:
            value (float): value to select 'b_factor's.
            invert (bool; defaults to False): whether to invert the selection.
            epsilon (float; defaults to 0.001): atoms |b_factor - value| <= epsilon are selected when invert=False.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[abs(self.b_factor - value) > epsilon]
        return self[abs(self.b_factor - value) <= epsilon]

    def segment_ids(self, ids: list[str], invert: bool = False):
        """Filter by 'segment_id'.

        Args:
            ids (list[str]): a list of 'segment_id's.
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.segment_id.isin(ids)]
        return self[self.segment_id.isin(ids)]

    def element_symbols(self, symbols: list[str], invert: bool = False):
        """Filter by 'element_symbol'.

        Args:
            symbols (list[str]): a list of 'element_symbol's.
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.element_symbol.str.strip().isin(symbols)]
        return self[self.element_symbol.str.strip().isin(symbols)]

    def charges(self, charges: list[str], invert: bool = False):
        """Filter by 'charge'.

        Args:
            charges (list[str]): a list of 'charge's.
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if invert:
            return self[~self.charge.str.strip().isin(charges)]
        return self[self.charge.str.strip().isin(charges)]

    def nmr_models(self, models: list[int], invert=False):
        """Filter by 'nmr_model'.

        Args:
            models (list[int]): a list of 'nmr_model's.
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        if "nmr_model" in self.columns:
            if invert:
                return self[~self.nmr_model.isin(models)]
            return self[self.nmr_model.isin(models)]
        elif models == [1]:
            if invert:
                return self[self.record_name == "IMPOSSIBLE"]
            return self

    def distances(
        self,
        other: tuple[float, float, float] | PdbDataFrame,
        cut_off: float = np.inf,
        to: str | None = None,
        invert=False,
    ):
        """Filter by distance to a reference point or group of atoms.

        Args:
            other (tuple[float, float, float] | list[float] | PdbDataFrame): the other group's coordinate(s).
            cut_off (float; defaults to numpy.inf): the distance cutoff to filter.
            to (str|None; defaults to None): if 'other' is a group atoms, using which method to determine
                whether the atoms meet the cut_off distance. If None, "COM" or center of mass is used.
            invert (bool; defaults to False): whether to invert the selection.

        Returns:
            PdbDataFrame: sub DataFrame after the filtering
        """
        cut_off = cut_off**2
        if to is None and isinstance(other, PdbDataFrame):
            to = "COM"
        other_center: tuple[float, float, float] | _lru_cache_wrapper[
            tuple[float, float, float]
        ] | None = None
        if to == "COM":
            if isinstance(other, PdbDataFrame):
                other_center = other.center_of_mass
            else:
                other_center = tuple(other)  # type: ignore
        elif to == "COG":
            if isinstance(other, PdbDataFrame):
                other_center = other.center_of_geometry
            else:
                other_center = tuple(other)  # type: ignore
        elif to in ["Any", "All"]:
            other_center = None
        else:
            raise NotImplementedError(
                f"Only center of mass (COM) and center of geometry (COG) is implemented, not {to}."
            )

        indices = []
        if to == "All":
            indices = list(self.atoms.index)
        for index, self_x, self_y, self_z in zip(
            self.atoms.index, self.atoms.x_coord, self.atoms.y_coord, self.atoms.z_coord
        ):
            if isinstance(other_center, tuple):
                distance = (
                    (self_x - other_center[0]) ** 2
                    + (self_y - other_center[1]) ** 2
                    + (self_z - other_center[2]) ** 2
                )  # type: ignore
                if distance <= cut_off and not invert:
                    indices.append(index)
                elif distance > cut_off and invert:
                    indices.append(index)
            elif isinstance(other, PdbDataFrame):
                for other_x, other_y, other_z in zip(
                    other.atoms.x_coord, other.atoms.y_coord, other.atoms.z_coord
                ):
                    distance = (
                        (self_x - other_x) ** 2
                        + (self_y - other_y) ** 2
                        + (self_z - other_z) ** 2
                    )
                    if to == "Any":
                        if (
                            distance <= cut_off
                            and (not invert)
                            and index not in indices
                        ):
                            indices.append(index)
                        elif distance > cut_off and invert and index not in indices:
                            indices.append(index)
                    else:  # 'All'
                        if distance > cut_off and (not invert) and index in indices:
                            indices.remove(index)
                        elif distance <= cut_off and invert and index in indices:
                            indices.remove(index)
        return self.iloc[indices]

    @staticmethod
    def get_residues(pdb_df: PdbDataFrame | pd.DataFrame) -> list[tuple]:
        """Function to get the list of residues.

        Args:
            pdb_df (PdbDataFrame|pd.DataFrame): a PdbDataFrame object.

        Returns:
            list[tuple]: a list of residues as (chain_id, residue_name, residue_number).
        """
        all_residues: list[tuple] = []
        for chain, residue_name, residue_number in zip(
            pdb_df["chain_id"], pdb_df["residue_name"], pdb_df["residue_number"]
        ):
            if residue_name.strip() not in AMINO_ACIDS:
                continue
            residue = (chain.strip(), residue_name.strip(), residue_number)
            if len(all_residues) == 0:
                all_residues.append(residue)
            elif residue != all_residues[-1]:
                all_residues.append(residue)
        return all_residues

    @functools.lru_cache()
    @staticmethod
    def get_distance_matrix(pdb_df: PdbDataFrame, use_r2: bool = True) -> np.ndarray:
        """Function to calculate the distance matrix.

        Args:
            pdb_df (PdbDataFrame): a PdbDataFrame object.
            use_r2 (bool; defaults to True): whether to use r2 or r for distance matrix.

        Returns:
            np.ndarray: distance matrix.
        """
        cols = ["x_coord", "y_coord", "z_coord"]
        if use_r2:
            return squareform(pdist(pdb_df[cols], "sqeuclidean"))
        else:
            return squareform(pdist(pdb_df[cols]))
