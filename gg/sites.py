""" Recognizing sites to apply modifier on """

from typing import Optional
from pandas import DataFrame
from ase import Atoms
from ase.constraints import FixAtoms
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
from gg.utils_graph import atoms_to_graph


class Sites:
    """Base class for sites"""

    def __init__(
        self,
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
    ):
        """
        Args: All the variables help in making graphs

            max_bond_ratio (float, optional): While making bonds how much error is allowed.
            Defaults to 1.2.

            max_bond (float, optional): Fixed bond distance to use, any distance above is ignored.
            Defaults to 0. If 0 , it is ignored

            contact_error (float, optional): Error allowed if atoms are too close to each other.
            Defaults to 0.2.

        """
        self.graph = None
        self.max_bond_ratio = max_bond_ratio
        self.max_bond = max_bond
        self.contact_error = contact_error

    @property
    def graph(self) -> nx.Graph:
        """
        Returns:
            nx.Graph:
        """
        return self.g

    @graph.setter
    def graph(self, g):
        self.g = g

    def get_graph(
        self, atoms: Atoms, self_interaction: bool = False, bothways: bool = True
    ) -> nx.Graph:
        """
        Args:
            atoms (_type_): _description_
            self_interaction (bool, optional): _description_. Defaults to False.
            both ways (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """

        nl = NeighborList(
            natural_cutoffs(atoms), self_interaction=self_interaction, bothways=bothways
        )
        nl.update(atoms)
        g = atoms_to_graph(
            atoms, nl, max_bond_ratio=self.max_bond_ratio, max_bond=self.max_bond
        )
        self.graph = g
        return g

    def get_sites(self, atoms: Atoms) -> list:
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class FlexibleSites(Sites):
    """A class which returns sites identified by index or constariants"""

    def __init__(
        self,
        constraints: bool = False,
        index: list = None,
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
    ):
        """_summary_

        Args:
            constraints (bool, optional): If true, only atoms which arent constrained considered.
            Defaults to False.

            index (list, optional): If list if indices is give, it will be used as it is.
            Defaults to None.

        Raises:
            RuntimeError: _description_
        """
        super().__init__(max_bond_ratio, max_bond, contact_error)

        if index is not None or constraints is True:
            self.index = index
            self.constraints = constraints
        else:
            raise RuntimeError("Specify either index or constraits")

    def get_sites(self, atoms: Atoms) -> list:
        if self.index:
            index = self.index

        else:
            index = range(len(atoms))

        if self.constraints:
            constrained_indices = set()
            for constraint in atoms.constraints:
                if isinstance(constraint, FixAtoms):
                    constrained_indices.update(constraint.index)

            unconstrained_indices = [i for i in index if i not in constrained_indices]
            return unconstrained_indices
        else:
            return list(index)


class SurfaceSites(Sites):
    """A class that figures out the surface atom by coordination number"""

    def __init__(
        self,
        max_coord: dict,
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
        com: Optional[bool] = True,
    ):
        """_summary_

        Args:
            max_coord (dict): Dictionary of the maximum coordination of each element used.
            Only atoms with coordination less than this value will be considered.

            com (Optional[bool], optional): If true atoms below the center of mass are ignored
            Defaults to True.
        """
        super().__init__(max_bond_ratio, max_bond, contact_error)
        self.max_coord = max_coord
        self.com = com
        self.df = None

    def get_sites(
        self, atoms: Atoms, self_interaction: bool = False, bothways: bool = True
    ) -> list:
        """
        Args:
            atoms (ase.Atoms):
            self_interaction (bool, optional): Defaults to False.
            both ways (bool, optional): Defaults to True.

        Returns:
            _type_: _description_
        """
        for sym in atoms.symbols:
            if sym not in list(self.max_coord.keys()):
                raise RuntimeError(f"Incomplete max_coord: Missing {sym}")

        _ = self.get_graph(atoms, self_interaction=self_interaction, bothways=bothways)
        sites = []
        for node in self.graph.nodes():
            cord = len([edge for edge in self.graph[node]])
            index = self.graph.nodes[node]["index"]
            symbol = atoms[index].symbol
            diff_cord = self.max_coord[symbol] - cord
            sites.append(
                {
                    "ind": index,
                    "symbol": symbol,
                    "cord": cord,
                    "diff_cord": diff_cord,
                    "z_coord": atoms[index].position[2],
                }
            )

        df = DataFrame(sites)

        if self.com:
            df = df[df.z_coord > atoms.get_center_of_mass()[2]]

        df = df[df.diff_cord > 0]
        df = df.sort_values(by=["cord", "z_coord"])
        self.df = df
        return df["ind"].to_list()
