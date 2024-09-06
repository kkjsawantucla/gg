""" Recognizing sites to apply modifier on """

from typing import Optional
import pandas as pd
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
from gg.utils_graph import atoms_to_graph


class Sites:
    """Base class for sites"""

    def __init__(self):
        self.graph = None

    @property
    def graph(self) -> nx.Graph:
        """
        Returns:
            nx.Graph: _description_
        """
        return self.g

    @graph.setter
    def graph(self, g):
        self.g = g

    def get_sites(self, atoms: Atoms):
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class SurfaceSites(Sites):
    """A class that figures out the surface atoms"""

    def __init__(
        self,
        max_coord: dict,
        max_bond_ratio: Optional[float] = 0,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
        com: Optional[bool] = True,
    ):
        super().__init__()
        self.max_coord = max_coord
        self.max_bond_ratio = max_bond_ratio
        self.max_bond = max_bond
        self.contact_error = contact_error
        self.com = com

    def get_graph(
        self, atoms: Atoms, self_interaction: bool = False, bothways: bool = True
    ):
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

    def get_sites(
        self, atoms: Atoms, self_interaction: bool = False, bothways: bool = True
    ):
        """
        Args:
            atoms (_type_): _description_
            self_interaction (bool, optional): _description_. Defaults to False.
            both ways (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        for sym in atoms.symbols:
            if sym not in list(self.max_coord.keys()):
                raise RuntimeError(f"Incomplete max_coord: Missing {sym}")
        if not self.graph:
            self.get_graph(atoms, self_interaction=self_interaction, bothways=bothways)
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

        df = pd.DataFrame(sites)

        if self.com:
            df = df[df.z_coord > atoms.get_center_of_mass()[2]]

        df = df[df.diff_cord > 0]
        df = df.sort_values(by=["cord", "z_coord"])
        return df["ind"].to_list()
