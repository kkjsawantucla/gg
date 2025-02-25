""" Recognizing sites to apply modifier on """

from typing import Optional, Callable, List
from itertools import product
from pandas import DataFrame
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.neighborlist import NeighborList, natural_cutoffs
import networkx as nx
from gg.utils_graph import atoms_to_graph

try:
    from scipy.spatial import Voronoi

    scipy_installed = True
except ImportError:
    scipy_installed = False


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
        return self.graph

    def get_sites(self, atoms: Atoms) -> list:
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class RuleSites(Sites):
    """_summary_"""

    def __init__(
        self,
        index_parser: Optional[Callable[[Atoms], list]] = None,
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
    ):
        super().__init__(max_bond_ratio, max_bond, contact_error)
        if index_parser is None:
            raise RuntimeError("A valid index_parser function must be provided.")
        self.index_parser = index_parser

    def get_sites(self, atoms: Atoms) -> list:
        return self.index_parser(atoms)


# Rule Defined
def get_unconstrained_sites(atoms: Atoms) -> list:
    """
    Returns a list of indices of atoms that are not constrained.
    """
    constrained_indices = set()
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            constrained_indices.update(constraint.index)
    return [i for i in range(len(atoms)) if i not in constrained_indices]


def get_tagged_sites(atoms: Atoms, tag: int = -1) -> list:
    """
    Returns a list of indices of atoms that have the specified tag.
    """
    return [i for i in range(len(atoms)) if atoms[i].tag == tag]


def get_above_com_sites(atoms: Atoms, perc: float = 1.0) -> list:
    """
    Returns a list of indices of atoms that are above the center of mass in the z-direction.
    The `percentage` param determines the fraction of the z-dist above the COM that is considered.

    Args:
        atoms (Atoms): The ASE Atoms object.
        percentage (float, optional): Fraction of the z-distance above the COM to consider.
                    Must be between 0 and 1. Defaults to 1.0 (all sites above COM).

    Returns:
        list: List of atom indices above the specified z-threshold.
    """
    if not 0 <= perc <= 1:
        raise ValueError("Percentage must be between 0 and 1.")

    # Get the center of mass (COM) z-coordinate
    com_z = atoms.get_center_of_mass()[2]

    # Get the z-coordinates of all atoms
    z_values = atoms.get_positions()[:, 2]

    # Calculate the maximum z-coordinate
    z_max = np.max(z_values)

    # Determine the threshold z-coordinate
    threshold_z = com_z + (1-perc) * (z_max - com_z)

    # Select atoms above the threshold z-coordinate
    above_threshold_indices = [i for i, z in enumerate(z_values) if z > threshold_z]

    return above_threshold_indices


def get_surface_sites_by_coordination(
    atoms: Atoms,
    max_coord: dict,
    max_bond_ratio: float = 1.2,
    max_bond: float = 0,
    self_interaction: bool = False,
    bothways: bool = True,
) -> list:
    """
    Identifies surface sites based on coordination numbers using a graph-based approach.

    Args:
        atoms (Atoms): The ASE Atoms object.
        max_coord (Dict[str, int]): Dictionary of maximum coordination numbers for each element.
        max_bond_ratio (float, optional): Tolerance for bond distances. Defaults to 1.2.
        max_bond (float, optional): Maximum bond distance. Defaults to 0 (no limit).
        contact_error (float, optional): Tolerance for atoms being too close. Defaults to 0.2.
        com (float, optional): Fraction of the z-range above the center of mass to consider.
                              Defaults to 0.1.
        self_interaction (bool, optional): Whether to include self-interactions in the graph.
                                           Defaults to False.
        bothways (bool, optional): Whether to consider bonds in both directions. Defaults to True.

    Returns:
        list: List of atom indices identified as surface sites.
    """
    # Validate max_coord
    for sym in atoms.symbols:
        if sym not in max_coord:
            raise RuntimeError(f"Incomplete max_coord: Missing {sym}")

    # Create the graph
    nl = NeighborList(
        natural_cutoffs(atoms), self_interaction=self_interaction, bothways=bothways
    )
    nl.update(atoms)
    graph = atoms_to_graph(atoms, nl, max_bond_ratio=max_bond_ratio, max_bond=max_bond)

    # Calculate coordination numbers and filter surface sites
    sites = []
    for node in graph.nodes():
        coord = len(list(graph[node]))  # Coordination number
        index = graph.nodes[node]["index"]
        symbol = atoms[index].symbol
        diff_coord = max_coord[symbol] - coord
        sites.append(
            {
                "ind": index,
                "symbol": symbol,
                "coord": coord,
                "diff_coord": diff_coord,
                "z_coord": atoms[index].position[2],
            }
        )

    # Convert to DataFrame for easier filtering
    df = DataFrame(sites)

    # Filter based on coordination number
    df = df[df.diff_coord > 0]

    # Sort by coordination number and z-coordinate
    df = df.sort_values(by=["coord", "z_coord"])

    # Return the list of indices
    return df["ind"].to_list()


def get_surface_sites_by_voronoi_pbc(atoms: Atoms) -> List[int]:
    """
    Identifies surface atoms using Voronoi tessellation with periodic boundary conditions.

    Args:
        atoms (Atoms): ASE Atoms object with PBC settings.

    Returns:
        List[int]: Indices of surface atoms.
    """
    cell = atoms.cell
    pbc = atoms.pbc
    positions = atoms.get_positions()
    extended_positions = []
    tags = []
    original_indices = []

    # Generate offsets based on PBC settings (e.g., [-1, 0, 1] for periodic dimensions)
    offsets = []
    for dim in range(3):
        if pbc[dim]:
            offsets.append([-1, 0, 1])
        else:
            offsets.append([0])

    # Create all offset combinations (e.g., 3x3x1 for a 2D-periodic slab)
    offset_combinations = product(*offsets)

    # Replicate atoms in neighboring cells based on PBC
    for n_x, ny, nz in offset_combinations:
        shift = n_x * cell[0] + ny * cell[1] + nz * cell[2]
        for i, pos in enumerate(positions):
            extended_positions.append(pos + shift)
            tags.append((n_x, ny, nz))
            original_indices.append(i)  # Track original atom index

    # Compute Voronoi tessellation
    vor = Voronoi(extended_positions)

    surface_indices = set()

    # Check for unbounded Voronoi cells in non-periodic directions
    for i, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        # If the region is unbounded (-1 in region), it's a surface atom
        if -1 in region:
            # Check if the atom is in the central cell
            if tags[i] == (0, 0, 0):
                surface_indices.add(original_indices[i])

    return sorted(surface_indices)
