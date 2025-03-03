""" Recognizing sites to apply modifier on """

from typing import Optional, Callable, List, Union
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

    SCIPY_INST = True
except ImportError:
    SCIPY_INST = False


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
    """A subclass of Sites that uses multiple rules to identify sites in an atomic structure."""

    def __init__(
        self,
        index_parsers: Optional[
            Union[Callable[[Atoms], list], List[Callable[[Atoms], list]]]
        ] = None,
        combine_rules: str = "union",
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
    ):
        """
        Args:
            index_parsers (Union[Callable[[Atoms], list], List[Callable[[Atoms], list]]], optional):
            A single rule or a list of rules (functions) that take an Atoms object
            of indices representing the sites of interest. Defaults to a function that returns all indices.
            combine_rules (str, optional):How to combine the results of multiple rules. Options are:
                - "union": Combine results using set union (default).
                - "intersection": Combine results using set intersection.
                Defaults to "union".
            max_bond_ratio (float, optional): While making bonds, how much error is allowed.
                Defaults to 1.2.
            max_bond (float, optional): Fixed bond distance to use, any distance above is ignored.
                Defaults to 0. If 0, it is ignored.
            contact_error (float, optional): Error allowed if atoms are too close to each other.
                Defaults to 0.2.
        """
        super().__init__(max_bond_ratio, max_bond, contact_error)

        # Default index_parser function that returns all indices
        if index_parsers is None:
            self.index_parsers = [lambda atoms: list(range(len(atoms)))]
        else:
            # Ensure index_parsers is always a list, even if a single function is provided
            self.index_parsers = (
                [index_parsers] if callable(index_parsers) else index_parsers
            )

        if combine_rules not in ["union", "intersection"]:
            raise ValueError("combine_rules must be 'union' or 'intersection'.")
        self.combine_rules = combine_rules

    def get_sites(self, atoms: Atoms) -> list:
        """
        Args:
            atoms (Atoms): The atomic structure to analyze.

        Returns:
            list: A list of indices representing the sites of interest.
        """
        # Apply each rule to the atoms object
        results = [parser(atoms) for parser in self.index_parsers]

        # Combine results based on the specified logic
        if self.combine_rules == "union":
            # Use set union to combine results
            return list(set().union(*results))
        elif self.combine_rules == "intersection":
            # Use set intersection to combine results
            return list(set(results[0]).intersection(*results[1:]))


# Rules Defined
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


def get_com_sites(
    atoms: Atoms, fraction: float = 1.0, direction: str = "above"
) -> list:
    """

    Args:
        atoms (Atoms): The ASE Atoms object.
        fraction (float, optional): Fraction of the z-distance from the COM to consider.
                                    Must be between 0 and 1. Defaults to 1.0.
        direction (str, optional): Whether to return atoms "above", "below", or "both".
                                   Defaults to "above".

    Returns:
        list: A list of atom indices based on the specified direction.
    """
    if not 0 <= fraction <= 1:
        raise ValueError("Fraction must be between 0 and 1.")
    if direction not in ["above", "below", "both"]:
        raise ValueError("Direction must be 'above', 'below', or 'both'.")

    if len(atoms) == 0:
        return []

    # Get the center of mass (COM) z-coordinate
    com_z = atoms.get_center_of_mass()[2]

    # Get the z-coordinates of all atoms
    z_values = atoms.get_positions()[:, 2]

    # Compute max and min z-values
    z_max = max(z_values)
    z_min = min(z_values)

    # Calculate z-thresholds for above and below
    threshold_z_above = com_z + (1 - fraction) * (z_max - com_z)
    threshold_z_below = com_z - (1 - fraction) * (com_z - z_min)

    above_indices = {i for i, z in enumerate(z_values) if z > threshold_z_above}
    below_indices = {i for i, z in enumerate(z_values) if z < threshold_z_below}

    if direction == "above":
        return list(above_indices)
    elif direction == "below":
        return list(below_indices)
    else:
        return list(above_indices | below_indices)


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
        if diff_coord > 0:
            sites.append(
                {
                    "ind": index,
                    "coord": coord,
                    "diff_coord": diff_coord,
                    "z_coord": atoms[index].position[2],
                }
            )

    # Convert to DataFrame for easier filtering
    df = DataFrame(sites)

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

    if not SCIPY_INST:
        print("Scipy isnt installed; get_surface_sites_by_voronoi_pbc wont work")
        return list(range(len(atoms)))

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
    for n_x, n_y, n_z in offset_combinations:
        shift = n_x * cell[0] + n_y * cell[1] + n_z * cell[2]
        for i, pos in enumerate(positions):
            extended_positions.append(pos + shift)
            tags.append((n_x, n_y, n_z))
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


def get_surface_by_normals(
    atoms: Atoms,
    surface_normal: float = 0.5,
    tolerance: float = 1e-5,
    normalize_final: bool = True,
    self_interaction: bool = False,
    bothways: bool = True,
) -> List:
    """
    Compute surface normals for an ASE Atoms object, considering periodic boundaries.

    Args:
        atoms (Atoms): ASE Atoms object.
        surface_normal (float): Threshold for identifying surface atoms based on normal magnitude.
        normalize_final (bool): Whether to normalize output normals.
        adsorbate_atoms (list): Indices of adsorbate atoms to exclude.

    Returns:
        np.ndarray: Surface normals for each atom.
        list: Indices of detected surface atoms.
    """
    atoms = atoms.copy()

    # Create the graph
    nl = NeighborList(
        natural_cutoffs(atoms), self_interaction=self_interaction, bothways=bothways
    )
    nl.update(atoms)

    normals = np.zeros((len(atoms), 3), dtype=float)
    for index in range(len(atoms)):
        normal = np.zeros(3, dtype=float)
        atom_pos = atoms.positions[index]

        for neighbor, offset in zip(*nl.get_neighbors(index)):
            neighbor_pos = atoms.positions[neighbor] + np.dot(offset, atoms.cell)
            normal += atom_pos - neighbor_pos  # Vector sum of neighbor directions

        # Store normal vector if it's above threshold
        if np.linalg.norm(normal) > surface_normal:
            normals[index, :] = (
                normal / np.linalg.norm(normal) if normalize_final else normal
            )

    # Identify surface atoms based on normal magnitude
    surface = [
        index
        for index in range(len(atoms))
        if np.linalg.norm(normals[index]) > tolerance
    ]
    return surface
