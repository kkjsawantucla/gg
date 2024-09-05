"""Importing modules for dealing with graph utilities"""

from typing import Optional, Tuple
from itertools import combinations
import numpy as np
import networkx as nx
from numpy.linalg import norm
import pandas as pd
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.io import read as read_atoms
from ase.data import covalent_radii
from ase.build import molecule
from gg.utils_graph import node_symbol, relative_position, atoms_to_graph
from gg.utils_graph import is_cycle, are_points_collinear_with_tolerance


__author__ = "Kaustubh Sawant"


class NoReasonableStructureFound(Exception):
    """ Custom error to handle touching atoms """


class SurfaceSites:
    """_summary_"""

    def __init__(
        self,
        max_coord: dict,
        surf_atom_sym: Optional[list] = None,
        max_bond_ratio: Optional[float] = 0,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
        com: Optional[bool] = True,
    ):
        self.max_coord = max_coord
        self.surf_atom_sym = surf_atom_sym
        self.max_bond_ratio = max_bond_ratio
        self.max_bond = max_bond
        self.contact_error = contact_error
        self.com = com

    def get_surface_sites(self, atoms, self_interaction=False, bothways=True):
        """
        Args:
            atoms (_type_): _description_
            self_interaction (bool, optional): _description_. Defaults to False.
            both ways (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        if not self.surf_atom_sym:
            self.surf_atom_sym = list(set(atoms.symbols))
        for sym in atoms.symbols:
            if sym not in list(self.max_coord.keys()):
                raise RuntimeError(f"Incomplete max_coord: Missing {sym}")

        nl = NeighborList(
            natural_cutoffs(atoms), self_interaction=self_interaction, bothways=bothways
        )
        nl.update(atoms)
        g = atoms_to_graph(
            atoms, nl, max_bond_ratio=self.max_bond_ratio, max_bond=self.max_bond
        )
        sites = []
        for node in g.nodes():
            cord = len([edge for edge in g[node]])
            index = g.nodes[node]["index"]
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

        df = df[df.diff_cord > 0].sort_values(by=["symbol", "cord"])
        if isinstance(self.surf_atom_sym, str):
            df = df[df["symbol"] == self.surf_atom_sym]
        else:
            df = df[df["symbol"].isin(self.surf_atom_sym)]
        df = df.sort_values(by=["cord", "z_coord"])
        return df["ind"].to_list(), g


def custom_copy(atoms: Atoms) -> Atoms:
    """copy atoms with calc
    Args:
        atoms (ase.Atoms):

    Returns:
        ase.Atoms:
    """
    atoms_copy = atoms.copy()
    if atoms.get_calculator():
        # calc_type = type(atoms.calc)
        # calc_params = atoms.calc.parameters
        # atoms_copy.calc = calc_type(**calc_params)
        atoms_copy.calc = atoms.calc
    return atoms_copy


# Function to ad adsorbate to atoms object
def add_ads(atoms: Atoms, ads: Atoms, offset: float) -> Atoms:
    """ Add adsorbate on a substrate with a particular offset
    Args:
        atoms (Atoms Object): Substrate
        ads (Atoms Object): Adsorbate to add
        offset (_type_): ref distance in the substrate where adsorbate will be added

    Returns:
        atoms (Atoms Object): _description_
    """

    if isinstance(ads, str):
        _ads = read_atoms(ads)
    elif isinstance(ads, Atoms):
        _ads = ads.copy()
    else:
        print("Please provide proper adsorbate file")
    for atom in _ads:
        atom.position[0] = atom.position[0] + offset[0]
        atom.position[1] = atom.position[1] + offset[1]
        atom.position[2] = atom.position[2] + offset[2]
    sub = custom_copy(atoms)
    sub += _ads
    return sub


def get_normals(index: list, atoms: Atoms, g: nx.Graph) -> Tuple[np.array, np.array]:
    """
    Args:
        index (_type_):
        atoms (ase.Atoms):
        G (nx.Graph):

    Returns:
        normal vector, reference position:
    """
    # Initially, get the right position vectors for the adsorption cluster
    ads_pos = np.zeros((len(index), 3))
    initial = index[0]
    ads_pos[0] = np.array([0, 0, 0])
    for i, j in enumerate(index[1:]):
        atom = atoms[j]
        edge_data = g.get_edge_data(node_symbol(atoms[initial]), node_symbol(atom))
        vector = edge_data["weight"]
        start = edge_data["start"]
        if start == initial:
            ads_pos[i + 1] = vector
        else:
            ads_pos[i + 1] = -vector

    ads_pos_sum = np.sum(-ads_pos, axis=0) / len(index)

    # Construct a Matrix with vectors surrounding the adsorption cluster
    normals = []
    for i, j in enumerate(index):
        atom = atoms[j]
        for neighbor in g.neighbors(node_symbol(atom)):
            n_index = g.nodes[neighbor]["index"]
            if n_index not in index:
                edge_data = g.get_edge_data(node_symbol(atom), neighbor)
                vector = edge_data["weight"]
                start = edge_data["start"]
                if start == j:
                    normal = vector
                else:
                    normal = -vector
                normal = normal + ads_pos_sum + ads_pos[i]
                normals.append(normal / norm(normal))
    v = np.array(normals)

    # Find the direction that best represents the empty space around the adsorption cluster
    _, _, vt = np.linalg.svd(v)
    vt = vt[-1]
    # The best normal could be +Vt or -Vt
    matrix1 = np.array([vt, -vt])
    dot_products = np.dot(matrix1, v.T)
    sums_per_vector = np.sum(dot_products, axis=1)
    max_index = np.argmax(sums_per_vector)
    vector_with_smallest_sum = matrix1[max_index]

    ref_pos = ads_pos_sum + atoms[initial].position
    return vector_with_smallest_sum, ref_pos


def move_along_normal(index: int, atoms: Atoms, g: nx.Graph) -> Atoms:
    """ Helper function to move atoms along the free normal """
    normal, ref_pos = get_normals([index], atoms, g)
    atom = atoms[index]
    offset = 0
    i = 0
    for neighbor in g.neighbors(node_symbol(atom)):
        atom_n = atoms[g.nodes[neighbor]["index"]]
        offset += (
            covalent_radii[atom_n.number]
            + covalent_radii[atom.number]
            - atoms.get_distance(index, atom_n.index, mic=True)
        )
        i += 1
    offset = offset / i
    atom.position = ref_pos + normal * offset / norm(normal)
    return atoms


def generate_sites(
    atoms: Atoms,
    ads: Atoms,
    graph: nx.Graph,
    index: list,
    coordination: int,
    ad_dist: float = 1.7,
    contact_error: float = 0.2,
):
    """
    Args:
        atoms (ase.Atoms):
        ads (ase.Atoms):
        graph (nx.Graph):
        index (list):
        coordination (int):
        ad_dist (float):  Defaults to 1.7.
        contact_error (float): _ Defaults to 0.2.

    Returns:
       ase.Atoms:
    """
    possible = list(combinations(index, coordination))
    valid = []

    for cycle in possible:
        if coordination == 1:
            valid.append(list(cycle))

        if coordination == 2:
            if graph.has_edge(
                node_symbol(atoms[cycle[0]]), node_symbol(atoms[cycle[1]])
            ):
                valid.append(list(cycle))

        if coordination == 3:
            nodes = [node_symbol(atoms[i]) for i in cycle]
            pos = [atoms[i].position for i in cycle]
            if is_cycle(graph, nodes):
                if not are_points_collinear_with_tolerance(
                    pos[0], pos[1], pos[2], tolerance=0.01
                ):
                    valid.append(list(cycle))

    movie = []

    for cycle in valid:
        normal, ref_pos = get_normals(cycle, atoms, graph)
        offset = normal * ad_dist / norm(normal)
        ads_copy = ads.copy()
        ads_copy.rotate([0, 0, 1], normal, center=[0, 0, 0])
        atoms_copy = add_ads(atoms, ads_copy, offset=offset + ref_pos)
        if check_contact(atoms_copy, error=contact_error):
            continue
        else:
            movie.append(atoms_copy)
    return movie


def formula_to_graph(formula, max_bond_ratio=1.2, max_bond=0):
    """
    Args:
        formula (str) or (ase.Atoms)

    Returns:
        _type_: _description_
    """
    if isinstance(formula, str):
        atoms = molecule(formula)
    elif isinstance(formula, Atoms):
        atoms = read_atoms(formula)
    else:
        raise RuntimeError("Issue in reading formula")

    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    g = atoms_to_graph(atoms, nl, max_bond_ratio=max_bond_ratio, max_bond=max_bond)

    return g


def check_contact(atoms, error=0.1, print_contact=False):
    """ Check if atoms touch within an error tolerance
    Args:
        atoms (ase.Atoms):
        error (float, optional): Error tolerated. Defaults to 0.1.
        print_contact (bool, optional): Print contact information. Defaults to False.

    Returns:
        Boolean:
    """
    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    close_contact = []
    for index, atom in enumerate(atoms):
        for neighbor, offset in zip(*nl.get_neighbors(index)):
            if sorted((index, neighbor)) not in close_contact:
                atom2 = atoms[neighbor]
                distance = np.linalg.norm(
                    atom.position - relative_position(atoms, neighbor, offset)
                )
                eqm_radii = covalent_radii[atom.number] + covalent_radii[atom2.number]
                if distance < (1 - error) * eqm_radii:
                    if print_contact:
                        print(
                            "Close Contact:",
                            node_symbol(atom),
                            node_symbol(atom2),
                            round(distance, 2),
                        )
                    close_contact.append(sorted((index, neighbor)))
    if close_contact:
        return True
    else:
        return False
