""" General utilities for the modifiers """

from typing import Tuple, Union
import numpy as np
import networkx as nx
from numpy.linalg import norm
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.data import covalent_radii, chemical_symbols
from ase.build import molecule
from ase.collections import g2
from gg.utils_graph import node_symbol, relative_position, atoms_to_graph
from gg.data import adsorbates
from pandas import read_fwf


__author__ = "Kaustubh Sawant"


class NoReasonableStructureFound(Exception):
    """Custom error to handle unrealistic or empty atoms"""


def custom_copy(atoms: Atoms) -> Atoms:
    """copy atoms with calc
    Args:
        atoms (ase.Atoms):

    Returns:
        ase.Atoms
    """
    atoms_copy = atoms.copy()
    if atoms.calc:
        atoms_copy.calc = atoms.calc
    return atoms_copy


# Function to ad adsorbate to atoms object
def add_ads(atoms: Atoms, ads: Atoms, offset: float, tag: bool = False) -> Atoms:
    """Add adsorbate on a substrate with a particular offset
    Args:
        atoms (Atoms Object): Substrate
        ads (Atoms Object): Adsorbate to add
        offset (_type_): ref distance in the substrate where adsorbate will be added

    Returns:
        atoms (Atoms Object): _description_
    """
    if isinstance(ads, str):
        if ads in adsorbates:
            _ads = adsorbates[ads]
        elif ads in g2.names:
            _ads = molecule(ads)
        elif ads in chemical_symbols:
            _ads = Atoms(ads, positions=[(0, 0, 0)])
        else:
            raise RuntimeError(f"Cannot convert string to Formula {ads}")
    elif isinstance(ads, Atoms):
        _ads = ads.copy()
    else:
        raise RuntimeError(f"Please provide proper file {ads}")

    if tag:
        for atom in _ads:
            atom.tag = -1
    for atom in _ads:
        atom.position[0] = atom.position[0] + offset[0]
        atom.position[1] = atom.position[1] + offset[1]
        atom.position[2] = atom.position[2] + offset[2]
    sub = custom_copy(atoms)
    sub += _ads
    return sub


def get_normals(
    index: list, atoms: Atoms, g: nx.Graph, method="svd"
) -> Tuple[np.array, np.array]:
    """
    Args:
        index (_type_):
        atoms (ase.Atoms):
        G (nx.Graph):
        method:
    Returns:
        normal vector, reference position:
    """
    # Initially, get the right position vectors for the adsorption cluster
    # This is important if the adsorption site is between 2+ atoms
    ads_pos = np.zeros((len(index), 3))
    initial = index[0]
    ads_pos[0] = np.array([0, 0, 0])
    for i, j in enumerate(index[1:]):
        atom = atoms[j]
        if g.has_edge(node_symbol(atoms[initial]), node_symbol(atom)):
            edge_data = g.get_edge_data(node_symbol(atoms[initial]), node_symbol(atom))
            vector = edge_data["weight2"]
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
                vector = edge_data["weight2"]
                start = edge_data["start"]
                if start == j:
                    normal = vector
                else:
                    normal = -vector
                normal = normal + ads_pos_sum + ads_pos[i]
                normals.append(normal)

    v = np.array(normals)
    # Find the direction that best represents the empty space around the adsorption cluster
    try:
        _, s, vt = np.linalg.svd(v, full_matrices=False)
    except np.linalg.LinAlgError:
        return np.array([0, 0, 1]), ads_pos_sum + atoms[initial].position
    else:
        rank = np.sum(s > 1e-8)

        if rank < 3 or method == "svd":
            vt = vt[-1]
        elif method == "mean":
            vt = np.mean(v, axis=0)

        # The best normal could be +Vt or -Vt
        matrix1 = np.array([vt, -vt])
        dot_products = np.dot(matrix1, v.T)
        sums_per_vector = np.sum(dot_products, axis=1)
        max_index = np.argmax(sums_per_vector)
        vector_with_smallest_sum = matrix1[max_index]
        # Return the vector direction and position from where to start
        ref_pos = ads_pos_sum + atoms[initial].position
        return vector_with_smallest_sum, ref_pos


def move_along_normal(index: int, atoms: Atoms, g: nx.Graph) -> Atoms:
    """Helper function to move atoms along the free normal"""
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


def formula_to_graph(formula, max_bond_ratio=1.2, max_bond=0) -> nx.Graph:
    """
    Args:
        formula (str) or (ase.Atoms)

    Returns:
        _type_: _description_
    """
    if isinstance(formula, str):
        if formula in adsorbates:
            atoms = adsorbates[formula]
        elif formula in g2.names:
            atoms = molecule(formula)
        elif formula in chemical_symbols:
            atoms = Atoms(formula, positions=[(0, 0, 0)])
        else:
            raise RuntimeError(f"Cannot convert string to Formula {formula}")
    elif isinstance(formula, Atoms):
        atoms = formula
    else:
        raise RuntimeError("Issue in reading formula")

    # Make the ase.NeighborList
    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)

    # Make graph from atoms
    g = atoms_to_graph(atoms, nl, max_bond_ratio=max_bond_ratio, max_bond=max_bond)

    return g


def check_contact(atoms, error=0.1, print_contact=False) -> bool:
    """Check if atoms touch within an error tolerance
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


def replace(atoms: Atoms, replace_with: Union[str, Atoms], offset: np.array) -> Atoms:
    """
    Args:
        atoms (Atoms):
        replace_with (str or atoms):
        offset (array):

    Returns:
        Atoms:
    """
    if isinstance(replace_with, str):
        if replace_with in adsorbates:
            rep_atoms = adsorbates[replace_with]
        elif replace_with in g2.names:
            rep_atoms = molecule(replace_with)
        elif replace_with in chemical_symbols:
            rep_atoms = Atoms(replace_with, positions=[(0, 0, 0)])
        else:
            raise RuntimeError(f"Cannot convert string to Formula {replace_with}")
        atoms = add_ads(atoms, rep_atoms, offset)

    elif isinstance(replace_with, Atoms):
        rep_atoms = replace_with
        rep_atoms_copy = rep_atoms.copy()
        atoms = add_ads(atoms, rep_atoms_copy, offset)
    return atoms


def is_within_tolerance(value: float, target: float, tolerance: float) -> bool:
    """
    Check if a value is within a specified tolerance of a target value.

    Parameters:
    value (float): The number to check.
    target (float): The target number.
    tolerance (float): The tolerance range.

    Returns:
    bool: True if value is within tolerance of target, False otherwise.
    """
    return abs(value - target) <= tolerance


def get_ref_pos_index(index: list, atoms: Atoms, g: nx.Graph) -> np.array:
    """
    Args:
        index (list):
        atoms (Atoms):
        g (nx.Graph):

    Returns:
        np.array:
    """
    # Initially, get the right position vectors for the adsorption cluster
    # This is important if the adsorption site is between 2+ atoms
    ads_pos = np.zeros((len(index), 3))
    initial = index[0]
    ads_pos[0] = np.array([0, 0, 0])
    for i, j in enumerate(index[1:]):
        atom = atoms[j]
        edge_data = g.get_edge_data(node_symbol(atoms[initial]), node_symbol(atom))
        vector = edge_data["weight2"]
        start = edge_data["start"]
        if start == initial:
            ads_pos[i + 1] = vector
        else:
            ads_pos[i + 1] = -vector

    ads_pos_sum = np.sum(-ads_pos, axis=0) / len(index)
    return ads_pos_sum + atoms[initial].position


def get_area(atoms: Atoms) -> float:
    """Return xy area of ase.Atoms"""
    a = np.array(
        [[atoms.cell[0][0], atoms.cell[0][1]], [atoms.cell[1][0], atoms.cell[1][1]]]
    )
    area = abs(np.linalg.det(a))
    return area


def extract_lowest_energy_from_oszicar(file_path):
    """Extract the lowest energy value from an OSZICAR file."""
    lowest_energy = float("inf")

    with open(file_path, "r") as file:
        for line in file:
            if "E0=" in line:
                parts = line.split("E0=")
                try:
                    energy = float(parts[1].split()[0])
                    lowest_energy = min(lowest_energy, energy)
                except (IndexError, ValueError):
                    print(
                        f"Skipping invalid energy format in {file_path}: {line.strip()}"
                    )

    return lowest_energy if lowest_energy != float("inf") else None


def extract_lowest_energy_from_outlog(file_path):
    """Extract the lowest energy value from an OSZICAR file."""

    df = read_fwf(file_path)
    return float(df["Energy"].iloc[-1])


def sort_atoms(atoms):
    """Return atoms in sorted symbol list"""
    #makes initial atom list
    symbols = []
    for atom in atoms:
        if atom.symbol not in symbols:
            symbols.append(atom.symbol)

    #Using the list above sorts the ase atom object
    sort = []
    for symbol in symbols:
        for m, atom in enumerate(atoms):
            if atom.symbol == symbol:
                sort.append(m)

    resort = list(range(len(sort)))
    for n in range(len(resort)):
        resort[sort[n]] = n

    atoms_sorted = atoms[sort]
    return atoms_sorted
