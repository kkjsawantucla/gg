"""Importing modules for dealing with graph utilities"""

import networkx as nx
import numpy as np
from numpy.linalg import norm
from ase.data import covalent_radii


__author__ = "Kaustubh Sawant"


def node_symbol(atom):
    """
    Args:
       atom (Atoms Object): atoms to convert

    Returns:
       str: "symbol_index"
    """
    return f"{atom.symbol}_{atom.index}"


def relative_position(atoms, neighbor, offset):
    """
    Args:
         atoms (ase.Atoms):
         neighbor (int): Index of the neighbor
         offset (array):

     Returns:
       np.array: position of neighbor wrt to offset
    """
    return atoms[neighbor].position + np.dot(offset, atoms.get_cell())


def node_match(n1, n2):
    """_summary_
    Args:
        n1 (str):
        n2 (str):
    Returns:
        Boolean: _description_
    """
    return n1["symbol"] == n2["symbol"]


def is_cycle(g, nodes):
    """Check if the nodes in graph G form a cycle
    Args:
       G (networkx Graph):
       nodes ([list of networkx nodes]):

    Returns:
       Boolean: True if they form cycle
    """
    start_node = next(iter(nodes))  # Get any node as starting point
    subgraph = g.subgraph(nodes)
    try:
        nx.find_cycle(subgraph, source=start_node)
        return True
    except nx.NetworkXNoCycle:
        return False


def are_points_collinear_with_tolerance(p1, p2, p3, tolerance=1e-7):
    """Check if three points are collinear with some tolerance
    Args:
        p1 (list or np_array):
        p2 (list or np_array):
        p3 (list or np_array):
        tolerance (_type_, optional): Defaults to 1e-7.

    Returns:
        Boolean: True if collinear
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    cross_product = np.cross(p2 - p1, p3 - p1)
    norm_cycle = norm(cross_product)

    return norm_cycle < tolerance


def atoms_to_graph(atoms, nl, max_bond=0, max_bond_ratio=0):
    """
    Args:
        atoms (_type_): _description_
        nl (_type_): _description_
        max_bond (int, optional): _description_. Defaults to 0.
        max_bond_ratio (int, optional): _description_. Defaults to 0.

    Returns:
        _type_: _description_
    """
    if max_bond == 0 and max_bond_ratio == 0:
        raise RuntimeError("Please Specify bond information")

    g = nx.Graph()
    for index, atom in enumerate(atoms):
        if not g.has_node(node_symbol(atom)):
            g.add_node(node_symbol(atom), index=atom.index, symbol=atom.symbol)
        for neighbor, offset in zip(*nl.get_neighbors(index)):
            atom2 = atoms[neighbor]
            vector = atom.position - relative_position(atoms, neighbor, offset)
            distance = np.linalg.norm(vector)
            eqm_radii = covalent_radii[atom.number] + covalent_radii[atom2.number]
            check = max(max_bond, eqm_radii * max_bond_ratio)
            if distance > check:
                continue
            if not g.has_node(node_symbol(atom2)):
                g.add_node(node_symbol(atom2), index=atom2.index, symbol=atom2.symbol)
            if not g.has_edge(node_symbol(atom), node_symbol(atom2)):
                g.add_edge(
                    node_symbol(atom), node_symbol(atom2), weight=vector, start=index
                )
    return g
