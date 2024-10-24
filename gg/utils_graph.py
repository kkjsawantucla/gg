""" Utilities for graph manipulation """

from typing import Union
import networkx as nx
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ase.data import covalent_radii, chemical_symbols
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.data.colors import jmol_colors

__author__ = "Kaustubh Sawant"


def node_symbol(atom: Atoms) -> str:
    """
    Args:
       atom (Atoms Object): atoms to convert

    Returns:
       str: "symbol_index"
    """
    return f"{atom.symbol}_{atom.index}"


def relative_position(atoms: Atoms, neighbor: int, offset: np.array) -> np.array:
    """
    Args:
         atoms (ase.Atoms):
         neighbor (int): Index of the neighbor
         offset (array):

     Returns:
       np.array: position of neighbor wrt to offset
    """
    return atoms[neighbor].position + np.dot(offset, atoms.get_cell())

def node_match(n1: str, n2: str) -> bool:
    """
    Args:
        n1 (str):
        n2 (str):
    Returns:
        Boolean:
    """
    return n1["symbol"] == n2["symbol"]


def is_cycle(g: nx.Graph, nodes: list) -> bool:
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


def are_points_collinear_with_tolerance(
    p1: Union[list, np.array],
    p2: Union[list, np.array],
    p3: Union[list, np.array],
    tolerance: float = 1e-7,
) -> bool:
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


def atoms_to_graph(
    atoms: Atoms, nl, max_bond: float = 0, max_bond_ratio: float = 0
) -> nx.Graph:
    """
    Args:
        atoms (ase.Atoms): an ase atoms object
        nl (ase.nl): Ase neighborlist
        max_bond (int, optional): . Defaults to 0.
        max_bond_ratio (int, optional): . Defaults to 0.

    Returns:
       (nx.Graph):
    """
    if max_bond == 0 and max_bond_ratio == 0:
        raise RuntimeError("Please Specify bond information")

    g = nx.Graph()
    for index, atom in enumerate(atoms):
        index_n = []
        if not g.has_node(node_symbol(atom)):
            g.add_node(node_symbol(atom), index=atom.index, symbol=atom.symbol)
        for neighbor, offset in zip(*nl.get_neighbors(index)):
            index_n.append(neighbor)
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
                    node_symbol(atom),
                    node_symbol(atom2),
                    weight=distance,
                    weight2=vector,
                    start=index,
                )
        if len(index_n) != len(set(index_n)):
            raise RuntimeError(
                "Two atoms connected multiple times! unit cell is too small to make graphs"
            )
    return g


def get_graph_hash(graph: nx.Graph) -> hash:
    """
    Args:
        graph (nx.Graph): Graph to be hashed

    Returns:
        hash : get hashed representation of the graph
    """
    # Create a sorted tuple of node attributes and edges to generate a unique hash
    node_attrs = tuple(
        sorted((node, data["symbol"]) for node, data in graph.nodes(data=True))
    )
    edges = tuple(sorted(graph.edges()))
    return hash((node_attrs, edges))


def get_unique_graph_indexes(graph_list: list) -> list:
    """
    Args:
        graph_list (list): list[nx.Graph]

    Returns:
        list: get list of unique graph indexes
    """
    unique_graphs = []
    seen_hashes = set()
    unique_indexes = []
    for index, graph in enumerate(graph_list):
        graph_hash = get_graph_hash(graph)
        if graph_hash not in seen_hashes:
            seen_hashes.add(graph_hash)
            # Perform a full isomorphism check to confirm the uniqueness
            if not any(
                nx.algorithms.isomorphism.GraphMatcher(
                    graph, unique_graph, node_match=node_match
                ).is_isomorphic()
                for unique_graph in unique_graphs
            ):
                unique_graphs.append(graph)
                unique_indexes.append(index)
    return unique_indexes


def list_atoms_to_graphs(
    list_atoms: list, max_bond: float = 0, max_bond_ratio: float = 0
) -> list:
    """
    Args:
        list_atoms (list[Atoms]): list of atoms to convert to graphs
        max_bond (int, optional): . Defaults to 0.
        max_bond_ratio (int, optional): . Defaults to 0.

    Returns:
        list[nx.Graph]: converted list of graphs
    """
    graph_list = []
    for atoms in list_atoms:
        nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
        nl.update(atoms)
        graph_list.append(
            atoms_to_graph(atoms, nl, max_bond=max_bond, max_bond_ratio=max_bond_ratio)
        )
    return graph_list


def get_unique_atoms(
    movie: list, max_bond: float = 0, max_bond_ratio: float = 0
) -> list:
    """
    Args:
        movie (list[Atoms]): _description_
        max_bond (float, optional): _description_. Defaults to 0.
        max_bond_ratio (float, optional): _description_. Defaults to 0.

    Returns:
        list[Atoms]: _description_
    """
    graph_list = list_atoms_to_graphs(
        movie, max_bond=max_bond, max_bond_ratio=max_bond_ratio
    )
    unique_indexes = get_unique_graph_indexes(graph_list)
    return [movie[i] for i in unique_indexes]


def is_unique_graph(graph: nx.Graph, graph_list: list) -> bool:
    """Check if the given graph is not isomorphic to any graph in the list.

    Args:
    graph (nx.Graph): The graph to check.
    graph_list (list of nx.Graph): The list of graphs to compare against.

    Returns:
    bool: True if the graph is unique, False otherwise.
    """
    for unique_graph in graph_list:
        if nx.algorithms.isomorphism.GraphMatcher(
            graph, unique_graph, node_match=node_match
        ).is_isomorphic():
            return False
    return True


def draw_graph(graph: nx.Graph, graph_type: str = "none", **kwargs) -> plt.figure:
    """Draw atoms graph

    Args:
        graph (nx.Graph): Graph to draw
        graph_type (str, optional): Defaults to "none".
    """
    color = []
    edgecolors = []
    for node in graph.nodes(data=True):
        symbol = node[1]["symbol"]
        color.append(jmol_colors[chemical_symbols.index(symbol)])
        edgecolors.append("black")
    if graph_type == "circular":
        nx.draw_circular(graph, node_color=color, edgecolors=edgecolors, **kwargs)
    if graph_type == "kamada_kawai":
        nx.draw_kamada_kawai(graph, node_color=color, edgecolors=edgecolors, **kwargs)
    else:
        layout = nx.spring_layout(graph, seed=4)
        plt.figure(figsize=(10, 10))
        # labels = nx.get_edge_attributes(graph, "weight")
        # nx.draw_networkx_edge_labels(graph, pos = layout, edge_labels=labels)
        nx.draw(graph, pos=layout, node_color=color, edgecolors=edgecolors, **kwargs)
    plt.draw()
    plt.show()
