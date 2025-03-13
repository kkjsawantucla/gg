""" Utilities for graph manipulation """

from typing import Union
import networkx as nx
from networkx.algorithms import weisfeiler_lehman_graph_hash
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
            atom2 = atoms[neighbor]
            vector = atom.position - relative_position(atoms, neighbor, offset)
            distance = np.linalg.norm(vector)
            eqm_radii = covalent_radii[atom.number] + covalent_radii[atom2.number]
            check = max(max_bond, eqm_radii * max_bond_ratio)
            if distance > check:
                continue
            index_n.append(neighbor)
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
    """Generate a hash for the graph using node degrees and chemical symbols."""
    degrees = sorted(dict(graph.degree()).values())  # Degree sequence
    elements = sorted([data["symbol"] for _, data in graph.nodes(data=True)])
    return hash((tuple(degrees), tuple(elements)))


def get_wl_hash(graph: nx.Graph, iterations: int = 3) -> hash:
    """Get weisfeiler_lehman_graph_hash"""
    return weisfeiler_lehman_graph_hash(
        graph, node_attr="symbol", iterations=iterations
    )


def get_unique_graph_indexes(
    graph_list: list, unique_method: str = "fullgraph", depth: int = 3
) -> list:
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
            unique_graphs.append(graph)
            unique_indexes.append(index)
        else:
            if unique_method == "fullgraph":
                # Perform a full isomorphism check to confirm the uniqueness
                if not any(
                    compare_fullgraph_uniqueness(graph, unique_graph)
                    for unique_graph in unique_graphs
                ):
                    unique_graphs.append(graph)
                    unique_indexes.append(index)
            if unique_method in chemical_symbols:
                # Perform a full isomorphism check to confirm the uniqueness
                if not any(
                    compare_subgraph_uniqueness(
                        graph, unique_graph, center_symbol=unique_method, depth=depth
                    )
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
    movie: list,
    max_bond: float = 0,
    max_bond_ratio: float = 0,
    unique_method: str = "fullgraph",
    depth: int = 3,
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
    unique_indexes = get_unique_graph_indexes(
        graph_list, unique_method=unique_method, depth=depth
    )
    return [movie[i] for i in unique_indexes]


def is_unique_graph(
    graph: nx.Graph, graph_list: list, comp_type: str = "fullgraph"
) -> bool:
    """Check if the given graph is not isomorphic to any graph in the list.

    Args:
    graph (nx.Graph): The graph to check.
    graph_list (list of nx.Graph): The list of graphs to compare against.

    Returns:
    bool: True if the graph is unique, False otherwise.
    """
    graph_hash = get_graph_hash(graph)
    for unique_graph in graph_list:
        unique_graph_hashes = get_graph_hash(unique_graph)
        if graph_hash != unique_graph_hashes:
            continue
        else:
            if comp_type == "fullgraph":
                if compare_fullgraph_uniqueness(graph, unique_graph):
                    return False
            elif comp_type in chemical_symbols:
                if compare_subgraph_uniqueness(
                    graph, unique_graph, center_symbol=comp_type
                ):
                    return False
    return True


def draw_graph(
    graph: nx.Graph, graph_type: str = "none", atoms: Atoms = None, **kwargs
) -> plt.Figure:
    """Draw an atomic graph with different layout options, displaying node degrees.

    Args:
        graph (nx.Graph): Graph to draw.
        graph_type (str, optional): Layout type ("none", "circular", "kamada_kawai"). Defaults to "none".
        **kwargs: Additional keyword arguments for the NetworkX draw functions.

    Returns:
        plt.Figure: The created matplotlib figure.
    """
    color = []
    edgecolors = []

    for node in graph.nodes(data=True):
        symbol = node[1].get("symbol")  # Default to carbon if symbol is missing
        color.append(jmol_colors[chemical_symbols.index(symbol)])
        edgecolors.append("black")

    fig, ax = plt.subplots(figsize=(10, 10))

    if graph_type == "circular":
        pos = nx.circular_layout(graph)
    elif graph_type == "kamada_kawai":
        pos = nx.kamada_kawai_layout(graph)
    elif graph_type == "atoms":
        if isinstance(atoms, Atoms):
            pos = {
                i[0]: (
                    atoms.positions[i[1]["index"]][0],
                    atoms.positions[i[1]["index"]][2],
                )
                for i in graph.nodes(data=True)
            }
        else:
            raise RuntimeError("Please provide ase.Atoms object")
    else:
        pos = nx.spring_layout(graph, seed=4)

    nx.draw(graph, pos, node_color=color, edgecolors=edgecolors, ax=ax, **kwargs)

    # Add node degrees as labels
    node_degrees = dict(graph.degree())
    for node, (x, y) in pos.items():
        ax.text(
            x,
            y,
            str(node_degrees[node]),
            fontsize=12,
            ha="center",
            va="center",
            color="black",
        )
    plt.show()
    return fig


def get_connecting_nodes(graph: nx.Graph, cluster_ind: list, atoms: Atoms) -> list:
    """
    Find nodes that connect the cluster to the rest of the graph and return their 'index' attribute.

    Parameters:
    graph (nx.Graph): The input graph.
    cluster_nodes (list): The list of nodes in the cluster.

    Returns:
    list: A list of 'index' are part of the connecting edges to the cluster.
    """
    cluster_nodes = [node_symbol(atoms[i]) for i in cluster_ind]
    cluster = set(cluster_nodes)
    connecting_nodes = set()

    for node in cluster:
        for neighbor in graph.neighbors(node):
            if neighbor not in cluster:
                connecting_nodes.add(neighbor)

    # Get the 'index' attribute of each connecting node
    connecting_indices = [graph.nodes[node]["index"] for node in connecting_nodes]
    return connecting_indices


def generate_centered_subgraphs(
    graph: nx.Graph, center_symbol: str, depth: int
) -> list[nx.Graph]:
    """
    Generate subgraphs centered around nodes with `center_symbol`,
    up to `depth` levels using ego_graph.
    """
    # Find all center nodes with the specified symbol
    center_nodes = [
        node for node, data in graph.nodes(data=True) if data["symbol"] == center_symbol
    ]

    # Generate ego graphs (subgraphs) for each center node
    subgraphs = [nx.ego_graph(graph, node, radius=depth) for node in center_nodes]

    return subgraphs


def compare_subgraph_uniqueness(
    graph1: nx.Graph, graph2: nx.Graph, center_symbol: str, depth: int = 3
) -> bool:
    """Compare uniqueness of two graphs by checking their subgraphs first."""

    subgraphs1 = generate_centered_subgraphs(graph1, center_symbol, depth)
    subgraphs2 = generate_centered_subgraphs(graph2, center_symbol, depth)

    # If different number of subgraphs, they are not identical
    if len(subgraphs1) != len(subgraphs2):
        return False

    # Compare subgraphs for isomorphism
    matched = set()
    for sg1 in subgraphs1:
        for i, sg2 in enumerate(subgraphs2):
            if (
                i not in matched
                and nx.algorithms.isomorphism.GraphMatcher(
                    sg1, sg2, node_match=node_match
                ).is_isomorphic()
            ):
                matched.add(i)
                break

    # If all subgraphs from graph1 found a match in graph2, they are equivalent
    return len(matched) == len(subgraphs1)


def compare_fullgraph_uniqueness(
    graph1: nx.Graph,
    graph2: nx.Graph,
) -> bool:
    """_summary_"""

    result = nx.algorithms.isomorphism.GraphMatcher(
        graph1, graph2, node_match=node_match
    ).is_isomorphic()

    return result
