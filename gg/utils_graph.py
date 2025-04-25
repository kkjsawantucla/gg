""" Utilities for graph manipulation """

from typing import Union
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from functools import partial
import networkx as nx
from networkx import Graph
from networkx.algorithms import weisfeiler_lehman_graph_hash
from networkx.algorithms.isomorphism import GraphMatcher
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ase.data import covalent_radii, chemical_symbols
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.data.colors import jmol_colors
from ase.build import molecule
from ase.collections import g2
from gg.data import adsorbates

__author__ = "Kaustubh Sawant"


def node_symbol(atom: Atoms) -> str:
    """Ensure symbol and index are converted to proper string types."""
    return f"{str(atom.symbol)}_{int(atom.index)}"


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


def node_match(n1: dict, n2: dict) -> bool:
    """
    Args:
        n1 (str):
        n2 (str):
    Returns:
        Boolean:
    """
    return n1["symbol"] == n2["symbol"]


def is_cycle(g: Graph, nodes: list) -> bool:
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
) -> Graph:
    """
    Args:
        atoms (ase.Atoms): an ase atoms object
        nl (ase.nl): Ase neighborlist
        max_bond (int, optional): . Defaults to 0.
        max_bond_ratio (int, optional): . Defaults to 0.

    Returns:
       (Graph):
    """
    if max_bond == 0 and max_bond_ratio == 0:
        raise RuntimeError("Please Specify bond information")
    g = Graph()
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


def get_graph_hash(graph: Graph) -> hash:
    """Generate a hash for the graph using node degrees and chemical symbols."""
    degrees = sorted(dict(graph.degree()).values())  # Degree sequence
    elements = sorted([data["symbol"] for _, data in graph.nodes(data=True)])
    return hash((tuple(degrees), tuple(elements)))


def get_wl_hash(graph: Graph, iterations: int = 3) -> hash:
    """Get weisfeiler_lehman_graph_hash"""
    return weisfeiler_lehman_graph_hash(
        graph, node_attr="symbol", iterations=iterations
    )


def get_unique_graph_indexes(
    graph_list: list, unique_method: str = "fullgraph", depth: int = 3
) -> list:
    """
    Args:
        graph_list (list): list[Graph]

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
            else:
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
        list[Graph]: converted list of graphs
    """
    graph_list = []
    for atoms in list_atoms:
        graph_list.append(_process_single_atoms(atoms, max_bond, max_bond_ratio))
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
    graph: Graph, graph_list: list, comp_type: str = "fullgraph"
) -> bool:
    """Check if the given graph is not isomorphic to any graph in the list.

    Args:
    graph (Graph): The graph to check.
    graph_list (list of Graph): The list of graphs to compare against.

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
            else:
                if compare_subgraph_uniqueness(
                    graph, unique_graph, center_symbol=comp_type
                ):
                    return False
    return True


def draw_graph(
    graph: Graph, graph_type: str = "none", atoms: Atoms = None, **kwargs
) -> plt.Figure:
    """Draw an atomic graph with different layout options, displaying node degrees.

    Args:
        graph (Graph): Graph to draw.
        graph_type (str, optional): Layout type ("none", "circular", "kamada_kawai").
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


def get_connecting_nodes(graph: Graph, cluster_ind: list, atoms: Atoms) -> list:
    """
    Find nodes that connect the cluster to the rest of the graph and return their 'index' attribute.

    Parameters:
    graph (Graph): The input graph.
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


def replace_subgraph_with_node(
    graph: nx.Graph, subgraph: nx.Graph, new_node_label: str
) -> nx.Graph:
    """
    Find instances of a subgraph and replace them with a single node while preserving connectivity.

    Args:
        graph (nx.Graph): The input graph.
        subgraph (nx.Graph): The subgraph to be replaced.
        new_node_label (str): The label of the new node that replaces the subgraph.

    Returns:
        nx.Graph: A modified graph with the subgraph replaced by a new node.
    """
    temp_graph = graph.copy()
    matcher = GraphMatcher(temp_graph, subgraph, node_match=node_match)
    matched_subgraphs = []

    # Step 1: Collect all subgraph matches before making modifications
    for match in matcher.subgraph_isomorphisms_iter():
        matched_subgraphs.append(set(match.keys()))  # Store node sets

    nodes_to_remove = set()
    replacement_mapping = {}

    for subgraph_nodes in matched_subgraphs:
        # Create the new replacement node
        new_node = f"{new_node_label}_{len(temp_graph.nodes)}"
        temp_graph.add_node(new_node, symbol=new_node_label)

        # Identify external edges (edges connecting subgraph to the rest of the graph)
        external_edges = []
        for node in subgraph_nodes:
            for neighbor in temp_graph.neighbors(node):
                if neighbor not in subgraph_nodes:
                    external_edges.append((new_node, neighbor))

        # Store the mapping and mark nodes for removal later
        replacement_mapping[frozenset(subgraph_nodes)] = new_node
        nodes_to_remove.update(subgraph_nodes)

        # Add external edges to the new node
        temp_graph.add_edges_from(external_edges)

    # Step 3: Remove all matched subgraph nodes after replacements
    for nodes in replacement_mapping:
        temp_graph.remove_nodes_from(nodes)

    return temp_graph


def generate_centered_subgraphs(
    graph: Graph, center_symbol: str, depth: int
) -> list[Graph]:
    """
    Generate subgraphs centered around nodes with `center_symbol`,
    up to `depth` levels using ego_graph.
    """
    if isinstance(center_symbol, str):
        center_symbol = [center_symbol]

    center_symbol = set(center_symbol)

    if not all(item in chemical_symbols for item in center_symbol):
        for symbol in center_symbol:
            if symbol in adsorbates:
                atoms = adsorbates[symbol]
            elif symbol in g2.names:
                atoms = molecule(symbol)
            else:
                raise RuntimeError("Cannot generate graphs")

            # Make the ase.NeighborList
            nl = NeighborList(
                natural_cutoffs(atoms), self_interaction=False, bothways=True
            )
            nl.update(atoms)
            subgraph = atoms_to_graph(atoms, nl, max_bond_ratio=1.2)
            graph = replace_subgraph_with_node(graph, subgraph, symbol)

    center_nodes = [
        node for node, data in graph.nodes(data=True) if data["symbol"] in center_symbol
    ]
    subgraphs = [nx.ego_graph(graph, node, radius=depth) for node in center_nodes]

    return subgraphs


def compare_subgraph_uniqueness(
    graph1: Graph, graph2: Graph, center_symbol: str, depth: int = 3
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
                and GraphMatcher(sg1, sg2, node_match=node_match).is_isomorphic()
            ):
                matched.add(i)
                break

    # If all subgraphs from graph1 found a match in graph2, they are equivalent
    return len(matched) == len(subgraphs1)


def compare_fullgraph_uniqueness(
    graph1: Graph,
    graph2: Graph,
) -> bool:
    """Compare Full Graphs"""

    return GraphMatcher(graph1, graph2, node_match=node_match).is_isomorphic()


def _process_single_atoms(atoms, max_bond, max_bond_ratio):
    """Top-level function for parallel processing (must be picklable)."""
    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    return atoms_to_graph(atoms, nl, max_bond=max_bond, max_bond_ratio=max_bond_ratio)


def list_atoms_to_graphs_parallel(
    list_atoms: list,
    max_bond: float = 0,
    max_bond_ratio: float = 0,
    num_workers: int = 1,
) -> list:
    """Parallel version of atoms-to-graph conversion."""
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        graph_list = list(
            executor.map(
                _process_single_atoms,
                list_atoms,
                repeat(max_bond),
                repeat(max_bond_ratio),
            )
        )
    return graph_list


def get_unique_atoms_parallel(
    movie: list,
    max_bond: float = 0,
    max_bond_ratio: float = 0,
    unique_method: str = "fullgraph",
    depth: int = 3,
    num_workers: int = 1,
) -> list:
    """ """
    graph_list = list_atoms_to_graphs_parallel(
        movie, max_bond=max_bond, max_bond_ratio=max_bond_ratio, num_workers=num_workers
    )
    unique_indexes = get_unique_graph_indexes_parallel(
        graph_list, unique_method=unique_method, depth=depth, num_workers=num_workers
    )
    return [movie[i] for i in unique_indexes]


def get_unique_graph_indexes_parallel(
    graph_list: list,
    unique_method: str = "fullgraph",
    depth: int = 3,
    num_workers: int = 1,
) -> list:
    """Parallelized version using batched hash precomputation and parallel isomorphism checks."""
    # Validate all graphs have string node keys before processing
    for g in graph_list:
        for node in g.nodes():
            if not isinstance(node, str):
                raise ValueError(
                    f"Invalid node key {node} (type {type(node)}) found in graph"
                )

    # Precompute all graph hashes in parallel
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        graph_hashes = list(executor.map(get_graph_hash, graph_list))

    unique_graphs = []  # Stores (graph, graph_hash) tuples
    seen_hashes = set()
    unique_indexes = []

    for index, (graph, gh) in enumerate(zip(graph_list, graph_hashes)):
        if gh not in seen_hashes:
            # New unique graph found
            seen_hashes.add(gh)
            unique_graphs.append((graph, gh))
            unique_indexes.append(index)
        else:
            # Find all candidates with matching hash
            candidates = [ug for ug, ug_h in unique_graphs if ug_h == gh]
            # Parallel isomorphism checks against all candidates
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                checks = list(
                    executor.map(
                        _check_isomorphism,
                        candidates,
                        repeat(graph),
                        repeat(unique_method),
                        repeat(depth),
                    )
                )
                print(checks)

            if not any(checks):
                # No matches found - add to unique
                unique_graphs.append((graph, gh))
                unique_indexes.append(index)
    return unique_indexes


# Worker function for parallel isomorphism checks
def _check_isomorphism(candidate_graph, current_graph, method, depth):
    if method == "fullgraph":
        return compare_fullgraph_uniqueness(current_graph, candidate_graph)
    else:
        return compare_subgraph_uniqueness(
            current_graph, candidate_graph, center_symbol=method, depth=depth
        )
