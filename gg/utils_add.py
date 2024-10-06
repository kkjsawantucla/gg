""" General utilities for the modifiers """

from typing import Union
from itertools import combinations
import numpy as np
import networkx as nx
from numpy.linalg import norm
from ase import Atoms
from ase.data import covalent_radii, atomic_numbers
from gg.utils import (
    get_normals,
    add_ads,
    check_contact,
    is_within_tolerance,
    get_ref_pos_index,
)
from gg.utils_graph import node_symbol
from gg.utils_graph import is_cycle, are_points_collinear_with_tolerance


__author__ = "Kaustubh Sawant"


def generate_surf_sites(
    atoms: Atoms,
    graph: nx.Graph,
    index: list,
    surface_coord: Union[list, int],
) -> list:
    """
    Args:
        atoms (ase.Atoms):
        ads (ase.Atoms):
        graph (nx.Graph):
        index (list):
        surface_coord (list or int):

    Returns:
       valid surface site coord
    """
    if isinstance(surface_coord, int):
        surface_coord = [surface_coord]

    total_valid = []
    for s in surface_coord:
        # If the coordination > 1, enumerate different combinations of 2
        possible = list(combinations(index, s))
        valid = []

        # The for loop to ensure the combination of indices are valid
        for cycle in possible:
            if s == 1:
                valid.append(list(cycle))

            # Check of the 2 atoms selected are connected
            if s == 2:
                if graph.has_edge(
                    node_symbol(atoms[cycle[0]]), node_symbol(atoms[cycle[1]])
                ):
                    valid.append(list(cycle))

            # Check if the atoms are connected to each other in a cyclic manner
            if s == 3:
                nodes = [node_symbol(atoms[i]) for i in cycle]
                pos = [atoms[i].position for i in cycle]
                if is_cycle(graph, nodes):
                    # Some small unit cell can give collinear atoms as 3 adsorbate cycle
                    if not are_points_collinear_with_tolerance(
                        pos[0], pos[1], pos[2], tolerance=0.01
                    ):
                        valid.append(list(cycle))
        total_valid = total_valid + valid

    return total_valid


def generate_add_mono(
    atoms: Atoms,
    ads: Atoms,
    graph: nx.Graph,
    index: list,
    surf_coord: Union[list, int],
    ad_dist: Union[float, str] = 1.7,
    contact_error: float = 0.2,
) -> list:
    """
    Args:
        atoms (ase.Atoms):
        ads (ase.Atoms):
        graph (nx.Graph):
        index (list):
        surf_coord (list,int):
        ad_dist (float or str):  Defaults to 1.7.
        contact_error (float): _ Defaults to 0.2.

    Returns:
       ase.Atoms:
    """
    valid = generate_surf_sites(atoms, graph, index, surf_coord)
    movie = []
    # This for loops find the best position to add and generate atoms
    for cycle in valid:
        normal, ref_pos = get_normals(cycle, atoms, graph)
        offset = ref_pos
        unit_normal = normal / norm(normal)

        # If the adsorbae distance is a chemical symbol
        # It will try to find the best position according to covalent radii
        # Not very accurate
        if isinstance(ad_dist, str):
            if ad_dist in ads.get_chemical_symbols():
                for i in cycle:
                    current_dist = norm(offset - atoms[i].position)
                    req_dist = (
                        covalent_radii[atoms[i].number]
                        + covalent_radii[atomic_numbers[ad_dist]] * 0.9
                    )
                    if current_dist == 0:
                        offset += req_dist * unit_normal / len(cycle)
                    elif current_dist < req_dist:
                        proj_len = distance_point_to_line(
                            offset, unit_normal, atoms[i].position
                        )
                        offset += proj_len * unit_normal / len(cycle)

        elif isinstance(ad_dist, float):
            if ad_dist < np.average([covalent_radii[atoms[i].number] for i in cycle]):
                print("Issue in distance of adsorbate and substrate")
                continue
            else:
                offset += ad_dist * unit_normal
        ads_copy = ads.copy()
        ads_copy.rotate([0, 0, 1], normal, center=[0, 0, 0])
        atoms_copy = add_ads(atoms, ads_copy, offset=offset)

        # Make a final check if atoms are too close to each other
        if check_contact(atoms_copy, error=contact_error):
            continue
        else:
            movie.append(atoms_copy)
    return movie


def generate_add_bi(
    atoms: Atoms,
    ads: Atoms,
    graph: nx.Graph,
    surf_index: list,
    surf_coord: Union[list, int],
    ads_index: list,
    ad_dist: Union[float, str] = 1.7,
    ads_add_error: float = 0.2,
    contact_error: float = 0.2,
):
    """
    Args:
        atoms (ase.Atoms):
        ads (ase.Atoms):
        graph (nx.Graph):
        index (list):
        surf_coord (list,int):
        ad_dist (float or str):  Defaults to 1.7.
        contact_error (float): _ Defaults to 0.2.

    Returns:
       ase.Atoms:
    """
    movie = []
    # Get all possible sites
    valid_sites = generate_surf_sites(atoms, graph, surf_index, surf_coord)
    ads_bi_dist = norm(ads[ads_index[0]].position - ads[ads_index[1]].position)

    # Of all the sites , consider of pair of sites which can adsorb bidentate ligand/ads
    valid_bi_sites = []
    possible = list(combinations(valid_sites, 2))

    # Of all the pair of sites, make sure distance between them is similar to the ligand/ads dist
    for sites in possible:
        ref1 = get_ref_pos_index(sites[0], atoms, graph)
        ref2 = get_ref_pos_index(sites[1], atoms, graph)
        diff = norm(ref1 - ref2)
        if is_within_tolerance(diff, ads_bi_dist, ads_bi_dist * ads_add_error):
            valid_bi_sites.append(sites)

    # Adding the ligand/ads to the sites
    for sites in valid_bi_sites:
        ads_copy = ads.copy()
        cycle_1, cycle_2 = sites

        # Get the site information like the normal direction
        normal_1, offset_1 = get_normals(cycle_1, atoms, graph)
        normal_2, offset_2 = get_normals(cycle_2, atoms, graph)
        unit_normal_1 = normal_1 / norm(normal_1)
        unit_normal_2 = normal_2 / norm(normal_2)

        diff = norm(offset_1 - offset_2)

        if isinstance(ad_dist[0], str):
            offset_1 = get_offset(
                offset_1, ad_dist[0], cycle_1, unit_normal_1, atoms, ads_copy
            )
        elif isinstance(ad_dist[0], float):
            if ad_dist < np.average([covalent_radii[atoms[i].number] for i in cycle_1]):
                print("Issue in distance of adsorbate and substrate")
                continue
            else:
                offset_1 += ad_dist[0] * unit_normal_1

        if isinstance(ad_dist[1], str):
            offset_2 = get_offset(
                offset_2, ad_dist[1], cycle_2, unit_normal_2, atoms, ads_copy
            )
        elif isinstance(ad_dist[1], float):
            if ad_dist < np.average([covalent_radii[atoms[i].number] for i in cycle_2]):
                print("Issue in distance of adsorbate and substrate")
                continue
            else:
                offset_2 += ad_dist[1] * unit_normal_1
        target_vector = offset_2 - offset_1
        if diff > ads_bi_dist:
            target_pos = (
                offset_1
                + target_vector / norm(target_vector)
                + (diff - ads_bi_dist) / 2
            )
        else:
            target_pos = offset_2

        ads_copy = rotate_atoms_bi_along_vector(
            ads_copy, ads_index, target_vector, target_pos
        )
        atoms_copy = add_ads(atoms, ads_copy, offset=[0, 0, 0])
        # Make a final check if atoms are too close to each other
        if check_contact(atoms_copy, error=contact_error):
            continue
        else:
            movie.append(atoms_copy)
    return movie


def distance_point_to_line(p1: np.array, d: np.array, p0) -> float:
    """
    Args:
        p1 (np.array): point on the line
        d (np.array): direction of line (unit normal)
        p0 (np.array): fixed point

    Returns:
        float: distance of point from line
    """
    # Calculate the vector from the point on the line to the point
    p1_to_p0 = p0 - p1

    # Calculate the cross product of the direction vector and the vector from the line to the point
    cross_product = np.cross(d, p1_to_p0)
    distance = norm(cross_product) / norm(d)

    return distance


def rotate_bi(atoms: Atoms, index: list) -> Atoms:
    """
    Args:
        atoms (ase.Atoms):
        index (list): list of index of atoms that generate

    Returns:
        ase.Atoms:
    """
    # Get the position of the first atom
    first_atom_position = atoms.positions[index[0]]

    # Translate all atoms so that the first atom is at the origin
    new_positions = atoms.positions - first_atom_position
    atoms.set_positions(new_positions)

    # Calculate the vector from the first to the second atom
    second_atom_position = atoms.positions[index[1]]
    vector = second_atom_position - first_atom_position

    # Calculate the angle to rotate around the z-axis
    angle_z = np.degrees(np.arctan2(vector[1], vector[0]))
    atoms.rotate(-angle_z, "z")

    # Calculate the angle to rotate around the y-axis
    angle_y = np.degrees(np.arctan2(vector[2], vector[0]))
    atoms.rotate(-angle_y, "y")

    # Calculate the angle to rotate around the x-axis
    y_positions = atoms.positions[:, 1]
    z_positions = atoms.positions[:, 2]
    angle_x = np.degrees(np.arctan2(np.mean(z_positions), np.mean(y_positions)))
    atoms.rotate(-angle_x + 90, "x")

    if atoms[index[1]].position[0] < 0:
        atoms.translate((-atoms[index[1]].position[0], 0, 0))

    return atoms


def get_offset(offset, ad_dist, cycle, unit_normal, atoms, ads):
    """
    Args:
    Returns:
    """
    if ad_dist in ads.get_chemical_symbols():
        for i in cycle:
            current_dist = norm(offset - atoms[i].position)
            req_dist = (
                covalent_radii[atoms[i].number]
                + covalent_radii[atomic_numbers[ad_dist]] * 0.9
            )
            if current_dist == 0:
                offset += req_dist * unit_normal / len(cycle)
            elif current_dist < req_dist:
                proj_len = distance_point_to_line(
                    offset, unit_normal, atoms[i].position
                )
                offset += proj_len * unit_normal / len(cycle)
    return offset


def rotate_atoms_bi_along_vector(
    atoms: Atoms,
    index: list,
    target_vector: np.array,
    target_pos: np.array,
) -> Atoms:
    """
    Args:
        atoms (Atoms):
        index (list):
        target_vector (np.array):
        target_pos (np.array):

    Returns:
        ase.Atoms:
    """

    # Get the positions of the two atoms to be aligned
    pos1 = atoms.get_positions()[index[0]]
    pos2 = atoms.get_positions()[index[1]]

    current_vector = pos2 - pos1
    current_vector = current_vector / norm(current_vector)
    target_vector = target_vector / norm(target_vector)

    # Calculate the rotation angle along z axis
    rotation_angle = np.degrees(
        safe_arccos(np.dot(current_vector, [target_vector[0], target_vector[1], 0]))
    )
    if target_vector[1] < 0:
        rotation_angle = -rotation_angle

    if rotation_angle > 0.001:
        atoms.rotate(rotation_angle, "z")

    # Calculate the rotation angle along xy plane
    rotation_angle = np.degrees(
        safe_arccos(np.dot([target_vector[0], target_vector[1], 0], target_vector))
    )
    if target_vector[2] < 0:
        rotation_angle = -rotation_angle

    rotation_vector = [-target_vector[1], target_vector[0], 0]
    if rotation_angle > 0.001 and norm(rotation_vector) > 0.001:
        atoms.rotate(-rotation_angle, rotation_vector)

    # Move the atom to the target position
    first_atom_position = atoms.positions[index[0]]
    new_positions = atoms.positions + target_pos - first_atom_position
    atoms.set_positions(new_positions)

    return atoms


def safe_arccos(x: float):
    """
    Compute the arccosine of x if x is within the range [-1, 1].
    If x is greater than or equal to 1, return 0.
    If x is less than or equal to -1, return pi.
    """
    if x >= 1:
        return 0
    elif x <= -1:
        return np.pi
    elif -1 < x < 1:
        return np.arccos(x)
