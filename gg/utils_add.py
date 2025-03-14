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
                        pos[0], pos[1], pos[2], tolerance=0.05
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
    method: str = "svd",
    tag: bool = True,
) -> list:
    """
    Args:
        atoms (ase.Atoms):
        ads (ase.Atoms):
        graph (nx.Graph):
        index (list):
        surf_coord (list,int):
        ad_dist (float or str): Defaults to 1.7.
        contact_error (float): Defaults to 0.2.
        method (str): Defaults to "svd".

    Returns:
       ase.Atoms:
    """
    valid = generate_surf_sites(atoms, graph, index, surf_coord)
    movie = []
    # This for loops find the best position to add and generate atoms

    for cycle in valid:
        normal, ref_pos = get_normals(cycle, atoms, graph, method=method)
        offset = ref_pos
        unit_normal = normal / norm(normal)

        # If the adsorbae distance is a chemical symbol
        # It will try to find the best position according to covalent radii
        # Not very accurate
        if isinstance(ad_dist, str):
            offset = get_offset(offset, ad_dist, cycle, unit_normal, atoms, ads)

        elif isinstance(ad_dist, float) or isinstance(ad_dist, int):
            if ad_dist < np.average([covalent_radii[atoms[i].number] for i in cycle]):
                print("Issue in distance of adsorbate and substrate")
                continue
            else:
                offset = get_average_cycle_position(cycle, atoms)
                if norm(offset - atoms[cycle[0]].position) < 0.01 and len(cycle) != 1:
                    continue
                offset += ad_dist * unit_normal
        ads_copy = ads.copy()
        ads_copy.rotate([0, 0, 1], normal, center=[0, 0, 0])
        atoms_copy = add_ads(atoms, ads_copy, offset=offset, tag=tag)

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
    method: str = "svd",
    tag: bool = True,
):
    """
    Args:
        atoms (ase.Atoms):
        ads (ase.Atoms):
        graph (nx.Graph):
        index (list):
        surf_coord (list,int):
        ad_dist (float or str):  Defaults to 1.7.
        contact_error (float): Defaults to 0.2.

    Returns:
       ase.Atoms:
    """
    movie = []
    # Get all possible sites
    valid_sites = generate_surf_sites(atoms, graph, surf_index, surf_coord)
    ads_bi_dist = norm(ads[ads_index[0]].position - ads[ads_index[1]].position)
    cell = atoms.cell[:]

    # Of all the sites , consider of pair of sites which can adsorb bidentate ligand/ads
    valid_bi_sites = []
    possible = list(combinations(valid_sites, 2))

    # Of all the pair of sites, make sure distance between them is similar to the ligand/ads dist
    for sites in possible:
        ref1 = get_ref_pos_index(sites[0], atoms, graph)
        ref2 = get_ref_pos_index(sites[1], atoms, graph)
        diff = minimum_image_distance(ref1, ref2, cell)
        if is_within_tolerance(diff, ads_bi_dist, ads_bi_dist * ads_add_error):
            valid_bi_sites.append(sites)

    # Adding the ligand/ads to the sites
    for sites in valid_bi_sites:
        ads_copy = ads.copy()
        cycle_1, cycle_2 = sites

        # Get the site information like the normal direction
        normal_1, offset_1 = get_normals(cycle_1, atoms, graph, method=method)
        normal_2, offset_2 = get_normals(cycle_2, atoms, graph, method=method)
        normal_t, _ = get_normals(
            list(set(cycle_2 + cycle_1)), atoms, graph, method=method
        )
        unit_normal_1 = normal_1 / norm(normal_1)
        unit_normal_2 = normal_2 / norm(normal_2)

        diff = minimum_image_distance(offset_1, offset_2, cell)
        diff_og = norm(offset_2 - offset_1)

        if isinstance(ad_dist[0], str):
            offset_1 = get_offset(
                offset_1, ad_dist[0], cycle_1, unit_normal_1, atoms, ads_copy
            )

        elif isinstance(ad_dist[0], float) or isinstance(ad_dist[0], int):
            if ad_dist[0] < np.average([covalent_radii[atoms[i].number] for i in cycle_1]):
                print("Issue in distance of adsorbate and substrate")
                continue
            else:
                offset_1 = get_average_cycle_position(cycle_1, atoms)
                if (
                    norm(offset_1 - atoms[cycle_1[0]].position) < 0.01
                    and len(cycle_1) != 1
                ):
                    continue
                offset_1 += ad_dist[0] * unit_normal_1

        if isinstance(ad_dist[1], str):
            offset_2 = get_offset(
                offset_2, ad_dist[1], cycle_2, unit_normal_2, atoms, ads_copy
            )
        elif isinstance(ad_dist[1], float) or isinstance(ad_dist[1], int):
            if ad_dist[1] < np.average([covalent_radii[atoms[i].number] for i in cycle_2]):
                print("Issue in distance of adsorbate and substrate")
                continue
            else:
                offset_2 = get_average_cycle_position(cycle_2, atoms)
                if (
                    norm(offset_2 - atoms[cycle_2[0]].position) < 0.01
                    and len(cycle_2) != 1
                ):
                    continue
                offset_2 += ad_dist[1] * unit_normal_1
        if abs(diff - diff_og) > 0.001:
            target_vector = minimum_image_distance(
                offset_1, offset_2, cell, get_norm=False
            )
        else:
            target_vector = offset_2 - offset_1
        if diff > ads_bi_dist:
            target_pos = (
                offset_1
                + target_vector / norm(target_vector) * (diff - ads_bi_dist) / 2
            )
        else:
            target_pos = offset_1
        ads_copy = rotate_atoms_bi_along_vector(
            ads_copy, ads_index, target_vector, target_pos, normal_t
        )
        atoms_copy = add_ads(atoms, ads_copy, offset=[0, 0, 0], tag=tag)
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
    vector = atoms.positions[index[1]] - atoms.positions[index[0]]

    # Calculate the angle to rotate around the z-axis
    angle_z = angle_between([vector[0], vector[1]], [1, 0])
    if vector[1] > 0:
        atoms.rotate(-angle_z, "z")
    else:
        atoms.rotate(angle_z, "z")

    vector = atoms.positions[index[1]] - atoms.positions[index[0]]
    # Calculate the angle to rotate around the y-axis
    angle_y = angle_between([vector[0], vector[2]], [1, 0])
    if vector[2] > 0:
        atoms.rotate(angle_y, "y")
    else:
        atoms.rotate(-angle_y, "y")

    # Calculate the angle to rotate around the x-axis
    y_positions = atoms.positions[:, 1]
    z_positions = atoms.positions[:, 2]
    angle_x = np.degrees(np.arctan2(np.mean(z_positions), np.mean(y_positions)))
    atoms.rotate(-angle_x + 90, "x")

    return atoms


def get_offset(offset, ad_dist, cycle, unit_normal, atoms, ads):
    """
    Args:
    Returns:
    """
    if ad_dist in ads.get_chemical_symbols():
        for i in cycle:
            current_dist = minimum_image_distance(
                offset, atoms[i].position, atoms.cell
            )
            req_dist = (
                covalent_radii[atoms[i].number]
                + covalent_radii[atomic_numbers[ad_dist]] * 0.9
            )
            if current_dist < 0.001:
                offset += req_dist * unit_normal / len(cycle)
            elif current_dist < req_dist:
                clos_pos = closest_periodic_image(
                    offset, atoms[i].position, atoms.cell[0], atoms.cell[1]
                )
                proj_len = distance_point_to_line(offset, unit_normal, clos_pos)
                offset += proj_len * unit_normal / len(cycle)
    return offset


def rotate_atoms_bi_along_vector(
    atoms: Atoms,
    index: list,
    target_vector: np.array,
    target_pos: np.array,
    new_norm: np.array,
) -> Atoms:
    """
    Args:
        atoms (Atoms): Input adsorbate to rotate
        index (list): Index of the atoms in adosrbate that form bond
        target_vector (np.array): Target vector along which rotated
        target_pos (np.array): Move to the new x,y,x position
        new_norm (np.array): Nprmal of atoms

    Returns:
        ase.Atoms: rotated atoms
    """

    # Get the positions of the two atoms to be aligned
    pos1 = atoms.get_positions()[index[0]]
    pos2 = atoms.get_positions()[index[1]]

    current_vector = pos2 - pos1
    current_vector = current_vector / norm(current_vector)
    target_vector = target_vector / norm(target_vector)
    if norm(new_norm) > 0.001:
        new_norm = new_norm / norm(new_norm)

    # Calculate the rotation angle along z axis
    if abs(target_vector[0]) > 0.001 or abs(target_vector[1]) > 0.001:
        rotation_angle = angle_between(
            current_vector, [target_vector[0], target_vector[1], 0]
        )
        if target_vector[1] < 0:
            rotation_angle = -rotation_angle

        if abs(rotation_angle) > 0.001:
            atoms.rotate(rotation_angle, "z")

    # Calculate the rotation angle along xy plane
    if abs(target_vector[2]) > 0.001:
        rotation_angle = angle_between(
            [target_vector[0], target_vector[1], 0], target_vector
        )
        if target_vector[2] < 0:
            rotation_angle = -rotation_angle
        if abs(target_vector[0]) < 0.001 and abs(target_vector[1]) < 0.001:
            rotation_vector = [0, 1, 0]
        else:
            rotation_vector = [-target_vector[1], target_vector[0], 0]
        if abs(rotation_angle) > 0.001:
            atoms.rotate(-rotation_angle, rotation_vector)

    if norm(new_norm) > 0.001:
        positions = [0, 0, 1]
        # Rotate along the target vector
        rotation_angle = angle_between(positions, new_norm)
        cross_product = np.cross(positions, new_norm)
        if norm(cross_product) < 0:
            rotation_angle = 360 - rotation_angle
        atoms.rotate(rotation_angle, target_vector)

    # Move the atom to the target position
    first_atom_position = atoms.positions[index[0]]
    new_positions = atoms.positions + target_pos - first_atom_position
    atoms.set_positions(new_positions)

    return atoms


def rotate_mono(atoms: Atoms, index: int) -> Atoms:
    """
    Args:
        atoms (ase.Atoms):
        index (list): list of index of atoms that generate

    Returns:
        ase.Atoms:
    """
    # Get the position of the first atom
    first_atom_position = atoms.positions[index]

    # Translate all atoms so that the first atom is at the origin
    new_positions = atoms.positions - first_atom_position
    atoms.set_positions(new_positions)

    y_positions = atoms.positions[:, 1]
    z_positions = atoms.positions[:, 2]
    angle_x = np.degrees(np.arctan2(np.mean(z_positions), np.mean(y_positions)))
    atoms.rotate(-angle_x + 90, "x")

    x_positions = atoms.positions[:, 0]
    z_positions = atoms.positions[:, 2]
    angle_x = np.degrees(np.arctan2(np.mean(z_positions), np.mean(x_positions)))
    atoms.rotate(-90 + angle_x, "y")

    return atoms


def angle_between(v1, v2):
    """_summary_
    Args:
        v1 (array): _description_
        v2 (array): _description_

    Returns:
        _type_: angle in degrees
    """
    # Compute dot product and magnitudes
    dot_product = np.dot(v1, v2)
    magnitude_v1 = norm(v1)
    magnitude_v2 = norm(v2)

    # Compute cosine of the angle
    cos_angle = dot_product / (magnitude_v1 * magnitude_v2)

    # Ensure the cosine value is in range [-1, 1] to avoid numerical errors
    cos_angle = np.clip(cos_angle, -1.0, 1.0)

    # Compute angle in radians and then convert to degrees
    angle_radians = np.arccos(cos_angle)
    angle_degrees = np.degrees(angle_radians)

    return angle_degrees


def minimum_image_distance(coord1, coord2, cell, get_norm=True):
    """
    Calculate the minimum image distance between two 3D coordinates in a periodic cell.

    Parameters:
    coord1 (tuple): A tuple representing the first coordinate (x1, y1, z1).
    coord2 (tuple): A tuple representing the second coordinate (x2, y2, z2).
    cell (array-like): A 3x3 array representing the periodic cell vectors.

    Returns:
    float: The minimum image distance between the two coordinates.
    """
    coord1 = np.array(coord1, dtype=np.float64)
    coord2 = np.array(coord2, dtype=np.float64)
    cell = np.array(cell, dtype=np.float64)
    delta = coord2 - coord1
    delta -= np.round(delta @ np.linalg.inv(cell)) @ cell
    if get_norm:
        return norm(delta)
    else:
        return delta


def get_average_cycle_position(cycle, atoms):
    """_summary_
    Args:
        cycle (list(int)):
        atoms (ase.Atoms):
    Returns:
        np.array:
    """
    init_pos = atoms[cycle[0]].position
    corr = 0
    for i in cycle[1:]:
        corr += minimum_image_distance(
            init_pos, atoms[i].position, cell=atoms.cell, get_norm=False
        )

    return init_pos + corr / 3


def closest_periodic_image(reference, point, x_vec, y_vec):
    """
    Finds the closest periodic image of 'point' to 'reference' within a 2D periodic cell in 3D space.

    Parameters:
        reference (array-like): The fixed reference point [x, y, z].
        point (array-like): The point whose closest periodic image is sought [x, y, z].
        x_vec (array-like): The x-direction periodic cell vector [x, y, z].
        y_vec (array-like): The y-direction periodic cell vector [x, y, z].

    Returns:
        np.ndarray: The coordinates of the closest periodic image.
    """
    reference = np.array(reference)
    point = np.array(point)
    x_vec = np.array(x_vec)
    y_vec = np.array(y_vec)

    # Generate translations along periodic directions
    translations = [
        i * x_vec + j * y_vec
        for i in range(-1, 2)  # Checking -1, 0, 1 in x direction
        for j in range(-1, 2)  # Checking -1, 0, 1 in y direction
    ]

    # Generate all periodic images
    periodic_images = [point + t for t in translations]

    # Find the closest image
    closest_image = min(
        periodic_images, key=lambda img: np.linalg.norm(img - reference)
    )

    return closest_image
