""" Cluster """

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from scipy.linalg import svd

def fit_plane(atoms: Atoms, indices: list) -> np.array:
    """
    Fit a plane to a set of points from an ASE Atoms object and return the normal vector.

    Parameters:
    atoms (ase.Atoms): The ASE Atoms object containing the points.
    indices (list): The list of indices of the atoms to use.

    Returns:
    numpy.ndarray: The normal vector to the fitted plane.
    """
    # Extract the points from the Atoms object based on the indices
    points = atoms.positions[indices]

    # Center the points
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid

    # Perform SVD on the centered points
    _, _, vh = svd(centered_points)

    # The normal vector is the last row of vh
    vector = vh[-1, :]

    return vector


def adjust_positions(positions: np.array, idx1: int, idx2: int, overlap: float):
    """Adjust positions of two atoms to avoid overlap."""
    pos1 = positions[idx1]
    pos2 = positions[idx2]
    direction = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
    positions[idx1] -= direction * overlap
    return positions


def minimum_image_distance(pos1: np.array, pos2: np.array, cell) -> float:
    """Calculate the minimum image distance between two positions."""
    delta = pos2 - pos1
    delta -= np.round(delta / cell) * cell
    return np.linalg.norm(delta)


def check_and_adjust(
    atoms: Atoms, positions: list, indices: list, tolerance: float
) -> Atoms:
    """Check for collisions and adjust positions."""
    check = False
    cell = atoms.get_cell().diagonal()
    for idx1 in indices:
        for idx2 in range(len(atoms)):
            if idx1 == idx2:
                continue
            pos1 = positions[idx1]
            pos2 = positions[idx2]
            distance = minimum_image_distance(pos1, pos2, cell)
            radius1 = covalent_radii[atoms.numbers[idx1]]
            radius2 = covalent_radii[atoms.numbers[idx2]]
            if distance < (radius1 + radius2) * (1 - tolerance):
                overlap = (radius1 + radius2) - distance
                positions = adjust_positions(positions, idx1, idx2, overlap)
                check = True
    return positions, check


def rotate_and_adjust_dynamically(
    atoms,
    indices,
    normal_vector,
    angle_degrees,
    step_degrees=1,
    tolerance=0.2,
    steps=100,
):
    """
    Rotate the specified atoms along the normal vector by a given angle,
    dynamically adjusting positions to avoid collisions.

    Parameters:
    atoms (ase.Atoms): The ASE Atoms object.
    indices (list): The list of indices of the atoms to rotate.
    normal_vector (numpy.ndarray): The normal vector to rotate along.
    angle_degrees (float): The rotation angle in degrees.
    step_degrees (float): The angle increment in degrees for each step.
    tolerance (float): The tolerance for collision adjustment.

    Returns:
    ase.Atoms: The modified ASE Atoms object with rotated and adjusted atoms.
    """

    # Convert angle from degrees to radians and calculate incremental angle
    incremental_angle = np.radians(step_degrees)

    # Normalize the normal vector
    normal_vector = normal_vector / np.linalg.norm(normal_vector)

    # Rotation matrix using incremental Rodrigues' rotation formula
    k = np.array(
        [
            [0, -normal_vector[2], normal_vector[1]],
            [normal_vector[2], 0, -normal_vector[0]],
            [-normal_vector[1], normal_vector[0], 0],
        ]
    )

    i = np.eye(3)
    incremental_r = (
        i
        + np.sin(incremental_angle) * k
        + (1 - np.cos(incremental_angle)) * np.dot(k, k)
    )

    # Calculate the center of mass of the specified indices
    subset_atoms = atoms[indices]
    center_of_mass = subset_atoms.get_center_of_mass()

    # Copy the positions to avoid modifying the original atom positions
    positions = atoms.positions.copy()

    # Rotate the specified atoms incrementally and adjust for collisions dynamically
    steps = int(angle_degrees / step_degrees)
    for _ in range(steps):
        for index in indices:
            positions[index] = (
                np.dot(incremental_r, positions[index] - center_of_mass)
                + center_of_mass
            )
        positions, check = check_and_adjust(atoms, positions, indices, tolerance)

    for i in range(steps):
        if check:
            positions, check = check_and_adjust(atoms, positions, indices, tolerance)
        else:
            break

    # Update the positions of the specified indices
    atoms.positions[indices] = positions[indices]

    return atoms


def translate_and_adjust(atoms, indices, translation_vector, tolerance=0.2, steps=100):
    """
    Translate the specified atoms by a given distance,
    dynamically adjusting positions to avoid collisions.

    Parameters:
    atoms (ase.Atoms): The ASE Atoms object.
    indices (list): The list of indices of the atoms to translate.
    translation_vector (numpy.ndarray): The translation vector.
    tolerance (float): The tolerance for collision adjustment.

    Returns:
    ase.Atoms: The modified ASE Atoms object with translated and adjusted atoms.
    """

    # Copy the positions to avoid modifying the original atom positions
    positions = atoms.positions.copy()

    # Translate the specified atoms
    for index in indices:
        positions[index] += translation_vector

    # Adjust for collisions dynamically
    positions, check = check_and_adjust(atoms, positions, indices, tolerance)

    for i in range(steps):
        if check:
            positions, check = check_and_adjust(atoms, positions, indices, tolerance)
        else:
            break

    # Update the positions of the specified indices
    atoms.positions[indices] = positions[indices]

    return atoms
