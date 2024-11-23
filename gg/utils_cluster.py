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


def center_and_rotate(
    atoms: Atoms, indices: list, normal_vector: np.array, angle_degrees: float
) -> Atoms:
    """
    Center the atoms around the center of mass of the specified indices and then rotate them along the normal vector.

    Parameters:
    atoms (ase.Atoms): The ASE Atoms object.
    indices (list): The list of indices of the atoms to use and rotate.
    normal_vector (numpy.ndarray): The normal vector to rotate along.
    angle_degrees (float): The rotation angle in degrees.

    Returns:
    ase.Atoms: The modified ASE Atoms object with atoms centered, rotated, and wrapped within the unit cell.
    """
    # Convert angle from degrees to radians
    angle_radians = np.radians(angle_degrees)

    subset_atoms = atoms[indices]
    center_of_mass = subset_atoms.get_center_of_mass()

    # Normalize the normal vector
    normal_vector = normal_vector / np.linalg.norm(normal_vector)

    # Rotation matrix using Rodrigues' rotation formula
    k = np.array(
        [
            [0, -normal_vector[2], normal_vector[1]],
            [normal_vector[2], 0, -normal_vector[0]],
            [-normal_vector[1], normal_vector[0], 0],
        ]
    )

    i = np.eye(3)
    r = i + np.sin(angle_radians) * k + (1 - np.cos(angle_radians)) * np.dot(k, k)

    # Rotate the specified atoms
    for index in indices:
        atoms.positions[index] = np.dot(r, atoms.positions[index] - center_of_mass) + center_of_mass

    return atoms

def rotate_and_adjust_dynamically(atoms, indices, normal_vector, angle_degrees, steps=50, tolerance = 0.2):
    """
    Rotate the specified atoms along the normal vector by a given angle,
    dynamically adjusting positions to avoid collisions.

    Parameters:
    atoms (ase.Atoms): The ASE Atoms object.
    indices (list): The list of indices of the atoms to rotate.
    normal_vector (numpy.ndarray): The normal vector to rotate along.
    angle_degrees (float): The rotation angle in degrees.
    steps (int): The number of steps to break the rotation into for dynamic adjustment.

    Returns:
    ase.Atoms: The modified ASE Atoms object with rotated and adjusted atoms.
    """
    def adjust_positions(positions, idx1, idx2, overlap):
        """Adjust positions of two atoms to avoid overlap."""
        pos1 = positions[idx1]
        pos2 = positions[idx2]
        direction = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)
        positions[idx1] -= direction * (overlap/2)
        positions[idx2] += direction * (overlap/2)

    def check_and_adjust(atoms, positions, indices):
        """Check for collisions and adjust positions."""
        for idx1 in indices:
            for idx2 in range(len(atoms)):
                if idx1 == idx2:
                    continue
                pos1 = positions[idx1]
                pos2 = positions[idx2]
                distance = np.linalg.norm(pos1 - pos2)
                radius1 = covalent_radii[atoms.numbers[idx1]]
                radius2 = covalent_radii[atoms.numbers[idx2]]
                if distance < (radius1 + radius2)*(1-tolerance):
                    overlap = (radius1 + radius2) - distance
                    adjust_positions(positions, idx1, idx2, overlap)

    # Convert angle from degrees to radians and calculate incremental angle
    total_angle_radians = np.radians(angle_degrees)
    incremental_angle = total_angle_radians / steps

    # Normalize the normal vector
    normal_vector = normal_vector / np.linalg.norm(normal_vector)

    # Rotation matrix using incremental Rodrigues' rotation formula
    K = np.array([
        [0, -normal_vector[2], normal_vector[1]],
        [normal_vector[2], 0, -normal_vector[0]],
        [-normal_vector[1], normal_vector[0], 0]
    ])

    i = np.eye(3)
    incremental_R = i + np.sin(incremental_angle) * K + (1 - np.cos(incremental_angle)) * np.dot(K, K)

    # Calculate the center of mass of the specified indices
    subset_atoms = atoms[indices]
    center_of_mass = subset_atoms.get_center_of_mass()

    # Copy the positions to avoid modifying the original atom positions
    positions = atoms.positions.copy()

    # Rotate the specified atoms incrementally and adjust for collisions dynamically
    for _ in range(steps):
        for index in indices:
            positions[index] = np.dot(incremental_R, positions[index] - center_of_mass) + center_of_mass
        check_and_adjust(atoms, positions, indices)

    # Update the positions of the specified indices
    atoms.positions[indices] = positions[indices]

    return atoms
