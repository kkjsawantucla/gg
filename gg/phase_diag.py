"""Plot Phase Diagram"""

import os
from functools import cmp_to_key
import yaml
import numpy as np
from scipy.spatial import HalfspaceIntersection
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read
from gg.utils import (
    extract_lowest_energy_from_oszicar,
    extract_lowest_energy_from_outlog,
)
from gg.reference import get_ref_coeff


class Phasediagramentry:
    """_summary_"""

    def __init__(self, enid, energy, n1, n2):
        """
        Args:
            entry (composition): Anentry object
            ref_entry (composition): Anentry object
        """
        self.energy = energy
        self.enid = enid
        self.n1 = n1
        self.n2 = n2

    def __repr__(self):
        return (
            f"Entry:{self.enid} with energy={self.energy:.4f}, x={self.n1}, y={self.n2}"
        )


# From Pymatgen
def get_phase_domains(entries, limits=None):
    """
    Args:
        entries (List of phasediagramentry class):
        limits (list, optional): _description_. Defaults to None.
    Returns:
        dict : domain vertices
    """
    if limits is None:
        limits = [[-4, 0], [-4, 0]]
    # Get hyperplanes
    # +nO+nH+E-E0 <= 0
    hyperplanes = [
        np.array([entry.n1, entry.n2, 1, -entry.energy]) for entry in entries
    ]
    hyperplanes = np.array(hyperplanes)

    max_contribs = np.max(np.abs(hyperplanes), axis=0)
    g_max = np.dot(
        -max_contribs, [8, 8, 0, 1]
    )  # Use to define minimum point which is the space away from minimum hypaerplane
    border_hyperplanes = [
        [-1, 0, 0, limits[0][0]],
        [1, 0, 0, -limits[0][1]],
        [0, -1, 0, limits[1][0]],
        [0, 1, 0, -limits[1][1]],
        [0, 0, -1, 2 * g_max],
    ]  # Define border planes
    hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])
    interior_point = np.average(limits, axis=1).tolist() + [
        g_max
    ]  # Define interior point
    hs_int = HalfspaceIntersection(
        hs_hyperplanes, np.array(interior_point)
    )  # get intersection points

    domains = {entry: [] for entry in entries}
    for intersection, facet in zip(hs_int.intersections, hs_int.dual_facets):
        for v in facet:
            if v < len(entries):
                this_entry = entries[v]
                domains[this_entry].append(intersection)

    # Remove entries with no pourbaix region
    domains = {k: v for k, v in domains.items() if v}
    domain_vertices = {}

    for entry, points in domains.items():
        points = np.array(points)[:, :2]  # drop the energy axis
        # Initial sort to ensure consistency
        points = points[np.lexsort(np.transpose(points))]
        center = np.average(points, axis=0)
        points_centered = points - center

        # Sort points by cross product of centered points,
        # isn't strictly necessary but useful for plotting tools
        points_centered = sorted(
            points_centered, key=cmp_to_key(lambda x, y: x[0] * y[1] - x[1] * y[0])
        )
        points = points_centered + center
        domain_vertices[entry] = points

    return domain_vertices


# Plotting Function
def phase_diagram_plot(stable_domain_vertices, limits, xlabel="1", ylabel="2"):
    """
    Args:
        stable_domain_vertices (dict): 
        limits (list): 
        xlabel (str, optional): . Defaults to "1".
        ylabel (str, optional): . Defaults to "2".
    """
    x_min = limits[0][0]
    x_max = limits[0][1]
    for entry, vertices in stable_domain_vertices.items():
        center = np.average(vertices, axis=0)
        x, y = np.transpose(np.vstack([vertices, vertices[0]]))
        plt.plot(x, y, "k-")
        plt.annotate(
            entry.enid, center, ha="center", va="center", fontsize=20, color="b"
        ).draggable()

    plt.xlim(x_min + 0.01, x_max - 0.01)
    plt.ylim(limits[1][0] + 0.01, limits[1][1] - 0.01)
    plt.xlabel(
        f"{xlabel} Chemical Potential \u03bc(O) [eV]",
        fontsize=24,

        labelpad=20,
    )
    plt.ylabel(
        f"{ylabel} Chemical Potential \u03bc(H) [eV]",
        fontsize=24,
        labelpad=20,
    )
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.tick_params(direction="in", length=10, width=2, top=True, right=True, pad=10)
    plt.tick_params(
        which="minor", direction="in", length=6, width=1, top=True, right=True
    )
    plt.axhline(y=-0.15, color="r", linestyle="--")
    fig = plt.gcf()
    fig.set_size_inches(15, 8)
    plt.savefig(f"{xlabel}_{ylabel}_fig")
    return


def read_chemical_potential(path):
    """
    Args:
        path (str): 

    Returns:
        dict: 
    """
    with open(path, "r", encoding="utf-8") as f:
        input_config = yaml.safe_load(f)
    return input_config["chemical_potential"]


def get_ref_potential(mu, atoms: Atoms, n1, n2):
    """
    Args:
        atoms (ase.Atoms):
    Returns:
        float: total ref value to substract
    """
    if n1 not in mu or n2 not in mu:
        raise RuntimeError(f"{n1} or {n2} not in mu")
    formula = atoms.get_chemical_formula()
    ref_sum = 0
    ref_coeff = get_ref_coeff(mu, formula)
    for key, value in mu.items():
        if key not in [n1, n2]:
            ref_sum += ref_coeff[key] * value
        elif key == n1:
            n1_slope = ref_coeff[key]
        elif key == n2:
            n2_slope = ref_coeff[key]
    return ref_sum, n1_slope, n2_slope


def get_entries_from_folders(
    n1, n2, base_folder="./", mu_path="./input.yaml", file_type=["OSZICAR", "CONTCAR"]
):
    mu = read_chemical_potential(mu_path)
    entries = []
    for root, _, files in os.walk(base_folder):
        if file_type[0] in files and file_type[1] in files:
            contcar_path = os.path.join(root, file_type[1])
            en_path = os.path.join(root, file_type[0])
            if file_type[0] == "OSZICAR":
                energy = extract_lowest_energy_from_oszicar(en_path)
            elif file_type[0] == "out.log":
                energy = extract_lowest_energy_from_outlog(en_path)
            else:
                raise RuntimeError("Unsupported energy file type.")

            if energy is None:
                continue

            atoms = read(contcar_path, format="vasp")
            ref_sum, n1_slope, n2_slope = get_ref_potential(mu, atoms, n1, n2)
            final_energy = energy - ref_sum
            entry_id = os.path.basename(root.replace("/", "_")[:-1])
            print(f"Adding entry: {entry_id} {n1}={n1_slope}, {n2}={n2_slope}")
            entry = Phasediagramentry(
                enid=entry_id, energy=final_energy, n1=n1_slope, n2=n2_slope
            )
            entries.append(entry)
    ref_entry = Phasediagramentry(enid="Reference", energy=0, n1=0, n2=0)
    entries.append(ref_entry)
    return entries


def plot_phase_diagram_from_run(
    n1,
    n2,
    limits=[[-3.5, -1.5], [-1, 0.5]],
    base_folder="./",
    mu_path="./input.yaml",
    file_type=["OSZICAR", "CONTCAR"],
):
    print(f"Generating entries for plotting from {base_folder} folder")
    entries = get_entries_from_folders(
        n1, n2, base_folder=base_folder, mu_path=mu_path, file_type=file_type
    )
    stable_vertices = get_phase_domains(entries, limits=limits)
    print(f"Plotting with xlabel: {n1} and ylabel: {n2}")
    phase_diagram_plot(stable_vertices, limits=limits, xlabel = n1, ylabel = n2)

    return
