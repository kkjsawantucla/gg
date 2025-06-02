"""Plot Phase Diagram with stoichiometry-based entry filtering"""

import os
import json
from dataclasses import dataclass, field, asdict
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
    get_area,
)
from gg.reference import get_ref_coeff


@dataclass(frozen=True)
class PhaseEntry:
    """Entry for phase diagram, with stoichiometry tracking"""

    enid: str  # Identifier for the entry
    energy: float  # Formation energy
    n1: float  # Slope coefficient for species 1
    n2: float  # Slope coefficient for species 2
    stoich: str = field(default=None, compare=False)  # Chemical formula string


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


def phase_diagram_plot(
    stable_domain_vertices,
    limits,
    xlabel="1",
    ylabel="2",
    annotate=False,
    mu=None,
    number_labels=True,  # NEW ─ turn numeric labelling on/off
):
    """
    Args:
        stable_domain_vertices (dict): PhaseEntry → list of vertices
        limits (list): [[xmin,xmax],[ymin,ymax]]
        xlabel (str): x-axis species
        ylabel (str): y-axis species
        annotate (bool): write labels inside the regions
        mu (dict|None): draw reference lines if given
        number_labels (bool): when annotate is True, write numbers
                              instead of full enid and print a legend
    """
    # --- plot each domain ---------------------------------------------------
    mapping = []
    for idx, (entry, vertices) in enumerate(stable_domain_vertices.items(), start=1):
        center = np.average(vertices, axis=0)
        x, y = np.transpose(np.vstack([vertices, vertices[0]]))
        plt.plot(x, y, "k-")

        if annotate:
            label = str(idx) if number_labels else entry.enid
            plt.annotate(
                label, center, ha="center", va="center", fontsize=18, color="b"
            ).draggable()
        if number_labels:
            mapping.append(f"{idx}: {entry.enid}")

    # --- axes limits / labels ----------------------------------------------
    plt.xlim(limits[0][0] + 0.01, limits[0][1] - 0.01)
    plt.ylim(limits[1][0] + 0.01, limits[1][1] - 0.01)
    plt.xlabel(
        f"{xlabel} chemical potential μ({xlabel}) [eV]", fontsize=24, labelpad=20
    )
    plt.ylabel(
        f"{ylabel} chemical potential μ({ylabel}) [eV]", fontsize=24, labelpad=20
    )
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.tick_params(direction="in", length=10, width=2, top=True, right=True, pad=10)
    plt.tick_params(
        which="minor", direction="in", length=6, width=1, top=True, right=True
    )

    # optional reference lines
    if mu:
        plt.axhline(y=mu[str(ylabel)], color='r', linestyle="--")
        plt.axvline(x=mu[str(xlabel)], color='r', linestyle="--")

    # --- print legend of numbers vs enid ------------------------------------
    if annotate and number_labels and mapping:
        legend_text = "\n".join(mapping)
        # place it just below the axes; tweak y-offset as needed
        plt.gcf().text(0.01, -0.08, legend_text, ha="left", va="top", fontsize=12)

    # --- figure size / save -------------------------------------------------
    plt.gcf().set_size_inches(15, 8)
    plt.savefig(f"{xlabel}_{ylabel}_phase_diag")


def read_chemical_potential(path):
    """
    Args:
        path (str): Path to YAML file containing chemical potentials

    Returns:
        dict: Mapping species -> chemical potential value
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
    n1,
    n2,
    base_folders=["./"],
    mu_path="./input.yaml",
    file_type=["OSZICAR", "CONTCAR"],
    reference=None,
):
    """
    Walk through subdirectories to collect entries, keeping only one per stoichiometry
    (lowest energy) plus a reference.

    Args:
        n1, n2 (str): Labels for chemical potentials
        base_folder (str): Root directory for VASP runs
        mu_path (str): Path to YAML file with chemical potentials
        file_type (list): [energy file, structure file]

    Returns:
        List[Phasediagramentry]: Filtered entries
    """
    mu = read_chemical_potential(mu_path)
    entries = []
    for base_folder in base_folders:
        for root, _, files in os.walk(base_folder):
            if file_type[0] in files and file_type[1] in files:
                contcar_path = os.path.join(root, file_type[1])
                en_path = os.path.join(root, file_type[0])

                # extract energy
                if file_type[0] == "OSZICAR":
                    energy = extract_lowest_energy_from_oszicar(en_path)
                elif file_type[0].endswith(".log"):
                    energy = extract_lowest_energy_from_outlog(en_path)
                else:
                    raise RuntimeError("Unsupported energy file type.")

                if energy is None:
                    continue

                atoms = read(contcar_path, format="vasp")
                area = get_area(atoms)
                ref_sum, n1_slope, n2_slope = get_ref_potential(mu, atoms, n1, n2)
                final_energy = (energy - ref_sum) / area
                n1_slope = n1_slope / area
                n2_slope = n2_slope / area
                stoich_formula = atoms.get_chemical_formula()
                entry_id = (
                    os.path.basename(root).replace("/", "_") + "_" + str(stoich_formula)
                )
                new_entry = PhaseEntry(
                    enid=entry_id,
                    energy=final_energy,
                    n1=n1_slope,
                    n2=n2_slope,
                    stoich=stoich_formula,
                )

                # check for existing same stoichiometry
                replaced = False
                for idx, existing in enumerate(entries):
                    if existing.stoich == new_entry.stoich:
                        # keep lower energy
                        if new_entry.energy < existing.energy:
                            print(
                                f"Replacing {existing.enid} with {new_entry.enid} for stoich={new_entry.stoich}"
                            )
                            entries[idx] = new_entry
                        else:
                            print(
                                f"Skipping {new_entry.enid} (higher energy than existing {existing.enid})"
                            )
                        replaced = True
                        break

                if not replaced:
                    entries.append(new_entry)

    # always include a zero-energy reference
    if reference:
        ref_entry = PhaseEntry(enid="Reference", energy=reference, n1=0.0, n2=0.0)
        entries.append(ref_entry)
    print(f"Number of unique entries found: {len(entries)}")
    return entries


def save_entries(entries, filename: str):
    """
    Save a list of PhaseEntry objects to a JSON file.
    """
    data = [asdict(entry) for entry in entries]
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def load_entries(filename: str):
    """
    Load a list of PhaseEntry objects from a JSON file.
    """
    with open(filename, "r", encoding="utf-8") as f:
        data = json.load(f)
    return [PhaseEntry(**d) for d in data]


def plot_phase_diagram_from_run(
    n1,
    n2,
    limits=[[-3.5, -1.5], [-1, 0.5]],
    base_folders=["./"],
    mu_path="./input.yaml",
    file_type=["OSZICAR", "CONTCAR"],
    read_from_file=False,
    annotate=True,
    number_labels=True,
):
    print(f"Generating entries for plotting from {base_folders} folder")
    mu = read_chemical_potential(mu_path)
    print(f"x-axis is {n1}: {mu[str(n1)]}eV")
    print(f"y-axis is {n2}: {mu[str(n2)]}eV")
    if read_from_file and os.path.isfile(read_from_file):
        print(f"Reading from {read_from_file}")
        entries = load_entries(read_from_file)
    else:
        entries = get_entries_from_folders(
            n1, n2, base_folders=base_folders, mu_path=mu_path, file_type=file_type
        )
        json_path = os.path.join("./", f"entries_{n1}_{n2}.json")
        save_entries(entries, json_path)
    print("Building 2D Hull")
    stable_vertices = get_phase_domains(entries, limits=limits)
    print(f"Plotting with xlabel:{n1} and ylabel:{n2}")
    phase_diagram_plot(
        stable_vertices,
        limits=limits,
        xlabel=n1,
        ylabel=n2,
        annotate=annotate,
        mu=mu,
        number_labels=number_labels,
    )

    return
