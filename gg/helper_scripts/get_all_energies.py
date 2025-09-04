#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from ase.io import read

import sys
sys.path.insert(0, "/jet/home/ksawant/apps/gg")
from gg.phase_diag import plot_phase_diagram_from_run

# --- from your file (phase_diag.py) ---
from gg.phase_diag import (
    read_chemical_potential,
    get_vib_correction,
    get_ref_coeff,
)
# --- from gg.utils used in your file ---
from gg.utils import (
    extract_lowest_energy_from_oszicar,
    extract_lowest_energy_from_outlog,
    get_area,
)

def get_ref_potential(mu, atoms):
    """
    Args:
        atoms (ase.Atoms):
    Returns:
        float: total ref value to substract
    """
    formula = atoms.get_chemical_formula()
    ref_sum = 0
    ref_coeff = get_ref_coeff(mu, formula)
    for key, value in mu.items():
        ref_sum += ref_coeff[key] * value
    return ref_sum

def collect_energies(
    base_folders,
    mu_path: str = "./input.yaml",
    energy_file: str = "out.log",      # or "out.log"
    structure_file: str = "CONTCAR",   # or "CONTCAR_mlp", etc.
    vib_corrections: dict | None = None,
):
    """
    Walk through base_folders, read energy + stoichiometry, compute formation energy per area (fe).
    Returns a pandas DataFrame with: path, energy, fe, stoichiometry
    """
    mu = read_chemical_potential(mu_path)
    rows = []

    for base in base_folders:
        for root, _, files in os.walk(base):
            if energy_file in files and structure_file in files:
                en_path = os.path.join(root, energy_file)
                struct_path = os.path.join(root, structure_file)

                # raw electronic energy
                if energy_file == "OSZICAR":
                    energy = extract_lowest_energy_from_oszicar(en_path)
                elif energy_file.endswith(".log"):
                    energy = extract_lowest_energy_from_outlog(en_path)
                else:
                    # Unsupported energy type; skip
                    continue
                if energy is None:
                    continue

                # structure + area + stoichiometry
                atoms = read(struct_path, format="vasp")
                area = get_area(atoms)
                stoich = atoms.get_chemical_formula()

                # reference + vib corrections (same logic as in your PD code)
                ref_sum = get_ref_potential(mu, atoms)
                vib_corr = get_vib_correction(atoms, vib_corrections)

                # fe: formation energy per surface area, consistent with your plotting workflow
                fe = (energy - ref_sum + vib_corr) / area

                rows.append(
                    {
                        "path": root,
                        "energy": energy,
                        "fe": fe,
                        "stoichiometry": stoich,
                    }
                )
                print(root,fe)

    df = pd.DataFrame(rows).sort_values(["fe"]).reset_index(drop=True)
    return df


def main():
    p = argparse.ArgumentParser(description="Collect energies to CSV: path, energy, fe, stoichiometry")
    p.add_argument("--folders", nargs="+", default=["./"], help="Base folders to search")
    p.add_argument("--mu", default="./input.yaml", help="YAML with chemical_potential")
    p.add_argument("--energy-file", default="opt.log", help="Energy file name: OSZICAR or out.log")
    p.add_argument("--structure-file", default="CONTCAR", help="Structure file name (e.g., CONTCAR)")
    p.add_argument("--csv", default="energies.csv", help="Output CSV path")
    args = p.parse_args()

    df = collect_energies(
        base_folders=args.folders,
        mu_path=args.mu,
        energy_file=args.energy_file,
        structure_file=args.structure_file,
        vib_corrections=None,  # plug your dict if you use it
    )
    df.to_csv(args.csv, index=False)
    print(f"Wrote {len(df)} rows to {args.csv}")


if __name__ == "__main__":
    main()

