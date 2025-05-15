import os
import argparse
from ase.io import read
from ase import Atoms
from gg.modifiers import Add, ModifierAdder
from gg.sites import RuleSites, get_com_sites, get_surface_sites_by_coordination
from gg.utils_graph import get_unique_atoms

parser = argparse.ArgumentParser(description="Command line utilities")
parser.add_argument(
    "--surface",
    "-surface",
    "-s",
    type=str,
    nargs="+",
    help="Path to surface atoms objects to adsorb",
)

parser.add_argument(
    "--dir",
    "-dir",
    "-d",
    type=str,
    default="above",
    help="Direction of the surface",
)

args = parser.parse_args()


adsSi = Atoms("Si", positions=[(0, 0, 0)])
adsO = Atoms("O", positions=[(0, 0, 0)])
max_coord = {"Al": 6, "Si": 4, "O": 3, "H": 1}

ss = RuleSites(
    index_parsers=[
        lambda atoms: get_com_sites(atoms, fraction=0.50, direction=str(args.dir)),
        lambda atoms: get_surface_sites_by_coordination(
            atoms, max_coord, max_bond=2, max_bond_ratio=0.0
        ),
    ],
    combine_rules="intersection",
    max_bond=2,
    max_bond_ratio=0.0,
    contact_error=0.2,
)

addSi = Add(
    ss,
    adsSi,
    1,
    ads_id="Si",
    surf_sym=["O"],
    print_movie=True,
    unique=False,
    unique_method="Si",
)
addO = Add(
    ss,
    adsO,
    1,
    ads_id="O",
    surf_sym=["Al", "Si"],
    print_movie=True,
    unique=False,
    unique_method="Si",
)
addSiO2_2 = ModifierAdder(
    [addO, addSi, addO],
    print_movie=True,
    unique=True,
    max_bond=2,
    max_bond_ratio=0.0,
    unique_method="Si",
)
addSiO2 = ModifierAdder(
    [addO, addO, addSi],
    print_movie=True,
    unique=True,
    max_bond=2,
    max_bond_ratio=0.0,
    unique_method="Si",
)

print("Starting Loop")
traj = []
for surface_name in args.surface:
    atoms = read(surface_name)
    print(f"read atoms at {surface_name}")
    sites = ss.get_sites(atoms)
    print(f"Found {len(sites)} sites")

    # AddH2O
    add_atoms = addSiO2.get_modified_atoms(atoms)
    add_atoms += addSiO2_2.get_modified_atoms(atoms)
    print(f"Found {len(add_atoms)} atoms")
    traj = traj + add_atoms

print("Generating Unique Sites")
add_atoms = get_unique_atoms(traj, max_bond=2, max_bond_ratio=0, unique_method="Si")
print(f"Found {len(add_atoms)} unique atoms")

i = 1
for atoms in add_atoms:
    os.makedirs(f"w_{i}")
    atoms.write(f"w_{i}/POSCAR", format="vasp")
    i += 1
