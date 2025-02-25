from ase import Atoms
from ase.build import fcc111
from ase.calculators.emt import EMT

adsH = Atoms("H", positions=[(0, 0, 0)])
atoms = fcc111("Pt", size=(3, 3, 4), vacuum=10.0)
atoms.calc = EMT()  # Add a calculator to the atoms object for geometric optimization

from gg.gcbh import Gcbh

G = Gcbh(atoms, config_file="input.yaml")

# Attach Modfiers to gcbh
max_coord = {"Pt": 12, "H": 2}

from gg.predefined_sites import SurfaceSites

ss = SurfaceSites(max_coord, max_bond_ratio=1.2)

from gg.modifiers import Add, Remove, ModifierAdder

addH = Add(
    ss,
    adsH,
    surf_coord=[1, 2, 3],
    ads_id=["H"],
    surf_sym=["Pt"],
    print_movie=False,
    weight=1.0,
)
remH = Remove(ss, adsH, max_bond_ratio=1.2, print_movie=False, weight=1.0)
swapH = ModifierAdder([addH, remH], weight=1)

G.add_modifier(addH, "Add_H")
G.add_modifier(remH, "Remove_H")
G.add_modifier(swapH, "Swap_H")

G.run(steps=25)
