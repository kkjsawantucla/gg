#! /usr/bin/env python
"""Example File to use nequip calculator"""

from gg.modifiers import Add,Remove,Swap,ModifierAdder
from gg.gcbh import Gcbh
from gg.sites import SurfaceSites
from nequip.ase import NequIPCalculator
from ase.io import read
from ase import Atoms

calc_path = './deployed_model.pth'
atom_dict = {'Si': 'Si', 'H':'H', 'O':'O', 'Al':'Al'}
calc = NequIPCalculator.from_deployed_model(calc_path,species_to_type_name=atom_dict,device='cuda')
atoms = read("POSCAR")
atoms.calc = calc
print(f"Successfully read atoms {atoms.get_chemical_formula()} and calculator {atoms.calc}")

adsH = Atoms("H", positions = [(0,0,0)])
adsOH = Atoms("OH", positions = [(0,0,0),(0.1,0.1,1.0)])
max_coord = {"Al": 6, "Si": 5, "O": 4, "H": 2}

ss = SurfaceSites(max_coord, max_bond_ratio=1.2)

#Build AddH2O
addH = Add(ss, adsH, 1, ads_dist=1.0, surf_sym = ["O"], print_movie=False, weight = 1.0)
addOH = Add(ss, adsOH, 1, ads_dist=1.8, surf_sym = ["Al","Si"],print_movie=False, weight = 1.0)
addH2O = ModifierAdder(1.0,[addH,addOH])

#Remove AddH2O
remH = Remove(ss, adsH, max_bond = 2, weight = 1.0)
remOH = Remove( ss, adsOH, max_bond = 2, weight = 1.0)
remH2O = ModifierAdder([remH,remOH])

#Swap Al,Si
swap_al_si = Swap(1.0,ss,["Al","Si"])

#Swap H2O
swap_H2O = ModifierAdder(1.0,[addH2O,remH2O])

G = Gcbh(atoms,config_file='input.yaml')

G.add_modifier(addH2O,'Add_H2O')
G.add_modifier(remH2O,'Remove_H2O')
G.add_modifier(swap_al_si,'Swap_Al_Si')
G.add_modifier(swap_H2O,'Swap_H2O')

G.run(steps = 1000)
