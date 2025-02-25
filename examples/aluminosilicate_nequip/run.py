#! /usr/bin/env python
"""Example File to use nequip calculator"""

from ase.io import read
from ase import Atoms
from gg.modifiers import Add,Remove,Swap,ModifierAdder
from gg.gcbh import Gcbh
from gg.predefined_sites import SurfaceSites
from nequip.ase import NequIPCalculator

#Read Input POSCAR
atoms = read("POSCAR")

#Define Calculator
CALC_PATH = './deployed_model.pth'
atom_dict = {'Si': 'Si', 'H':'H', 'O':'O', 'Al':'Al'}
calc = NequIPCalculator.from_deployed_model(CALC_PATH,species_to_type_name=atom_dict,device='cuda')
atoms.calc = calc

#Define Adsorbates to be used
adsH = Atoms("H", positions = [(0,0,0)])
adsOH = Atoms("OH", positions = [(0,0,0),(0.1,0.1,1.0)])

#Define Surface Sites
max_coord = {"Al": 6, "Si": 5, "O": 4, "H": 2}
ss = SurfaceSites(max_coord, max_bond_ratio=1.2)

#Build Add H2O
addH = Add(ss, adsH, 1, ads_dist=1.0, surf_sym = ["O"], print_movie=False, weight = 1.0)
addOH = Add(ss, adsOH, 1, ads_dist=1.8, surf_sym = ["Al","Si"],print_movie=False, weight = 1.0)
addH2O = ModifierAdder(1.0,[addH,addOH])

#Build Remove H2O
remH = Remove(ss, adsH, max_bond = 2, weight = 1.0)
remOH = Remove( ss, adsOH, max_bond = 2, weight = 1.0)
remH2O = ModifierAdder([remH,remOH])

#Build Swap Al and Si Atoms
swap_al_si = Swap(1.0,ss,["Al","Si"])

#Build Swap H2O
swap_H2O = ModifierAdder(1.0,[addH2O,remH2O])

#Initialize a GCBH calculator
G = Gcbh(atoms,config_file='input.yaml')

#Attach modifiers to the gcbh calculator
G.add_modifier(addH2O,'Add_H2O')
G.add_modifier(remH2O,'Remove_H2O')
G.add_modifier(swap_al_si,'Swap_Al_Si')
G.add_modifier(swap_H2O,'Swap_H2O')

#Run GCBH
G.run(steps = 1000)
