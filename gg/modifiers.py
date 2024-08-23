import networkx as nx
import numpy as np
import pandas as pd
from ase.neighborlist import NeighborList,natural_cutoffs
from ase.data import atomic_numbers
from ase.io import read as read_atoms
from numpy.linalg import norm
from ase import Atoms
from utils import *

class parent_modifier:
    def __init__(self,name,atoms,weight):
        self.name = name
        if isinstance(atoms, str):
            self.og_atomss = read_atoms(atoms)
        elif isinstance(atoms, Atoms):
            self.og_atoms = atoms.copy()
        else:
            print("Please provide proper atoms file")
        self._atoms = self.og_atoms.copy()
        self.og_weight = weight
        self.weight = weight

    @property
    def atoms(self):
        return self._atoms 

    @atoms.setter 
    def atoms(self, atoms):
        self._atoms = atoms.copy()
    
            
class rattle(parent_modifier):
    def __init__(self, name, atoms, weight, stdev=0.001, contact_error=0.2):
        self.stdev = stdev
        self.contact_error = contact_error
        super().__init__(name, atoms, weight)
        
    def get_modified_atoms(self):
        self.atoms.rattle(stdev=self.stdev)
        if check_contact(self.atoms,error=self.contact_error):
            print("Atoms touching")
        return self.atoms

class add(parent_modifier):
    def __init__(self, name, atoms, weight, surface_sites, ads, ads_coord, ad_dist=1.8, movie=False):
        super().__init__(name, atoms, weight)
        self.ss = surface_sites
        self.ads = ads
        self.ads_coord = ads_coord
        self.ad_dist = ad_dist
        self.movie = movie
        
    def get_modified_atoms(self):
        self.df, self.G = self.ss.get_surface_sites(self.atoms)
        movie = generate_sites(self.atoms, self.ads, self.G, self.df["ind"], self.ads_coord, ad_dist = self.ad_dist, contact_error = self.ss.contact_error)

        if self.movie:
            return movie
        else:
            return movie[0]

class remove(parent_modifier):
    def __init__(self, name, atoms, weight, surface_sites):
        super().__init__(name, atoms, weight)
        self.ss = surface_sites
        
    def get_modified_atoms(self):
        self.df, self.G = self.ss.get_surface_sites(self.atoms)
        ind_to_remove = int(self.df["ind"].to_list()[0])
        del self.atoms[ind_to_remove]
        return self.atoms