"""Import Modules for basic functioning"""

import random
from ase.io import read as read_atoms
from ase import Atoms
from networkx.algorithms import isomorphism
from gg.utils import (
    check_contact,
    generate_sites,
    NoReasonableStructureFound,
    custom_copy,
    formula_to_graph,
)
from gg.reference import is_element, is_chemical_formula

__author__ = "Kaustubh Sawant"


class ParentModifier:
    """Parent bare bones modifier which serves as basis for other modifiers
    Args:
        weight (float): Modifier Weight
    """

    def __init__(self, weight):
        self.og_weight = weight
        self.weight = weight

    @property
    def atoms(self):
        """
        Returns:
            ase.Atoms
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        if isinstance(atoms, str):
            self._atoms = read_atoms(atoms)
        elif isinstance(atoms, Atoms):
            self._atoms = custom_copy(atoms)
        else:
            print("Please provide proper atoms file")

    def get_modified_atoms(self, atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class ModifierAdder(ParentModifier):
    """Add modifiers together to create a new modifier"""

    def __init__(self, weight, modifier_instances):
        super().__init__(weight)
        if isinstance(modifier_instances, list):
            self.modifier_instances = modifier_instances
        else:
            raise RuntimeError("modifier_instances isnt a list. Please provide a list")

    def get_modified_atoms(self, atoms):
        """
        Returns:
            ase.Atoms:
        """
        for instance in self.modifier_instances:
            self.atoms = instance.get_modified_atoms(atoms)
        return self.atoms


class Rattle(ParentModifier):
    """Modifier that rattles the atoms with some stdev"""

    def __init__(self, weight, stdev=0.001, contact_error=0.2):
        self.stdev = stdev
        self.contact_error = contact_error
        super().__init__(weight)

    def get_modified_atoms(self, atoms):
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        self.atoms.rattle(stdev=self.stdev)
        if check_contact(self.atoms, error=self.contact_error):
            raise NoReasonableStructureFound("Atoms Touching")
        else:
            return self.atoms


class Add(ParentModifier):
    """Modifier that adds an adsorbate at certain specific sites"""

    def __init__(
        self,
        weight,
        surface_sites,
        ads,
        ads_coord,
        ad_dist=1.8,
        movie=False,
    ):
        """
        Args:
            weight (float): _description_
            surface_sites (gg.SurfaceSites): a class which figures out surface sites
            ads (str) or (ase.Atoms): Adsorbate to add
            ads_coord (int): Adsorbate coordination number for adding
            ad_dist (float, optional): Distance of adsorbate from surface site.
            Defaults to 1.8.
            movie (bool, optional): return a movie o best sites or one random site.
            Defaults to False.
        """
        super().__init__(weight)
        self.ss = surface_sites
        if isinstance(ads, str):
            self.ads = read_atoms(ads)
        elif isinstance(ads, Atoms):
            self.ads = custom_copy(ads)
        self.ads_coord = ads_coord
        self.ad_dist = ad_dist
        self.print_moview = movie

    def get_modified_atoms(self, atoms):
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind, g = self.ss.get_surface_sites(self.atoms)
        movie = generate_sites(
            self.atoms,
            self.ads,
            g,
            df_ind,
            self.ads_coord,
            ad_dist=self.ad_dist,
            contact_error=self.ss.contact_error,
        )
        if not movie:
            raise NoReasonableStructureFound(
                "Movie was empty, most likely due to issues with atoms touching"
            )
        if self.print_moview:
            return movie
        else:
            return movie[0]


class Remove(
    ParentModifier,
):
    """Modifier that randomly removes an atom or molecule"""

    def __init__(self, weight, surface_sites, to_del, max_bond_ratio=1.2, max_bond=0):
        """
        Args:
            weight (str):

            surface_sites (gg.SurfaceSites): a class which figures out surface sites

            to_del (str) or (ase.Atoms): atoms to delete. If a string is provided,
            it tries to make a molecule out of it.

            max_bond_ratio (float, optional): Defaults to 1.2.

            max_bond (int, optional):  Defaults to 0.
        """
        super().__init__(weight)

        self.to_del = to_del
        self.ads_g = formula_to_graph(
            self.to_del, max_bond_ratio=max_bond_ratio, max_bond=max_bond
        )
        self.ss = surface_sites

    def node_match(self, n1, n2):
        """node matching criteria
        Args:
            n1 (str):
            n2 (str):
        Returns:
            Boolean:
        """
        return n1["symbol"] == n2["symbol"]

    def get_modified_atoms(self, atoms):
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind, atoms_g = self.ss.get_surface_sites(self.atoms)
        del df_ind
        graph_match = isomorphism.GraphMatcher(
            atoms_g, self.ads_g, node_match=self.node_match
        )
        all_isomorphisms = list(graph_match.subgraph_isomorphisms_iter())
        if not all_isomorphisms:
            raise NoReasonableStructureFound("No adsorbate in the atoms to remove")

        ind_to_remove_list = []
        for mapping in all_isomorphisms:
            matched_nodes = list(mapping.keys())
            ind_to_remove = [atoms_g.nodes[node]["index"] for node in matched_nodes]
            ind_to_remove_list.append(ind_to_remove)

        del self.atoms[ind_to_remove_list[0]]
        return self.atoms


class Swap(
    ParentModifier,
):
    """Modifier that swaps two atoms"""

    def __init__(self, weight, surface_sites, to_swap):
        """
        Args:
            weight (str):

            surface_sites (gg.SurfaceSites): A class which figures out surface sites

            to_swap (list): List of atoms to swapt
        """
        super().__init__(weight)

        self.to_swap = to_swap
        self.ss = surface_sites

    def get_modified_atoms(self, atoms):
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind,_ = self.ss.get_surface_sites(self.atoms)

        random_elem = random.sample(self.to_swap, 2)
        ind_1 = []
        ind_2 = []
        for atom in self.atoms:
            if atom.index in df_ind:
                if atom.symbol == random_elem[0]:
                    ind_1.append(atom.index)
                elif atom.symbol == random_elem[1]:
                    ind_2.append(atom.index)

        swap_1 = random.sample(ind_1, 1)
        swap_2 = random.sample(ind_2, 1)

        self.atoms[swap_1].symbol = random_elem[1]
        self.atoms[swap_2].symbol = random_elem[0]
        return
   