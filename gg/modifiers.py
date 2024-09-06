"""Import Modules for basic functioning"""

import random
from itertools import product
from ase.io import read as read_atoms
from ase import Atoms
from ase.data import covalent_radii
from networkx.algorithms import isomorphism
from gg.utils import (
    check_contact,
    generate_sites,
    custom_copy,
    formula_to_graph,
    move_along_normal,
    NoReasonableStructureFound,
)
from gg.utils_graph import get_unique_atoms
from gg.sites import SurfaceSites

__author__ = "Kaustubh Sawant"


class ParentModifier:
    """Parent bare bones modifier which serves as the basis for other modifiers
    Args:
        weight (float): Modifier Weight
    """

    def __init__(self, weight: float):
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
    def atoms(self, atoms: Atoms):
        if isinstance(atoms, str):
            self._atoms = read_atoms(atoms)
        elif isinstance(atoms, Atoms):
            self._atoms = custom_copy(atoms)

    def get_modified_atoms(self, atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class ModifierAdder(ParentModifier):
    """Add modifiers together to create a new modifier"""

    def __init__(self, weight: float, modifier_instances: list[SurfaceSites]):
        super().__init__(weight)
        if isinstance(modifier_instances, list):
            self.modifier_instances = modifier_instances
        else:
            raise RuntimeError("modifier_instances isn't a list. Please provide a list")

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        for instance in self.modifier_instances:
            self.atoms = instance.get_modified_atoms(self.atoms)
        return self.atoms


class Rattle(ParentModifier):
    """Modifier that rattles the atoms with some stdev"""

    def __init__(self, weight: float, stdev: float = 0.001, contact_error: float = 0.2):
        self.stdev = stdev
        self.contact_error = contact_error
        super().__init__(weight)

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
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
        weight: float,
        surface_sites: SurfaceSites,
        ads: str,
        ads_coord: int,
        surf_sym: list,
        ad_dist: float = 1.8,
        print_movie: bool = False,
    ):
        """
        Args:
            weight (float): _description_
            surface_sites (gg.SurfaceSites): a class that figures out surface sites
            ads (str) or (ase.Atoms): Adsorbate to add
            ads_coord (int): Adsorbate coordination number for adding
            surf_sym (list): Surface elements where adsorbate can add
            ad_dist (float, optional): Distance of adsorbate from surface site.
            Defaults to 1.8.
            movie (bool, optional): return a movie of all unique sites or one random site.
            Defaults to False.
        """
        super().__init__(weight)
        self.ss = surface_sites
        if isinstance(ads, str):
            self.ads = read_atoms(ads)
        elif isinstance(ads, Atoms):
            self.ads = custom_copy(ads)
        self.surf_sym = surf_sym
        self.ads_coord = ads_coord
        self.ad_dist = ad_dist
        self.print_movie = print_movie

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        g = self.ss.graph
        index = [ind for ind in df_ind if self.atoms[ind].symbol in self.surf_sym]

        movie = generate_sites(
            self.atoms,
            self.ads,
            g,
            index,
            self.ads_coord,
            ad_dist=self.ad_dist,
            contact_error=self.ss.contact_error,
        )
        if not movie:
            raise NoReasonableStructureFound(
                "Movie was empty, most likely due to issues with atoms touching in Add Modifier"
            )
        if self.print_movie:
            return get_unique_atoms(
                movie, max_bond=self.ss.max_bond, max_bond_ratio=self.ss.max_bond_ratio
            )
        else:
            return random.sample(movie, 1)[0]


class Remove(
    ParentModifier,
):
    """Modifier that randomly removes an atom or molecule"""

    def __init__(
        self,
        weight: float,
        surface_sites: SurfaceSites,
        to_del: str,
        max_bond_ratio: float = 1.2,
        max_bond: float = 0,
        print_movie: bool = False,
    ):
        """
        Args:
            weight (str):

            surface_sites (gg.SurfaceSites): a class that figures out surface sites

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
        self.print_movie = print_movie

    def node_match(self, n1: str, n2: str):
        """node matching criteria
        Args:
            n1 (str):
            n2 (str):
        Returns:
            Boolean:
        """
        return n1["symbol"] == n2["symbol"]

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        atoms_g = self.ss.graph
        graph_match = isomorphism.GraphMatcher(
            atoms_g, self.ads_g, node_match=self.node_match
        )
        all_isomorphisms = list(graph_match.subgraph_isomorphisms_iter())
        if not all_isomorphisms:
            raise NoReasonableStructureFound(
                "No adsorbate in the atoms to remove in Remove Modifier"
            )

        ind_to_remove_list = []
        for mapping in all_isomorphisms:
            matched_nodes = list(mapping.keys())
            ind_to_remove = [atoms_g.nodes[node]["index"] for node in matched_nodes]
            if all(element in df_ind for element in ind_to_remove):
                ind_to_remove_list.append(ind_to_remove)
            else:
                continue

        if not ind_to_remove_list:
            raise NoReasonableStructureFound()

        if self.print_movie:
            movie = []
            for ind_to_remove in ind_to_remove_list:
                a = custom_copy(self.atoms)
                del a[ind_to_remove]
                movie.append(a)
            return movie
        else:
            random_remove = random.sample(ind_to_remove_list, 1)[0]
            del self.atoms[random_remove]
            return self.atoms


class Swap(
    ParentModifier,
):
    """Modifier that swaps two atoms"""

    def __init__(
        self,
        weight: float,
        surface_sites: SurfaceSites,
        swap_sym: list,
        swap_ind: list = None,
        print_movie: bool = False,
    ):
        """
        Args:
            weight (str):

            surface_sites (gg.SurfaceSites): A class which figures out surface sites

            swap_sym (list): List of atoms to swap

            swap_ind (list): List of indices to swap. Default to None
        """
        super().__init__(weight)
        self.swap_sym = swap_sym
        self.swap_ind = swap_ind
        self.ss = surface_sites
        self.print_movie = print_movie

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms

        if self.swap_ind:
            if len(self.swap_ind) == 2:
                swap_1 = self.swap_ind[0]
                swap_2 = self.swap_ind[1]
                random_elem = [self.atoms[swap_1].symbol, self.atoms[swap_2].symbol]
            else:
                raise RuntimeError("Multiple indices given to Swap Modifier")
        else:
            df_ind = self.ss.get_sites(self.atoms)
            random_elem = random.sample(self.swap_sym, 2)

            ind_1 = []
            ind_2 = []
            for atom in self.atoms:
                if atom.index in df_ind:
                    if atom.symbol == random_elem[0]:
                        ind_1.append(atom.index)
                    elif atom.symbol == random_elem[1]:
                        ind_2.append(atom.index)
                    else:
                        continue
        movie = []
        combinations = product(ind_1, ind_2)
        for comb in combinations:
            atoms = custom_copy(self.atoms)
            swap_1, swap_2 = comb
            symbols = atoms.get_chemical_symbols()
            symbols[swap_1] = random_elem[1]
            symbols[swap_2] = random_elem[0]
            atoms.set_chemical_symbols(symbols)
            if (
                covalent_radii[atoms[swap_1].number]
                >= covalent_radii[atoms[swap_2].number]
            ):
                index = swap_1
            else:
                index = swap_2
            self.ss.get_graph(atoms)
            g = self.ss.graph
            atoms = move_along_normal(index, atoms, g)

            if check_contact(atoms, error=self.ss.contact_error):
                del atoms
                continue
            else:
                movie.append(atoms)
                del atoms
        if self.print_movie:
            return get_unique_atoms(
                movie, max_bond=self.ss.max_bond, max_bond_ratio=self.ss.max_bond_ratio
            )
        else:
            return random.sample(movie, 1)[0]
