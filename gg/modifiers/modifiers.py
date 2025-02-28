"""Import Modules for basic functioning"""

import random
from typing import Union
from itertools import product
import numpy as np
from ase.io import read as read_atoms
from ase import Atoms
from ase.data import covalent_radii
from networkx.algorithms import isomorphism
from gg.utils import (
    check_contact,
    replace,
    custom_copy,
    formula_to_graph,
    move_along_normal,
    NoReasonableStructureFound,
)
from gg.utils_graph import get_unique_atoms
from gg.sites import Sites

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
        """Accept only string or ase.Atoms object"""
        if isinstance(atoms, str):
            self._atoms = read_atoms(atoms)
        elif isinstance(atoms, Atoms):
            self._atoms = custom_copy(atoms)
        else:
            raise RuntimeError("Cannot Set Atoms")

    @atoms.deleter
    def atoms(self):
        """Deletes the atoms object and resets it to None."""
        self._atoms = None

    def get_modified_atoms(self, atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class ModifierAdder(ParentModifier):
    """Add modifiers together to create a new modifier"""

    def __init__(
        self,
        modifier_instances: list,
        weight: float = 1,
        max_bond_ratio: float = 1.2,
        max_bond: float = 0,
        print_movie: bool = False,
        unique: bool = True,
    ):
        super().__init__(weight)
        if isinstance(modifier_instances, list):
            self.modifier_instances = modifier_instances
        else:
            raise RuntimeError("modifier_instances isn't a list. Please provide a list")
        self.max_bond_ratio = max_bond_ratio
        self.max_bond = max_bond
        self.print_movie = print_movie
        self.unique = unique

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        if self.print_movie:
            movie = [self.atoms]
            for instance in self.modifier_instances:
                new_movie = []
                for atoms in movie:
                    new_atoms = instance.get_modified_atoms(atoms)
                    if isinstance(new_atoms, Atoms):
                        new_movie.append(new_atoms)
                    elif isinstance(new_atoms, list):
                        new_movie = new_movie + new_atoms
                    else:
                        continue
                movie = new_movie
            if not movie:
                raise NoReasonableStructureFound("Movie was empty")
            if self.unique:
                return get_unique_atoms(
                    movie,
                    max_bond=self.max_bond,
                    max_bond_ratio=self.max_bond_ratio,
                )
            else:
                return movie
        else:
            atoms = self.atoms
            for instance in self.modifier_instances:
                atoms = instance.get_modified_atoms(atoms)
                if isinstance(atoms, list):
                    raise NoReasonableStructureFound(
                        f"Switch off Print Movie for instance {instance}"
                    )
            return atoms


class Rattle(ParentModifier):
    """Modifier that rattles the atoms with some stdev"""

    def __init__(
        self, weight: float = 1, stdev: float = 0.001, contact_error: float = 0.2
    ):
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


class Translate(ParentModifier):
    """Modifier that adds an monodentate adsorbate"""

    def __init__(
        self,
        surface_sites: Sites,
        translate: tuple = (True, True, True),
        max_translate: tuple = (0.2, 0.2, 0.2),
        surf_sym: list = None,
        pick_random: bool = False,
        contact_error: float = 0.2,
        weight: float = 1,
    ):
        """
        Args:

        """
        super().__init__(weight)
        self.ss = surface_sites
        self.translate = translate
        self.max_translate = max_translate
        self.surf_sym = surf_sym
        self.ran = pick_random
        self.contact_error = contact_error

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        if self.surf_sym:
            index = [ind for ind in df_ind if self.atoms[ind].symbol in self.surf_sym]
        else:
            index = df_ind

        if self.ran:
            index = [np.random.choice]

        disp = [
            np.random.uniform(-max_val, max_val) if truth else 0
            for truth, max_val in zip(self.translate, self.max_translate)
        ]

        for i in index:
            self.atoms[i].position += disp

        if check_contact(self.atoms, error=self.contact_error):
            raise NoReasonableStructureFound("Atoms Touching")
        else:
            return self.atoms


class Remove(
    ParentModifier,
):
    """Modifier that randomly removes an atom or molecule"""

    def __init__(
        self,
        surface_sites: Sites,
        to_del: Union[Atoms, str],
        max_bond_ratio: float = 1.2,
        max_bond: float = 0,
        print_movie: bool = False,
        unique: bool = False,
        weight: float = 1,
    ):
        """
        Args:
            surface_sites (gg.Sites): Class that figures out surface sites

            to_del (str) or (ase.Atoms): Atoms to delete. If a string is provided,
            it tries to make a molecule out of it.

            max_bond_ratio (float, optional): Max bond ratio to make graph of atoms to delete.
            Defaults to 1.2.

            max_bond (int, optional): Max bond ratio to make graph of atoms to delete.
            Defaults to 0.

            print_movie (bool, optional): Return a movie of all sites or one random site.
            Defaults to False.

            unique (bool, optional): Return only unique sites.

            weight (float): weight for gcbh.
            Defaults to 1.
        """
        super().__init__(weight)

        self.to_del = to_del

        # Make graph for the adsorbate
        self.ads_g = formula_to_graph(
            self.to_del, max_bond_ratio=max_bond_ratio, max_bond=max_bond
        )
        self.ss = surface_sites
        self.print_movie = print_movie
        self.unique = unique

    def node_match(self, n1: str, n2: str):
        """node matching criteria
        Args:
            n1 (str):
            n2 (str):
        Returns:
            Boolean:
        """
        return n1["symbol"] == n2["symbol"]

    def get_ind_to_remove_list(self, atoms: Atoms) -> list:
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        atoms_g = self.ss.get_graph(self.atoms)

        # Check if the adsorbate graph exists in atoms graph
        graph_match = isomorphism.GraphMatcher(
            atoms_g, self.ads_g, node_match=self.node_match
        )
        all_isomorphisms = list(graph_match.subgraph_isomorphisms_iter())
        if not all_isomorphisms:
            raise NoReasonableStructureFound(
                "No adsorbate in the atoms to remove in Remove Modifier"
            )

        # Figure out the indices of the atoms to remove
        ind_to_remove_list = []
        for mapping in all_isomorphisms:
            matched_nodes = list(mapping.keys())
            ind_to_remove = [atoms_g.nodes[node]["index"] for node in matched_nodes]
            if all(element in df_ind for element in ind_to_remove):
                ind_to_remove_list.append(ind_to_remove)
            else:
                continue
        # Check its not empty
        if not ind_to_remove_list:
            raise NoReasonableStructureFound(
                "Index of the atoms to be removed isnt in Site Class"
            )
        return ind_to_remove_list

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Args:
            atoms (ase.Atoms): The atoms object on which the adsorbate will be added
        Returns:
            ase.Atoms if print_movie = True
            list[ase.Atoms] if print_movie = False
        """
        ind_to_remove_list = self.get_ind_to_remove_list(atoms)

        if self.print_movie:
            movie = []
            for ind_to_remove in ind_to_remove_list:
                a = custom_copy(self.atoms)
                del a[ind_to_remove]
                movie.append(a)
            if self.unique:
                return get_unique_atoms(
                    movie,
                    max_bond=self.ss.max_bond,
                    max_bond_ratio=self.ss.max_bond_ratio,
                )
            else:
                return movie
        else:
            random_remove = random.sample(ind_to_remove_list, 1)[0]
            del self.atoms[random_remove]
            return self.atoms

    def get_n(self, atoms: Atoms) -> int:
        """
        Args:
            atoms (Atoms):

        Returns:
            int: Instances a particular subgraph is seen
        """
        ind_to_remove_list = self.get_ind_to_remove_list(atoms)
        return len(ind_to_remove_list)


class Swap(
    ParentModifier,
):
    """Modifier that swaps two atoms"""

    def __init__(
        self,
        surface_sites: Sites,
        swap_sym: list,
        swap_ind: list = None,
        print_movie: bool = False,
        unique: bool = True,
        weight: float = 1,
    ):
        """
        Args:
            surface_sites (gg.Sites): Class which figures out surface sites

            swap_sym (list): List of atom symbols that are allowed to swap

            swap_ind (list): List of indices to swap.
            Default to None.

            print_movie (bool, optional): Return a movie of all sites or one random site.
            Defaults to False.

            unique (bool, optional): Return only unique sites.
            Defaults to True

            weight (float): weight for gcbh.
            Defaults to 1.
        """
        super().__init__(weight)
        self.swap_sym = swap_sym
        self.swap_ind = swap_ind
        self.ss = surface_sites
        self.print_movie = print_movie
        self.unique = unique

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Args:
            atoms (ase.Atoms): The atoms object on which the adsorbate will be added
        Returns:
            ase.Atoms if print_movie = True
            list[ase.Atoms] if print_movie = False
        """
        self.atoms = atoms

        if self.swap_ind:
            if len(self.swap_ind) == 2:
                ind_1 = [self.swap_ind[0]]
                ind_2 = [self.swap_ind[1]]
                random_elem = [self.atoms[ind_1[0]].symbol, self.atoms[ind_2[0]].symbol]
            else:
                raise RuntimeError("Multiple indices given to Swap Modifier")
        else:
            df_ind = self.ss.get_sites(self.atoms)

            # Randomly select two elements to swap
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

        # select comination of indices for the two elements
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
            g = self.ss.get_graph(atoms)

            # Move atoms if there is significant diff in covalent radii
            atoms = move_along_normal(index, atoms, g)

            if check_contact(atoms, error=self.ss.contact_error):
                del atoms
                continue
            else:
                movie.append(atoms)
                del atoms
        if not movie:
            raise NoReasonableStructureFound("Movie was empty")
        if self.print_movie:
            if self.unique:
                return get_unique_atoms(
                    movie,
                    max_bond=self.ss.max_bond,
                    max_bond_ratio=self.ss.max_bond_ratio,
                )
            else:
                return movie
        else:
            return random.sample(movie, 1)[0]


class Replace(
    Remove,
):
    """Modifier that replaces one atoms object with another"""

    def __init__(
        self,
        surface_sites: Sites,
        to_del: Union[Atoms, str],
        with_replace: Union[Atoms, str],
        max_bond_ratio: float = 1.2,
        max_bond: float = 0,
        print_movie: bool = False,
        unique: bool = True,
        weight: float = 1,
    ):
        """
        Args:
            surface_sites (gg.Sites): Class that figures out surface sites

            to_del (str) or (ase.Atoms): Atoms to delete. If a string is provided,
            it tries to make a molecule out of it.

            with_replace (str) or (ase.Atoms): Atoms to replace with. If a string is provided,
            it tries to make a molecule out of it.

            max_bond_ratio (float, optional): Max bond ratio to make graph of atoms to delete.
            Defaults to 1.2.

            max_bond (int, optional): Max bond ratio to make graph of atoms to delete.
            Defaults to 0.

            print_movie (bool, optional): Return a movie of all sites or one random site.
            Defaults to False.

            unique (bool, optional): Return only unique sites.

            weight (float): weight for gcbh.
            Defaults to 1.
        """
        super().__init__(
            surface_sites=surface_sites,
            to_del=to_del,
            max_bond_ratio=max_bond_ratio,
            max_bond=max_bond,
            print_movie=print_movie,
            unique=unique,
            weight=weight,
        )

        self.with_rep = with_replace

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Args:
            atoms (ase.Atoms): The atoms object on which the adsorbate will be added
        Returns:
            ase.Atoms if print_movie = True
            list[ase.Atoms] if print_movie = False
        """
        ind_to_remove_list = self.get_ind_to_remove_list(atoms)

        if self.print_movie:
            movie = []
            for ind_to_remove in ind_to_remove_list:
                a = custom_copy(self.atoms)
                positions = a.get_positions()[ind_to_remove]
                offset = np.mean(positions, axis=0)
                del a[ind_to_remove]
                a = replace(a, self.with_rep, offset)
                movie.append(a)
            if self.unique:
                return get_unique_atoms(
                    movie,
                    max_bond=self.ss.max_bond,
                    max_bond_ratio=self.ss.max_bond_ratio,
                )
            else:
                return movie
        else:
            random_remove = random.sample(ind_to_remove_list, 1)[0]
            positions = self.atoms.get_positions()[random_remove]
            offset = np.mean(positions, axis=0)
            del self.atoms[random_remove]
            self.atoms = replace(self.atoms, self.with_rep, offset)
            return self.atoms
