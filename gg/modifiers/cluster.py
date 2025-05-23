"""Import Modules for basic functioning"""

import random
import numpy as np
from ase import Atoms
from gg.utils import (
    check_contact,
    NoReasonableStructureFound,
)
from gg.utils_graph import get_connecting_nodes, get_unique_atoms
from gg.utils_cluster import (
    fit_plane,
    rotate_and_adjust_dynamically,
    translate_and_adjust,
)
from gg.modifiers.modifiers import ParentModifier
from gg.sites import Sites

__author__ = "Kaustubh Sawant"


class ClusterRotate(ParentModifier):
    """Modifier that rotates a cluster of atoms"""

    def __init__(
        self,
        surface_sites: Sites,
        max_angle: int = 180,
        rotate_vector: tuple = None,
        contact_error: float = 0.2,
        weight: float = 1,
        nmovie: int = 1,
        print_movie: bool = False,
        unique: bool = True,
        unique_method: str = "fullgraph",
        unique_depth: int = 3,
    ):
        """
        Args:
            surface_sites (gg.Sites): Class that figures out surface sites.

            max_angle (float): Maximum rotation angle allowed (degrees).
            Defaults to 180.

            rotate_vector (list/tuple): The vector along which the atoms are rotated.
            If None, the vector normal to the surface is selected.
            Defaults to None.

            contact_error: Allowable tolerance in atoms touching
            Defaults to 0.2.

            weight (float): weight for gcbh.
            Defaults to 1.
        """
        super().__init__(weight)
        self.ss = surface_sites
        self.max_angle = max_angle
        self.contact_error = contact_error
        self.normal_vector = rotate_vector
        if self.normal_vector:
            assert len(self.normal_vector) == 3
        self.print_movie = print_movie
        self.nmovie = nmovie
        self.unique = unique
        self.unique_method = unique_method
        self.unique_depth = unique_depth

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        if not df_ind:
            raise NoReasonableStructureFound("No tagged sites found")
        atoms_g = self.ss.get_graph(self.atoms)

        movie = []
        for _ in range(self.nmovie):
            a = self.atoms.copy()
            angle = np.random.uniform(0, self.max_angle)
            touch_i = get_connecting_nodes(atoms_g, df_ind, a)
            if self.normal_vector:
                normal_vector = self.normal_vector
            normal_vector = fit_plane(a, touch_i)
            a = rotate_and_adjust_dynamically(a, df_ind, normal_vector, angle)
            if check_contact(a, error=self.contact_error):
                raise NoReasonableStructureFound("Atoms Touching")
            else:
                movie.append(a)

        if not movie:
            raise NoReasonableStructureFound(
                "Movie was empty, most likely due to issues with atoms touching in Add Modifier"
            )
        if self.print_movie:
            if self.unique:
                return get_unique_atoms(
                    movie,
                    max_bond=self.ss.max_bond,
                    max_bond_ratio=self.ss.max_bond_ratio,
                    unique_method=self.unique_method,
                    depth=self.unique_depth,
                )
            else:
                return movie
        else:
            return random.sample(movie, 1)[0]


class ClusterTranslate(ParentModifier):
    """Modifier to translate a cluster of atoms in the unit cell"""

    def __init__(
        self,
        surface_sites: Sites,
        max_displace: float = 5,
        allowed_direction: tuple = (True, True, False),
        contact_error: float = 0.2,
        weight: float = 1,
        nmovie: int = 1,
        print_movie: bool = False,
        unique: bool = True,
        unique_method: str = "fullgraph",
        unique_depth: int = 3,
    ):
        """
        Args:
            surface_sites (gg.Sites): Class that figures out surface sites.

            max_displace (float): Maximum displacement allowed (angstrom).
            Defaults to 5.

            allowed_direction (tuple(bool)): To allow displacement in x,y,x direction
            Defaults to (True, True, False), surface movement

            contact_error: Allowable tolerance in atoms touching
            Defaults to 0.2.

            weight (float): weight for gcbh.
            Defaults to 1.

        """
        super().__init__(weight)
        self.ss = surface_sites
        self.max_displace = max_displace
        self.contact_error = contact_error
        self.allowed_direction = allowed_direction
        self.print_movie = print_movie
        self.nmovie = nmovie
        self.unique = unique
        self.unique_method = unique_method
        self.unique_depth = unique_depth

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        if not df_ind:
            raise NoReasonableStructureFound("No tagged sites found")

        movie = []
        for _ in range(self.nmovie):
            a = self.atoms.copy()
            random_values = (
                random.uniform(0, self.max_displace),
                random.uniform(0, self.max_displace),
                random.uniform(0, self.max_displace),
            )
            displace_vector = tuple(
                value if flag else 0
                for value, flag in zip(random_values, self.allowed_direction)
            )
            a = translate_and_adjust(a, df_ind, displace_vector, tolerance=0.2)
            if check_contact(a, error=self.contact_error):
                raise NoReasonableStructureFound("Atoms Touching")
            else:
                movie.append(a)

        if not movie:
            raise NoReasonableStructureFound(
                "Movie was empty, most likely due to issues with atoms touching in Add Modifier"
            )
        if self.print_movie:
            if self.unique:
                return get_unique_atoms(
                    movie,
                    max_bond=self.ss.max_bond,
                    max_bond_ratio=self.ss.max_bond_ratio,
                    unique_method=self.unique_method,
                    depth=self.unique_depth,
                )
            else:
                return movie
        else:
            return random.sample(movie, 1)[0]
