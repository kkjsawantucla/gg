import random
from typing import Union
from itertools import combinations
from ase.collections import g2
from ase.build import molecule
from ase import Atoms
from ase.data import chemical_symbols
from gg.utils import (
    custom_copy,
    NoReasonableStructureFound,
)
from gg.utils_add import generate_add_mono, rotate_mono, rotate_bi, generate_add_bi
from gg.utils_graph import get_unique_atoms
from gg.sites import Sites
from gg.modifiers.modifiers import ParentModifier

from gg.data import adsorbates


class Add(ParentModifier):
    """Modifier that adds an monodentate adsorbate"""

    @staticmethod
    def _unique_symbols(symbols: list) -> list:
        return list(dict.fromkeys(symbols))

    def __init__(
        self,
        surface_sites: Sites,
        ads: str,
        surf_coord: int,
        surf_sym: list,
        ads_id: Union[str] = None,
        ads_dist: Union[float] = None,
        print_movie: bool = False,
        unique: bool = True,
        ads_rotate: bool = True,
        weight: float = 1,
        normal_method: str = "svd",
        unique_method: str = "fullgraph",
        unique_depth: int = 3,
        tag: bool = True,
    ):
        """
        Args:
            surface_sites (gg.Sites): Class that figures out surface sites

            ads (str) or (ase.Atoms): Adsorbate to add

            surf_coord (list[int]): How many bonds the adsorbate will make with the surface

            surf_sym (list[str]): Surface elements where adsorbate can add

            ads_id (list[float]): Strings denoting chemical symbol of adsorbate atom
            Defaults to None

            ads_dist (str, optional): Distance of adsorbate from substrate.
            Defaults to covalent radii of ads_id from ase database.

            print_movie (bool, optional): return a movie of all sites or one random site.
            Defaults to False.

            unique (bool, optional): return only unique sites.
            Defaults to True.

            normal_method (str): Determines how normals are calculated. It could be "svd" or "mean"
            Defaults to "svd"

            unique_method (str): Determines how uniqueness is calculated. User can specify atom symbols/mol
            to construct subgraphs e.g: unique_method = ["CO"]
            Defaults to "fullgraph"

            unique_depth (int): Determines the depth of subgraphs created to calculate uniqueness.
            Defaults to 3. If unique_method is "fullgraph" the value is ignored

            tag (bool): add to tag=-1 to the adsorbate (imp for clusters)
            Defaults to True.

            weight (float): weight for gcbh.
            Defaults to 1.
        """
        super().__init__(weight)
        self.ss = surface_sites
        if isinstance(ads, str):
            if ads in adsorbates:
                self.ads = adsorbates[ads]
            elif ads in g2.names:
                self.ads = molecule(ads)
            elif ads in chemical_symbols:
                self.ads = Atoms(ads, positions=[(0, 0, 0)])
            else:
                raise RuntimeError(f"Cannot convert string to Formula {ads}")
        elif isinstance(ads, Atoms):
            self.ads = custom_copy(ads)
        if isinstance(surf_sym, list):
            self.surf_sym = self._unique_symbols(surf_sym)
        else:
            self.surf_sym = surf_sym
        if isinstance(surf_coord, int):
            self.surf_coord = [surf_coord]
        elif isinstance(surf_coord, list):
            self.surf_coord = surf_coord
        else:
            raise NoReasonableStructureFound("Please enter proper value for surf_coord")
        if ads_id:
            if isinstance(ads_id, list) and all(isinstance(item, str) for item in ads_id):
                self.ads_id = self._unique_symbols(ads_id)
            else:
                self.ads_id = ads_id
        if ads_dist:
            self.ads_dist = ads_dist
        else:
            self.ads_dist = ads_id
        self.print_movie = print_movie
        self.unique = unique
        self.ads_rotate = ads_rotate
        self.method = normal_method
        self.tag = tag
        self.unique_method = unique_method
        self.unique_depth = unique_depth

    def get_all_adsorbates(self, atoms: Atoms, chem_symbol_list) -> list:
        """
        Args:
            atoms (ase.Atoms):
        Returns:
            list:
        """
        ads_list = []
        ads_dist_list = []
        for ind, atom in enumerate(atoms):
            if atom.symbol in chem_symbol_list:
                ads = rotate_mono(atoms.copy(), ind)
                ads_list.append(ads)
                if isinstance(self.ads_dist, float) or isinstance(self.ads_dist, int):
                    ads_dist_list.append(self.ads_dist)
                else:
                    ads_dist_list.append(atom.symbol)

        return ads_list, ads_dist_list

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Args:
            atoms (ase.Atoms): The atoms object on which the adsorbate will be added
        Returns:
            ase.Atoms if print_movie = True
            list[ase.Atoms] if print_movie = False
        """

        # Check multiple possibilities of adsorbate
        if isinstance(self.ads_id, list):
            if all(isinstance(item, str) for item in self.ads_id):
                ads_list, ad_dist_list = self.get_all_adsorbates(self.ads, self.ads_id)
        else:
            ads_list = [self.ads]
            ad_dist_list = [self.ads_dist]

        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        g = self.ss.get_graph(self.atoms)
        index = [ind for ind in df_ind if self.atoms[ind].symbol in self.surf_sym]

        if not index:
            raise NoReasonableStructureFound(
                "No surface sites found, check your Sites Class"
            )

        movie = []
        for i, ads in enumerate(ads_list):
            # Read gg.utils_add to understand the working
            movie += generate_add_mono(
                self.atoms,
                ads,
                g,
                index,
                self.surf_coord,
                ad_dist=ad_dist_list[i],
                contact_error=self.ss.contact_error,
                method=self.method,
                tag=self.tag,
            )
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


class AddBi(Add):
    """Modifier that adds an adsorbate at certain specific sites"""

    def __init__(
        self,
        surface_sites: Sites,
        ads: str,
        surf_coord: int,
        surf_sym: list,
        ads_id: list,
        ads_dist: Union[float, str] = None,
        print_movie: bool = False,
        unique: bool = True,
        ads_rotate: bool = True,
        add_ads_error: float = 0.5,
        normal_method: str = "mean",
        tag: bool = True,
        unique_method: str = "fullgraph",
        unique_depth: int = 3,
        weight: float = 1,
    ):
        """
        Args:

            surface_sites (gg.Sites): Class that figures out surface sites.

            ads (str) or (ase.Atoms): Adsorbate to add.

            surf_coord (list[int]): How many bonds the adsorbate will make with the surface.

            surf_sym (list[str]): Surface elements where adsorbate can add.

            ads_id (list of [int or str]): Strings denoting chemical symbol of adsorbate atom.

            ads_dist (list of [float or str, optional]): Distance of adsorbate from surface site.
            If its string denoting chemical symbol of adsorbate atom,
            then distance is set by atomic radii.
            Defaults to covalent radii of atoms mentioned in ads_id.

            print_movie (bool, optional): Return a movie of all sites or one random site.
            Defaults to False.

            unique (bool, optional): Return only unique sites.
            Defaults to True.

            ads_rotate (bool,optional): Rotate atoms such that they point in +z direction.
            Defaults to True.

            add_ads_error (float): The error in distance between bidentate adsorbate sites.
            Defaults to 0.5 (equivalent to 50%)

            normal_method (str): Determines how normals are calculated. It could be "svd" or "mean"
            Defaults to "mean"

            unique_method (str): Determines how uniqueness is calculated. User can specify atom symbols/mol
            to construct subgraphs e.g: unique_method = ["C"]
            Defaults to "fullgraph"

            unique_depth (int): Determines the depth of subgraphs created to calculate uniqueness.
            Defaults to 3. If unique_method is "fullgraph" the value is ignored

            tag (bool): add to tag=-1 to the adsorbate (imp for clusters)
            Defaults to 1.

            weight (float): weight for gcbh.
            Defaults to 1.
        """
        super().__init__(
            surface_sites=surface_sites,
            ads=ads,
            surf_coord=surf_coord,
            surf_sym=surf_sym,
            ads_dist=ads_dist,
            print_movie=print_movie,
            unique=unique,
            ads_rotate=ads_rotate,
            weight=weight,
            normal_method=normal_method,
            tag=tag,
            unique_method=unique_method,
            unique_depth=unique_depth,
        )

        if isinstance(ads_id, list) and all(isinstance(item, str) for item in ads_id):
            self.ads_id_list = self._unique_symbols(ads_id)
        else:
            self.ads_id_list = ads_id

        if all(isinstance(item, str) for item in self.ads_id_list):
            self.ads_id_list, self.ads, self.ads_dist_list = self.get_all_adsorbates(
                self.ads, self.ads_id_list
            )

        self.ads_add_error = add_ads_error

    def get_all_adsorbates(self, atoms: Atoms, chem_symbol_list) -> list:
        """
        Args:
            atoms (ase.Atoms):
        """
        list_ads = []
        ads_list = []
        ads_dist_list = []
        possible = list(combinations(range(len(atoms)), 2))
        for ind in possible:
            ind_1, ind_2 = ind
            ads_id = [ind_1, ind_2]
            if (
                atoms[ind_1].symbol in chem_symbol_list
                and atoms[ind_2].symbol in chem_symbol_list
            ):
                ads = rotate_bi(atoms.copy(), ads_id)

                list_ads.append(ads_id)
                ads_list.append(ads)

                if isinstance(self.ads_dist, float) or isinstance(self.ads_dist, int):
                    ads_dist_list.append([self.ads_dist, self.ads_dist])
                else:
                    ads_dist_list.append(
                        [atoms[ads_id[0]].symbol, atoms[ads_id[1]].symbol]
                    )

        return list_ads, ads_list, ads_dist_list

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        """
        Args:
            atoms (ase.Atoms): The atoms object on which the adsorbate will be added
        Returns:
            ase.Atoms if print_movie = True,
            list[ase.Atoms] if print_movie = False
        """
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        g = self.ss.get_graph(self.atoms)
        index = [ind for ind in df_ind if self.atoms[ind].symbol in self.surf_sym]

        if not index:
            raise NoReasonableStructureFound(
                "No surface sites found, check your Sites Class"
            )
        movie = []

        for i, ads_id in enumerate(self.ads_id_list):
            # Read gg.utils_add to understand the working
            movie += generate_add_bi(
                self.atoms,
                self.ads[i],
                g,
                index,
                self.surf_coord,
                ad_dist=self.ads_dist_list[i],
                ads_index=ads_id,
                contact_error=self.ss.contact_error,
                ads_add_error=self.ads_add_error,
                method=self.method,
                tag=self.tag,
            )
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
                    depth=self.unique_depth,
                )
            else:
                return movie
        else:
            return random.sample(movie, 1)[0]
