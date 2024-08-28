"""Import Modules for basic functioning"""

from ase.io import read as read_atoms
from ase import Atoms
from gg.utils import check_contact, generate_sites

__author__ = "Kaustubh Sawant"


class ParentModifier:
    """
    Args:
        name (str): Unique Name of the Modifier
        atoms (ase.Atoms): atoms object
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
            self._atoms = read_atoms(atoms.copy())
        elif isinstance(atoms, Atoms):
            self._atoms = atoms.copy()
        else:
            print("Please provide proper atoms file")

    def get_modified_atoms(self,atoms) -> Atoms:
        """
        Returns:
            ase.Atoms:
        """
        raise NotImplementedError


class Rattle(ParentModifier):
    """Modifier that rattles the atoms with some stdev"""

    def __init__(
        self,
        weight,
        stdev=0.001,
        contact_error=0.2
    ):
        self.stdev = stdev
        self.contact_error = contact_error
        super().__init__(weight)

    def get_modified_atoms(self,atoms):
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        self.atoms.rattle(stdev=self.stdev)
        if check_contact(self.atoms, error=self.contact_error):
            print("Atoms touching")
        return self.atoms


class Add(ParentModifier):
    """Modifier tha adds an adsorbate at certain specific sites"""

    def __init__(
        self,
        weight,
        surface_sites,
        ads,
        ads_coord,
        ad_dist=1.8,
        movie=False,
    ):
        super().__init__(weight)
        self.ss = surface_sites
        self.ads = ads
        self.ads_coord = ads_coord
        self.ad_dist = ad_dist
        self.movie = movie

    def get_modified_atoms(self,atoms):
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

        if self.movie:
            return movie
        else:
            return movie[0]


class Remove(ParentModifier):
    """Modifier that randomly removes an atom"""

    def __init__(self, weight, surface_sites):
        super().__init__(weight)
        self.ss = surface_sites

    def get_modified_atoms(self,atoms):
        """
        Returns:
            ase.Atoms:
        """
        self.atoms = atoms
        df_ind, g = self.ss.get_surface_sites(self.atoms)
        del g
        ind_to_remove = int(df_ind.to_list()[0])
        del self.atoms[ind_to_remove]
        return self.atoms
