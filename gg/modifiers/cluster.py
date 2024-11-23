"""Import Modules for basic functioning"""

import numpy as np
from ase import Atoms
from gg.utils import (
    check_contact,
    NoReasonableStructureFound,
)
from gg.utils_graph import get_connecting_nodes
from gg.utils_cluster import fit_plane, rotate_and_adjust_dynamically
from gg.modifiers.modifiers import ParentModifier
from gg.sites import SurfaceSites

__author__ = "Kaustubh Sawant"


class ClusterRotate(ParentModifier):
    """Modifier that adds an monodentate adsorbate"""

    def __init__(
        self,
        surface_sites: SurfaceSites,
        rotate_axis: tuple = (True, True, True),
        max_angle: int = 180,
        contact_error: float = 0.2,
        weight: float = 1,
    ):
        """
        Args:

        """
        super().__init__(weight)
        self.ss = surface_sites
        self.translate = rotate_axis
        self.max_angle = max_angle
        self.contact_error = contact_error

    def get_modified_atoms(self, atoms: Atoms) -> Atoms:
        self.atoms = atoms
        df_ind = self.ss.get_sites(self.atoms)
        atoms_g = self.ss.get_graph(self.atoms)

        angle = np.random.uniform(0, self.max_angle)
        touch_i = get_connecting_nodes(atoms_g,df_ind,self.atoms)
        normal_vector = fit_plane(self.atoms, touch_i)
        self.atoms = rotate_and_adjust_dynamically(self.atoms,df_ind,normal_vector,angle)

        if check_contact(self.atoms, error=self.contact_error):
            raise NoReasonableStructureFound("Atoms Touching")
        else:
            return self.atoms
