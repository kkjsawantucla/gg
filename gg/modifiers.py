"""Import Modules for basic functioning"""

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
        super().__init__(weight)
        self.ss = surface_sites
        self.ads = ads
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
    """Modifier that randomly removes an atom"""

    def __init__(self, weight, surface_sites, to_del, max_bond_ratio=1.2, max_bond=0):
        super().__init__(weight)

        self.to_del = to_del
        self.ads_g = formula_to_graph(
            self.to_del, max_bond_ratio=max_bond_ratio, max_bond=max_bond
        )
        self.ss = surface_sites

    def node_match(self, n1, n2):
        """_summary_
        Args:
            n1 (str):
            n2 (str):
        Returns:
            Boolean: _description_
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
