"""Predefined Surfaces"""

from typing import Optional
from ase import Atoms
from gg.sites import Sites
from gg.sites import (
    get_above_com_sites,
    get_surface_sites_by_coordination,
    get_tagged_sites,
    get_unconstrained_sites,
)


class FlexibleSites(Sites):
    """A class that identifies sites based on constraints, tags, or explicit indices."""

    def __init__(
        self,
        constraints: bool = False,
        index: Optional[list[int]] = None,
        tag: Optional[int] = None,
        opp_tag: bool = False,
        com: Optional[float] = None,
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
    ):
        """
        Args:
            constraints (bool, optional): If True, only unconstrained atoms are considered.
            Defaults to False.
            index (List[int], optional): A list of specific atom indices to consider.
            Defaults to None.

            tag (int, optional): If provided, only atoms with this tag are considered.
            Defaults to None.

            opp_tag (bool, optional): If True, considers atoms without the specified tag.
            Defaults to False.

            com (float, optional): If provided, considers atoms above the center of mass (COM)
                                   in the z-direction. The value determines the fraction of the
                                   z-range above the COM to consider. Defaults to None.

            max_bond_ratio (float, optional): Tolerance for bond distances. Defaults to 1.2.

            max_bond (float, optional): Maximum bond distance. Defaults to 0 (no limit).

            contact_error (float, optional): Tolerance for atoms being too close. Defaults to 0.2.
        """
        super().__init__(
            max_bond_ratio=max_bond_ratio,
            max_bond=max_bond,
            contact_error=contact_error,
        )

        # Validate inputs
        if index is None and not constraints and tag is None and com is None:
            raise RuntimeError(
                "Specify at least one of: index, constraints, tag, or com."
            )

        self.constraints = constraints
        self.index = index
        self.tag = tag
        self.opp_tag = opp_tag
        self.com = com

    def get_sites(self, atoms: Atoms) -> list[int]:
        """
        Identifies sites based on the specified rules.

        Args:
            atoms (Atoms): The ASE Atoms object.

        Returns:
            List[int]: List of atom indices identified as sites.
        """
        # Start with all indices if no specific index list is provided
        indices = self.index if self.index is not None else list(range(len(atoms)))

        # Apply constraints filter
        if self.constraints:
            indices = list(set(indices) & set(get_unconstrained_sites(atoms)))

        # Apply tag filter
        if self.tag is not None:
            if self.opp_tag:
                # Get indices of atoms without the specified tag
                tagged_indices = set(get_tagged_sites(atoms, tag=self.tag))
                indices = list(set(indices) - tagged_indices)
            else:
                # Get indices of atoms with the specified tag
                tagged_indices = set(get_tagged_sites(atoms, tag=self.tag))
                indices = list(set(indices) & tagged_indices)

        # Apply COM filter
        if self.com is not None:
            above_com_indices = set(get_above_com_sites(atoms, perc=self.com))
            indices = list(set(indices) & above_com_indices)

        return indices


class SurfaceSites(Sites):
    """A class that figures out the surface atom by coordination number"""

    def __init__(
        self,
        max_coord: dict,
        max_bond_ratio: Optional[float] = 1.2,
        max_bond: Optional[float] = 0,
        contact_error: Optional[float] = 0.2,
        com: Optional[bool] = 0.1,
    ):
        """
        Args:
        max_coord (dict): Dictionary of the maximum coordination of each element used.
        Only atoms with coordination less than this value will be considered.

        max_bond_ratio (float): Tolerance in the sum of covallent radii between two atoms to be considered a bond.
        Defaults to 1.2 (equivalent to 20% tolerance)

        max_bond (float): Maximum distance of a bond allowed, ignored if equal to zero.
        Defaults to 0

        Contact Error (float): To ensure atoms arent too close to each other, the fraction of tolerance allowed.
        Defaults to 0.2 (equivalent to 20% tolerance)

        com (Optional[bool], optional): If true atoms below the center of mass are ignored
        Defaults to True.
        """
        super().__init__(max_bond_ratio, max_bond, contact_error)
        self.max_coord = max_coord
        self.com = com
        self.df = None

    def get_sites(
        self, atoms: Atoms, self_interaction: bool = False, bothways: bool = True
    ) -> list:
        """
        Args:
            atoms (ase.Atoms): Atoms object to determine sites.

            self_interaction (bool): Input of ase.neighborlist.
            Defaults to True.

            bothways (bool): Input of ase.neighborlist.
            Defaults to False.

        Returns:
            list: list of atoms index considered for modifications.
        """
        indices = get_surface_sites_by_coordination(
            atoms,
            max_coord=self.max_coord,
            max_bond_ratio=self.max_bond_ratio,
            max_bond=self.max_bond,
            self_interaction=self_interaction,
            bothways=bothways,
        )

        # Apply COM filter
        if self.com:
            above_com_indices = set(get_above_com_sites(atoms, perc=self.com))
            indices = list(set(indices) & above_com_indices)

        return indices
