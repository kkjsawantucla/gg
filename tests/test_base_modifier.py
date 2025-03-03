"""test base modifiers"""

import pytest
from ase.build import fcc111, add_adsorbate
from ase import Atoms
from gg.modifiers import Remove, Replace
from gg.predefined_sites import FlexibleSites


@pytest.fixture
def slab_with_multiple_co():
    """Fixture: Au(111) slab with multiple CO adsorbates placed on top sites."""

    slab = fcc111(
        "Au", size=(3, 3, 4), vacuum=10.0
    )  # 3x3 unit cells, 4 layers, 10 Å vacuum
    top_sites = [atom.index for atom in slab if atom.tag == 1]

    # Define CO molecule
    CO = Atoms("CO", positions=[[0, 0, 0], [0, 0, 1.1]])

    # Add CO on top of the slab
    for i in top_sites[:2]:
        add_adsorbate(slab, CO, height=2.0, position=slab[i].position[:2])

    return slab


def test_remove_modifier_uniqueness(slab_with_multiple_co):
    """Test Remove modifier ensures uniqueness when removing multiple CO adsorbates."""

    sites = FlexibleSites(constraints=True)

    # Initialize Remove modifier with print_movie=True to get all possible removals
    remover = Remove(
        surface_sites=sites,
        to_del="CO",
        max_bond_ratio=1.2,
        print_movie=True,
        unique=True,
    )

    # Apply modifier
    unique_structures = remover.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert (
        len(unique_structures) == 1
    ), "Expected only one unique structure after removal"
    assert (
        len(slab_with_multiple_co) == len(unique_structures[0]) + 2
    ), "Number of atoms should reduce by 2"


def test_replace_modifier(slab_with_multiple_co):
    """Test Replace modifier swaps CO with OH."""

    sites = FlexibleSites(constraints=True)

    # Initialize Replace modifier
    replacer = Replace(
        surface_sites=sites,
        to_del="CO",
        with_replace="OH",
        max_bond_ratio=1.2,
        print_movie=True,
        unique=True,
    )

    # Apply modifier
    unique_structures = replacer.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert (
        len(unique_structures) == 1
    ), "Expected only one unique structure after removal"
    assert unique_structures[0].get_chemical_symbols().count("O") == 2
    assert unique_structures[0].get_chemical_symbols().count("H") == 1
    assert len(unique_structures[0]) == len(
        slab_with_multiple_co
    ), "Expected same number of atoms"
