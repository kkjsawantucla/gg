"""test base modifiers"""

import pytest
from ase.build import fcc111, add_adsorbate
from ase import Atoms
from gg.modifiers import Remove, Replace, Swap, Rattle, Translate, ModifierAdder
from gg.predefined_sites import FlexibleSites


@pytest.fixture
def slab_with_multiple_co():
    """Fixture: Au(111) slab with multiple CO adsorbates placed on top sites."""

    slab = fcc111(
        "Au", size=(3, 3, 4), vacuum=10.0
    )  # 3x3 unit cells, 4 layers, 10 Å vacuum
    top_sites = [atom.index for atom in slab if atom.tag == 1]
    slab[top_sites[3]].symbol = "Cu"

    # Define CO molecule
    co = Atoms("CO", positions=[[0, 0, 0], [0, 0, 1.1]])

    # Add CO on top of the slab
    for i in top_sites[:2]:
        add_adsorbate(slab, co, height=2.0, position=slab[i].position[:2])

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
    ), "Expected only one unique structure after replacement"
    assert unique_structures[0].get_chemical_symbols().count("O") == 2
    assert unique_structures[0].get_chemical_symbols().count("H") == 1
    assert len(unique_structures[0]) == len(
        slab_with_multiple_co
    ), "Expected same number of atoms"


def test_swap_modifier(slab_with_multiple_co):
    """Test Swap modifier which swaps Cu with Au."""

    sites = FlexibleSites(constraints=True, contact_error=0.25, com=0.5)

    # Initialize Swap modifier
    swapper = Swap(
        surface_sites=sites,
        swap_sym=["Au", "Cu"],
        print_movie=True,
        unique=True,
    )

    # Apply modifier
    unique_structures = swapper.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert len(unique_structures) == 5, "Expected five structures after swap"

    for i in unique_structures:
        assert len(i) == len(slab_with_multiple_co), "Expected same number of atoms"


def test_rattle_modifier(slab_with_multiple_co):
    """Test rattle modifier."""

    sites = FlexibleSites(constraints=True, contact_error=0.3, com=0.5)

    # Initialize Rattle modifier
    swapper = Rattle(surface_sites=sites, seed=32, stdev=0.1)

    # Apply modifier
    structure = swapper.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert len(structure) == len(
        slab_with_multiple_co
    ), "Expected only one unique structure after rattle"


def test_translate_modifier(slab_with_multiple_co):
    """Test translate modifier."""

    sites = FlexibleSites(constraints=True, contact_error=0.3, com=0.5)

    # Initialize Translate modifier
    trans = Translate(
        surface_sites=sites,
        translate=(True, True, True),
        surf_sym=["O", "H", "C"],
        max_translate=(0.5, 0.5, 0.5),
        seed=42,
    )

    # Apply modifier
    structure = trans.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert len(structure) == len(
        slab_with_multiple_co
    ), "Expected only one unique structure after translate"


def test_modifier_adder(slab_with_multiple_co):
    """Test translate modifier."""
    sites = FlexibleSites(constraints=True,contact_error=0.25, com=0.5)

    # Initialize Remove modifier
    remover = Remove(
        surface_sites=sites,
        to_del="CO",
        max_bond_ratio=1.2,
        print_movie=True,
        unique=True,
    )

    # Initialize Swap modifier
    swapper = Swap(
        surface_sites=sites,
        swap_sym=["Au", "Cu"],
        print_movie=True,
        unique=True,
    )

    # Initialize Modifier Adder
    modifier_add = ModifierAdder([remover, swapper], print_movie=True, unique=True)

    # Apply modifier
    unique_structures = modifier_add.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert len(unique_structures) == 4, "Expected four structures after swap"

    for i in unique_structures:
        assert len(i) == len(slab_with_multiple_co) - 2, "Expected less number of atoms"
