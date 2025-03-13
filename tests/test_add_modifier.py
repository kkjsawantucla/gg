"""test base modifiers"""

import pytest
import os
from ase.io import read
from ase.build import fcc111, add_adsorbate
from ase import Atoms
from gg.modifiers import Add, AddBi, ModifierAdder
from gg.predefined_sites import FlexibleSites, SurfaceSites


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


@pytest.fixture
def slab_with_si_alumina():
    """Read S-Al POSCAR"""
    # Get the directory of the current test file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to the file in the data directory
    file_path = os.path.join(current_dir, "POSCAR_SiAl")

    assert os.path.exists(file_path), f"File {file_path} does not exist"
    slab = read(file_path)

    return slab


def test_add_CO(slab_with_multiple_co):
    """Test Add modifier ensures uniqueness when adding multiple CO adsorbates."""

    # Define the Site Module
    sites = FlexibleSites(com=0.5)

    add_co = Add(
        sites,
        "CO",
        surf_coord=[1, 2, 3],
        ads_id=["C"],
        surf_sym=["Au"],
        print_movie=True,
        unique=True,
    )

    # Apply modifier
    unique_structures = add_co.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert len(unique_structures) == 27, "Expected 27 structures"
    for structure in unique_structures:
        assert (
            len(slab_with_multiple_co) == len(structure) - 2
        ), "Number of atoms should increase by 2"


def test_add_formate(slab_with_multiple_co):
    # Define the Site Module
    sites = FlexibleSites(com=0.5)

    add_oh = AddBi(
        sites,
        "HCOO",
        surf_coord=[1],
        ads_id=["O"],
        surf_sym=["Au"],
        print_movie=True,
        unique=True,
    )

    # Apply modifier
    unique_structures = add_oh.get_modified_atoms(slab_with_multiple_co)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert len(unique_structures) == 7, "Expected 7 structures"
    for structure in unique_structures:
        assert (
            len(slab_with_multiple_co) == len(structure) - 4
        ), "Number of atoms should increase by 4"


def test_add_H2O(slab_with_si_alumina):
    """Test Add modifier ensures uniqueness when adding H2O dissociatively"""

    max_coord = {"Al": 6, "Si": 4, "O": 4, "H": 1}
    # Initialize SurfaceSites, switch off com
    sites = SurfaceSites(max_coord=max_coord, com=0.5)

    add_oh = Add(
        sites,
        "OH",
        surf_coord=[1],
        ads_id=["O"],
        surf_sym=["Al", "Si"],
        print_movie=True,
        unique=True,
    )

    add_h = Add(
        sites,
        "H",
        surf_coord=[1],
        ads_id=["H"],
        surf_sym=["O"],
        print_movie=True,
        unique=True,
    )

    add_h2o = ModifierAdder([add_oh, add_h], print_movie=True, unique=True)

    # Apply modifier
    unique_structures = add_h2o.get_modified_atoms(slab_with_si_alumina)

    # Check results
    assert isinstance(unique_structures, list), "Expected a list of unique structures"
    assert len(unique_structures) == 10, "Expected structures"
    for structure in unique_structures:
        assert (
            len(slab_with_si_alumina) == len(structure) - 3
        ), "Number of atoms should increase by 2"
