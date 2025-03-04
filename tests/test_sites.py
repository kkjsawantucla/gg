"""test Site Class"""

import os
import pytest
from ase.io import read
from ase.build import fcc111, add_adsorbate
from ase import Atoms
from gg.sites import RuleSites
from gg.predefined_sites import FlexibleSites, SurfaceSites
from gg.sites import (
    get_com_sites,
    get_surface_by_normals,
    get_surface_sites_by_voronoi_pbc,
)


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


@pytest.fixture
def slab_with_cu():
    """Fixture: Cu(111) slab with CO."""

    slab = fcc111("Cu", size=(3, 3, 4), vacuum=10.0)
    top_sites = [atom.index for atom in slab if atom.tag == 1]

    # Define CO molecule
    co = Atoms("CO", positions=[[0, 0, 0], [0, 0, 1.1]])

    # Add CO on top of the slab
    add_adsorbate(slab, co, height=2.0, position=slab[top_sites[0]].position[:2])

    return slab


def test_FlexibleSites(slab_with_si_alumina):
    """Test base Flexible Sites"""
    ss = FlexibleSites(com=0.5)
    sites = ss.get_sites(slab_with_si_alumina)

    assert len(sites) == 20, "The number of sites should be 20"


def test_SurfaceSites(slab_with_si_alumina):
    """Test base Surface Sites"""
    max_coord = {"Al": 6, "Si": 4, "O": 2, "H": 1}

    # Initialize SurfaceSites
    ss = SurfaceSites(max_coord=max_coord)
    sites = ss.get_sites(slab_with_si_alumina)
    assert len(sites) == 0, "The number of sites should be 0"

    # Initialize SurfaceSites, switch off com
    ss = SurfaceSites(max_coord=max_coord, com=False)
    sites = ss.get_sites(slab_with_si_alumina)
    assert len(sites) == 8, "The number of sites should be 8"


def test_RuleSite_com(slab_with_si_alumina):
    ss = RuleSites(
        index_parsers=[
            lambda atoms: get_com_sites(atoms, fraction=0.2, direction="both")
        ],
    )
    sites = ss.get_sites(slab_with_si_alumina)
    assert len(sites) == 13, "The number of sites should be 13"


def test_RuleSite_normal(slab_with_cu):
    ss = RuleSites(
        index_parsers=[
            lambda atoms: get_surface_by_normals(atoms, rem_symbols=["C", "O"])
        ],
    )
    sites = ss.get_sites(slab_with_cu)
    assert len(sites) == 18, "The number of sites should be 18"


def test_RuleSite_voronoi(slab_with_cu):
    ss = RuleSites(
        index_parsers=[
            lambda atoms: get_surface_sites_by_voronoi_pbc(
                atoms, rem_symbols=["C", "O"]
            )
        ],
    )
    sites = ss.get_sites(slab_with_cu)
    assert len(sites) == 18, "The number of sites should be 18"


def test_RuleSite_com_normal(slab_with_cu):
    ss = RuleSites(
        index_parsers=[
            lambda atoms: get_surface_by_normals(atoms, rem_symbols=["C", "O"]),
            lambda atoms: get_com_sites(atoms, fraction=0.6, direction="above"),
        ],
        combine_rules="union",
    )
    sites = ss.get_sites(slab_with_cu)

    assert len(sites) == 20, "The number of sites should be 20"
