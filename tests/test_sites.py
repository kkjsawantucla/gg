"""test Site Class"""

import os
import pytest
from ase.io import read
from gg.sites import RuleSites, get_com_sites
from gg.predefined_sites import FlexibleSites, SurfaceSites


@pytest.fixture
def slab_with_si_alumina():
    """Read POSCAR"""
    # Get the directory of the current test file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to the file in the data directory
    file_path = os.path.join(current_dir,"POSCAR_SiAl")

    assert os.path.exists(file_path), f"File {file_path} does not exist"
    slab = read(file_path)

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


def test_RuleSitecom(slab_with_si_alumina):
    ss = RuleSites(
        index_parsers=[
            lambda atoms: get_com_sites(atoms, fraction=0.2, direction="both")
        ],
        max_bond=2,
        max_bond_ratio=0.0,
        contact_error=0.2,
    )
    sites = ss.get_sites(slab_with_si_alumina)
    assert len(sites) == 13, "The number of sites should be 13"
