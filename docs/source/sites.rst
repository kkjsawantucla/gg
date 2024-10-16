Sites
=========

The Sites class help make graphs for modifiers to work on and also determine the surface site for modifications.

FlexibleSites
-------------

.. code-block:: python

    from gg.sites import FlexibleSites

    FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface

    atoms = read('POSCAR')
    list_sites = FS.get_sites(atoms)

.. autoclass:: gg.sites.FlexibleSites
    :members: 
    :undoc-members:
    :show-inheritance:

SurfaceSites
------------

.. code-block:: python

    from gg.sites import SurfaceSites

    max_coord = {'Pt': 12, 'O': 4, 'H': 2}
    SS = SurfaceSites(max_coord=max_coord,max_bond_ratio=1.2) #Define class to figure out surface

    atoms = read('POSCAR')
    list_sites = FS.get_sites(atoms)

.. autoclass:: gg.sites.SurfaceSites
    :members: 
    :undoc-members:
    :show-inheritance: