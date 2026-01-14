Sites
=========

The Sites class help make graphs for modifiers to work on and also determine the surface site for modifications.

FlexibleSites
-------------

This is a simple sites class which returns either atoms which arent constrained as surface sites or you can specify specific index.
Hardcoding of indexes isnt advisable as the atoms object changes during gcbh runs.

.. code-block:: python

    from gg.predefined_sites import FlexibleSites

    FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface

    atoms = read('POSCAR')
    list_sites = FS.get_sites(atoms)

.. autoclass:: gg.predefined_sites.FlexibleSites
    :members: 
    :undoc-members:
    :show-inheritance:

SurfaceSites
------------

This class uses co-ordination number to determine the surface sites. However, we need to specify the maximum co-ordination allowed for each atom.

.. code-block:: python

    from gg.predefined_sites import SurfaceSites

    max_coord = {'Pt': 12, 'O': 4, 'H': 2}
    SS = SurfaceSites(max_coord=max_coord,max_bond_ratio=1.2) #Define class to figure out surface

    atoms = read('POSCAR')
    list_sites = FS.get_sites(atoms)

.. autoclass:: gg.predefined_sites.SurfaceSites
    :members: 
    :undoc-members:
    :show-inheritance:


RuleSites
------------

This class can define complex Sites based on user defined rules.

    .. code-block:: python

        from gg.sites import RuleSites, get_com_sites, get_surface_sites_by_coordination

        #Define maximum co-odrination each species can have
        max_coord = {"Al": 6, "Si": 4, "O": 4, "H": 1}

        ss = RuleSites(
            index_parsers=[
                lambda atoms: get_com_sites(atoms, fraction=0.50, direction="above"),
                lambda atoms: get_surface_sites_by_coordination(
                    atoms, max_coord, max_bond=2,
                ),
            ],
            combine_rules="intersection",
        )

    .. autoclass:: gg.sites.RuleSites
    :members: 
    :undoc-members:
    :show-inheritance:
