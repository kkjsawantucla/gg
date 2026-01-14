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

    atoms = read("POSCAR")
    list_sites = FS.get_sites(atoms)

Common usage patterns
^^^^^^^^^^^^^^^^^^^^^

* Use ``constraints=True`` when the slab is fixed with ``ase.constraints.FixAtoms``.
* Use ``tag``/``opp_tag`` to include or exclude atoms with a specific tag (e.g., tag surface atoms in your input file).
* Use ``com`` (e.g., ``com=0.6``) to focus on atoms above the center of mass for clusters or slabs.

.. autoclass:: gg.predefined_sites.FlexibleSites
    :members: 
    :undoc-members:
    :show-inheritance:

SurfaceSites
------------

This class uses co-ordination number to determine the surface sites. However, we need to specify the maximum co-ordination allowed for each atom.

.. code-block:: python

    from gg.predefined_sites import SurfaceSites

    max_coord = {"Pt": 12, "O": 4, "H": 2}
    SS = SurfaceSites(max_coord=max_coord,max_bond_ratio=1.2) #Define class to figure out surface

    atoms = read("POSCAR")
    list_sites = SS.get_sites(atoms)

Notes
^^^^^

* ``com`` expects a fraction (e.g., ``0.1`` keeps atoms above the COM by 10% of the z-range).
* Set ``com=None`` to keep all atoms regardless of height.

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
