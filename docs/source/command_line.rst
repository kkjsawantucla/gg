Command Line Tools
==================

We have created simple command line tools as a wrapper to the Add :ref:`Modifiers` so they can be used independently.

Add Mono
--------
Base: :ref:`Add Monodentate`

.. code-block:: bash

   add_mono -s <path_to_surface> -a <adsorbate|path> -sc <surface coordination> -aa <adsorbate identity>
   
   # Example
   add_mono -s POSCAR_Pt -a OH -sc 1 2 3 -aa O

Add Bi
------
Base: :ref:`Add Bidentate`

.. code-block:: bash

   add_bi -s <path_to_surface> -a <adsorbate|path> -sc <surface coordination> -aa <adsorbate identity>

   # Example using a name
   add_bi -s POSCAR_Pt -a HCOO -sc 1 2 3 -aa O

.. list-table:: Args for add_bi
   :widths: 20 10 70
   :header-rows: 0

   * - -surface / -s
     - str
     - Path to surface atoms objects to adsorb
   * - -adsorbate / -a
     - str
     - Adsorbate name or path to atoms object
   * - -ads_atom / -aa
     - str
     - Atom through which the adsorbate attaches to the surface
   * - -surf_coord / -sc
     - int
     - Coordination of surface atoms: e.g., 1-top, 2-bridge, ...
   * - -surf_sym / -ss
     - str
     - Symbols of atoms that can adsorb
   * - -unique / -u
     - str
     - Whether to check if the structures are repeated
   * - -max_bond_ratio / -mbr
     - float
     - Allowable tolerance between bonds
   * - -max_bond / -mb
     - float
     - Fixed allowable bond distance
   * - --add_ads_error / -ads_err
     - float
     - Allowable distance between two sites for adsorption
   * - --norm_method / -norm_method
     - str
     - Method to calculate normals
