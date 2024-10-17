Command Line Tools
==================

We have created simple command line tools as a wrapper to the Add :ref:`Modifiers` to use it independently.

Add Mono
--------
Base: :ref:`Add Monodentate`

.. code-block:: bash

    add_mono -s <path_to_surface> -a <path_to_adsorbate> -sc <surface coordination> -aa <adsorbate identity>

    #Example
    add_mono -s POSCAR_Pt -a OH.POSCAR -sc 1 2 3 -aa O

Args

.. list-table:: Args for add_mono
    :widths: 25 25 50

    * - -surface/-s
     -str
     -Path to surface atoms objects to adsorb
   * - -adsorbate/-a
     -str
     -Path to adsorbate atoms objects to process
   * - -ads_atom/-aa
     -str
     -Atom through which the adsorbate attaches to the surface
   * - -surf_coord/-sc
     -int
     -Coordination of surface atoms: eg. 1-top, 2-bridge,...

   * - -surf_sym/-ss
     -str
     -Symbols of atoms that can adsorb
   * - -unique/-u
     -str
     -Whether to check if the structures are repeated

   * - -max_bond_ratio/-mbr
     -float
     -Allowable tolerance between bonds

   * --max_bond/-mb 
     -float
     -Fixed allowable bond distance


Add Bi
--------
Base: :ref:`Add Bidentate`

.. code-block:: bash

    add_bi -s <path_to_surface> -a <path_to_adsorbate> -sc <surface coordination> -aa <adsorbate identity>

    #Example
    add_bi -s POSCAR_Pt -a OCHO.POSCAR -sc 1 2 3 -aa O
