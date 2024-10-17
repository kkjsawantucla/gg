Command Line Tools
==================

We have created simple command line tools as a wrapper to the Add modifier to use it independently.

Add Mono
--------

.. code-block:: bash

    add_mono -s <path_to_surface> -a <path_to_adsorbate> -sc <surface coordination> -aa <adsorbate identity>

    #Example
    add_mono -s POSCAR_Pt -a OH.POSCAR -sc 1 2 3 -aa O

Add Bi
--------

.. code-block:: bash

    add_bi -s <path_to_surface> -a <path_to_adsorbate> -sc <surface coordination> -aa <adsorbate identity>

    #Example
    add_bi -s POSCAR_Pt -a OCHO.POSCAR -sc 1 2 3 -aa O
