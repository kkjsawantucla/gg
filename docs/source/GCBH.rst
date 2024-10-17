GCBH
====

.. contents::
    :local:

Calculator
----------
The Gcbh calculator is an ase Dynamics Child on unhealthy steriods. It runs the grand canonical basin hopping, however certain functionalities are hard coded.

You can find an example for input.yaml in `gg/examples/ <https://github.com/kkjsawantucla/gg/blob/main/examples/aluminosilicate_nequip/input.yaml>`_ folder

.. code-block:: python

    from gg.gcbh import Gcbh
    from ase.io import read
    from ase.calculators.emt import EMT

    atoms = read('POSCAR')
    atoms.calc = EMT()
    G = Gcbh(atoms,config_file='input.yaml')

.. autoclass:: gg.gcbh.Gcbh
    :members: run
    :undoc-members:
    :show-inheritance:
