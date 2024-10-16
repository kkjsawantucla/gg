Usage
=====

.. contents::
   :local:

GCBH
----
The Gcbh calculator is an ase Dynamics Child on unhealthy steriods. It runs the grand canonical basin hopping, however certain functionalities are hard coded.

.. code-block:: python

   from gg.gcbh import Gcbh
   from ase.io import read
   from ase.calculators.emt import EMT
   
   atoms = read('POSCAR')
   atoms.calc = EMT()
   G = Gcbh(atoms,config_file='input.yaml')


Inputs
------
1. atoms (ase.Atoms): An `ase Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ object as a starting point. The object should have a `ase.calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html>`_ attached. 
2. logfile (str): path to a file that logs the calculator's output.
3. trajectory (str): path to a file that logs all the atoms structure files visited by the calculator.
4. config file (str): A yaml file that takes specific inputs for the Gcbh calculaor. In the future, more functionalities will be read from the config file. Please check the example folders to check the currently available functionalities.
5. restart (bool): To control restart from previous calculations.
6. optimizer (ase.optimizer): An `ase Optimizer <https://wiki.fysik.dtu.dk/ase/ase/optimize.html>`_ that controls geometric optimization of a given ase.atoms object and reduce forces. The default is BFGS.
