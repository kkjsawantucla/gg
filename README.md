
# g(raph) -g(gcbh)

gg is an open-source code for building graph-based grand canonical basin hopping calculator

[![Documentation Status](https://readthedocs.org/projects/nequip/badge/?version=latest)]()

**PLEASE NOTE:** The code is currently under active development and is still in beta versions 0. x.x.

### Requirements
- [Python](https://www.python.org/) 3.7  or later
- [NumPy](https://numpy.org/doc/stable/reference/)
- [ase](https://wiki.fysik.dtu.dk/ase/)
- [NetworkX](https://networkx.org/)
- [pandas](https://pandas.pydata.org/)
- [yaml](https://pyyaml.org/)

## Installation
Clone Directory
~~~bash
git clone https://github.com/kkjsawantucla/gg.git
~~~

Install using pip
~~~bash
cd gg/
pip install .
~~~

Alternatively, you can add ./fga to your $PYTHONPATH. (not recommended)
~~~bash
export PYTHONPATH=$PYTHONPATH:"<path_to_gg>"
~~~

## Usage

#### GCBH
The Gcbh calculator is an ase Dynamics Child on unhealthy steriods. It runs the grand canonical basin hopping, however certain functionalities are hard coded.

##### Inputs
1. atoms (ase.Atoms): An [ase Atoms](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) object as a starting point. The object should have a [ase.calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html) attached. 
2. logfile (str): path to a file that logs the calculator's output.
3. trajectory (str): path to a file that logs all the atoms structure files visited by the calculator.
4. config file (str): A yaml file that takes specific inputs for the Gcbh calculaor. In the future, more functionalities will be read from the config file. Please check the example folders to check the currently available functionalities.
5. restart (bool): To control restart from previous calculations.
6. optimizer (ase.optimizer): Optimizer that controls geometric optimization of a given ase.atoms object. The default is BFGS.

~~~bash
from gg.gcbh import Gcbh
from ase.io import read
from ase.calculators.emt import EMT

atoms = read('POSCAR')
atoms.calc = EMT()
G = Gcbh(atoms,config_file='input.yaml')
~~~

#### Modifiers
The modifiers form the building block of the code. They determine how the atoms are modified during each basin hopping step. The code provides basic modifiers as building blocks for more complex modifiers.

1. Add Modifier
The modifier can add an adsorbate, or moiety at specific sites on the parent atoms object. 
