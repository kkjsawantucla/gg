
# g(raph) -g(gcbh)

gg is an open-source code for building graph-based grand canonical basin hopping calculator

[![Documentation Status](https://readthedocs.org/projects/nequip/badge/?version=latest)]()

**PLEASE NOTE:** the code is currently under active development and is still in beta versions 0.x.x.

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

Alternatively, you can just add ./fga to your $PYTHONPATH. (not recommended)
~~~bash
export PYTHONPATH=$PYTHONPATH:"<path_to_gg>"
~~~

## Usage

#### Modifiers
The modifiers form the building block of the code. They determine how the atoms are modified during each basin hopping step. The code provides basic modifiers as building blocks for more complex modifiers.

1. Add Modifier
