
# g(raph) - g(cbh)

gg is an open-source code for building graph-based grand canonical basin hopping calculator

[![Documentation Status](https://readthedocs.org/projects/graph-gcbh/badge/?version=latest)](https://graph-gcbh.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

**PLEASE NOTE:** The code is currently under active development and is still in beta versions 0.x.x

## Requirements
- [python](https://www.python.org/) 3.7  or later
- [numPy](https://numpy.org/doc/stable/reference/)
- [ase](https://wiki.fysik.dtu.dk/ase/)
- [networkx](https://networkx.org/)
- [pandas](https://pandas.pydata.org/)
- [pyYAML](https://pyyaml.org/)
- [scipy](https://scipy.org/)

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

Alternatively, you can add ./gg to your $PYTHONPATH. (not recommended)
~~~bash
export PYTHONPATH=$PYTHONPATH:"<path_to_gg>"
~~~

## Usage
For usage check https://graph-gcbh.readthedocs.io/en/latest/index.html

## Calculator examples (MACE, FAIR chemistry, NequIP)
gg works with any ASE calculator. Below are short examples showing how to attach
popular ML potentials to an `ase.Atoms` object before passing it to `gg`.

### MACE
~~~python
from ase.io import read
from mace.calculators import mace_mp

atoms = read("POSCAR")
calc = mace_mp(model="medium-mpa-0", device="cuda")  # or device="cpu"
atoms.calc = calc
~~~

### Meta FAIR chemistry (Open Catalyst / fairchem)
~~~python
from ase.io import read
from fairchem.core.common.relaxation.ase_utils import OCPCalculator

atoms = read("POSCAR")
calc = OCPCalculator(
    model_name="eqv2_31M_omol",  # replace with your model name
    checkpoint_path="path/to/checkpoint.pt",  # optional if model is cached
    cpu=False,
)
atoms.calc = calc
~~~
Refer to the FAIR chemistry documentation for available model names and
checkpoint handling.

### NequIP
~~~python
from ase.io import read
from nequip.ase import NequIPCalculator

atoms = read("POSCAR")
calc = NequIPCalculator.from_deployed_model("path/to/nequip_model.pth")
atoms.calc = calc
~~~

## Cite
K. J. Sawant, P. Sautet, 2025, DOI 10.26434/chemrxiv-2025-t71bx
