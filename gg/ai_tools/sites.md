# Sites

The `Sites` class helps make graphs for modifiers to work on and also determines the surface site for modifications.

## FlexibleSites

This is a simple sites class which returns either atoms which aren't constrained as surface sites, or you can specify specific indices. Hardcoding of indices isn't advisable as the atoms object changes during GCBH runs.

```python
from gg.predefined_sites import FlexibleSites

FS = FlexibleSites(constraints=True, max_bond_ratio=1.2) # Define class to figure out surface

atoms = read('POSCAR')
list_sites = FS.get_sites(atoms)

```

### `class gg.sites.FlexibleSites(constraints: bool = False, index: list = None, max_bond_ratio: float | None = 1.2, max_bond: float | None = 0, contact_error: float | None = 0.2)`

**Bases:** `Sites`

**Args:**

* **`constraints`** (`bool`, optional): If true, only atoms which aren't constrained are considered. Defaults to `False`.
* **`index`** (`list`, optional): If a list of indices is given, it will be used as it is. Defaults to `None`.

#### `get_sites(atoms: Atoms) → list`

**Args:**

* **`atoms`** (`ase.Atoms`): Atoms object to determine sites.

**Returns:**

* **`list`**: List of atom indices considered for modifications.

---

## SurfaceSites

This class uses coordination number to determine the surface sites. However, we need to specify the maximum coordination allowed for each atom.

```python
from gg.predefined_sites import SurfaceSites

max_coord = {'Pt': 12, 'O': 4, 'H': 2} #Every atom in the POSCAR should be defined
SS = SurfaceSites(max_coord=max_coord, max_bond_ratio=1.2) # Define class to figure out surface

atoms = read('POSCAR')
list_sites = FS.get_sites(atoms)

```

### `class gg.sites.SurfaceSites(max_coord: dict, max_bond_ratio: float | None = 1.2, max_bond: float | None = 0, contact_error: float | None = 0.2, com: bool | None = True)`

**Bases:** `Sites`

**Args:**

* **`max_coord`** (`dict`): Dictionary of the maximum coordination of each element used. Only atoms with coordination less than this value will be considered.
* **`max_bond_ratio`** (`float`): Tolerance in the sum of covalent radii between two atoms to be considered a bond. Defaults to `1.2`.
* **`max_bond`** (`float`): Maximum distance of a bond allowed; ignored if equal to zero. Defaults to `0`.
* **`contact_error`** (`float`): To ensure atoms aren't too close to each other, the fraction of tolerance allowed. Defaults to `0.2`.
* **`com`** (`Optional[bool]`, optional): If true, atoms below the center of mass are ignored. Defaults to `True`.

#### `get_sites(atoms: Atoms, self_interaction: bool = False, bothways: bool = True) → list`

**Args:**

* **`atoms`** (`ase.Atoms`): Atoms object to determine sites.
* **`self_interaction`** (`bool`): Input of `ase.neighborlist`. Defaults to `True` (Note: Documentation text says "Defaults to True" in description but signature says `False`).
* **`bothways`** (`bool`): Input of `ase.neighborlist`. Defaults to `False` (Note: Documentation text says "Defaults to False" in description but signature says `True`).

**Returns:**

* **`list`**: List of atom indices considered for modifications.

## RuleSites

`RuleSites` is a composable `Sites` implementation that identifies candidate modification sites by running **one or more “rule” functions** (called `index_parsers`) and then **combining** their results via **union** or **intersection**. 

This is handy when “surface sites” aren’t captured by one heuristic alone—e.g., you might want atoms **above the center of mass** *and* **under-coordinated**.

### Example

```python
from gg.sites import RuleSites
from gg.sites import get_com_sites, get_surface_sites_by_coordination

max_coord = {"Al": 6, "Si": 6, "O": 4, "H": 2}

RS = RuleSites(
    index_parsers=[
        lambda atoms: get_com_sites(atoms, fraction=0.75, direction="above", axis="y"),
        lambda atoms: get_surface_sites_by_coordination(
            atoms, max_coord=max_coord, max_bond=2, max_bond_ratio=0.0
        ),
    ],
    combine_rules="intersection",
    max_bond=2,
    max_bond_ratio=0.0,
    contact_error=0.3,
)

atoms = read("POSCAR")
sites = RS.get_sites(atoms)          # returns a set of indices
sites_list = sorted(list(sites))     # optional: stable list form
```

---

### `class gg.sites.RuleSites(index_parsers: Callable | list[Callable] | None = None, combine_rules: str = "union", max_bond_ratio: float | None = 1.2, max_bond: float | None = 0, contact_error: float | None = 0.3)`

**Bases:** `Sites` 

**Args:**

* **`index_parsers`** (`Callable[[ase.Atoms], list] | list[Callable[[ase.Atoms], list]] | None`, optional):
  A single rule function or a list of rule functions. Each rule takes an `ase.Atoms` object and returns a list of atom indices.
  If `None`, a default rule is used that returns **all indices** (`range(len(atoms))`). 

* **`combine_rules`** (`str`, optional):
  How to combine the outputs of multiple rules:

  * `"union"`: include indices returned by **any** rule (default)
  * `"intersection"`: include only indices returned by **all** rules
    **Note:** the implementation actually uses only the **first character** of the string (`"u"` or `"i"`), so `"union"`/`"u"` and `"intersection"`/`"i"` behave equivalently. 

* **`max_bond_ratio`** (`float | None`, optional):
  Stored on the base `Sites` class for graph/bond construction helpers (e.g., for parsers that build graphs). Defaults to `1.2`. 

* **`max_bond`** (`float | None`, optional):
  Stored on the base `Sites` class as an optional absolute bond cutoff distance. If `0`, it is effectively ignored. Defaults to `0`. 

* **`contact_error`** (`float | None`, optional):
  Stored on the base `Sites` class as a tolerance for “too-close” contacts. Defaults to `0.3` in the code. 

---

#### `get_sites(atoms: Atoms) → set`

**Args:**

* **`atoms`** (`ase.Atoms`): Atoms object to analyze.

**Returns:**

* **`set`**: A set of atom indices representing the sites of interest after combining rule outputs.
  (If you need a list, use `sorted(list(result))`.) 

---

### Built-in rule helpers (optional to use)
The module also provides several ready-made rule functions you can plug into `index_parsers`, including: 
* `get_unconstrained_sites(atoms)` — indices not fixed by `ase.constraints.FixAtoms`
* `get_tagged_sites(atoms, tag=-1)` — indices whose `atoms[i].tag == tag`
* `get_com_sites(atoms, fraction=1.0, direction="above"|"below"|"both", axis="x"|"y"|"z"|0|1|2)` — indices selected relative to the center of mass along an axis
* `get_surface_sites_by_coordination(atoms, max_coord, ...)` — under-coordinated atoms based on a neighborlist/graph-derived coordination number

These can be mixed freely via `combine_rules` to build a robust “site definition” without hardcoding indices. 
