# Modifiers

The modifiers form the building block of the code. They determine how the atoms are modified during each basin hopping step. The code provides basic modifiers as building blocks for more complex modifiers.

## Add Monodentate

The modifier can add a monodentate adsorbate, or moiety at specific sites on the parent atoms object.

```python
from gg.modifiers import Add
from gg.predefined_sites import FlexibleSites
from ase.build import fcc111

atoms = fcc111("Pt", size=(3, 3 , 4), vacuum=10.0)

FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
add_OH = Add(FS, ads="OH", surf_coord=[1,2,3], ads_id=["O"], surf_sym=["Pt"], print_movie=True, unique=True)

modified_atoms = add_OH.get_modified_atoms(atoms) #Atoms with the adsorbate

```

### `gg.modifiers.add.Add`

**Bases:** `ParentModifier`

**Args:**

* `surface_sites` (**gg.Sites**): Class that figures out surface sites.
* `ads` (**str** or **ase.Atoms**): Adsorbate to add.
* `surf_coord` (**list[int]**): How many bonds the adsorbing atom will make with the surface.
* `surf_sym` (**list[str]**): Surface elements where adsorbate can add.
* `ads_id` (**list[float]**): Strings denoting chemical symbol of adsorbate atom. Defaults to `None`.
* `ads_dist` (**str**, optional): Distance of adsorbate from surface site. If `ads_id` is mentioned, this variable is ignored. Defaults to `1.8`.
* `print_movie` (**bool**, optional): Return a movie of all sites or one random site. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique sites. Defaults to `True`.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.

#### `Add.get_modified_atoms(atoms: Atoms) → Atoms`

**Args:**

* `atoms` (**ase.Atoms**): The atoms object on which the adsorbate will be added.

**Returns:**

* `ase.Atoms` if `print_movie = True`
* `list[ase.Atoms]` if `print_movie = False`

---

## Add Bidentate

The modifier can add a bidentate adsorbate, or moiety at specific sites on the parent atoms object.

```python
#Example for Adding Formate as Bidentate Molecule
from gg.modifiers import AddBi
from gg.predefined_sites import FlexibleSites
from ase.build import fcc111

atoms = fcc111("Pt", size=(3, 3 , 4), vacuum=10.0)

FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
add_formate = AddBi(FS, 
    ads="HCOO", #Read adsorbate from database
    surf_coord=[1,2,3], #Co-ordination of surface sites
    ads_id=["O"], #The symbol in adsorbate which can attach to surface
    surf_sym=["Pt"], 
    print_movie=True, 
    unique=True)

modified_atoms = add_formate.get_modified_atoms(atoms) #Atoms with the adsorbate

```

### `gg.modifiers.add.AddBi`

**Bases:** `Add`

**Args:**

* `surface_sites` (**gg.Sites**): Class that figures out surface sites.
* `ads` (**str** or **ase.Atoms**): Adsorbate to add.
* `surf_coord` (**list[int]**): How many bonds the adsorbing atom will make with the surface.
* `surf_sym` (**list[str]**): Surface elements where adsorbate can add.
* `ads_id` (**list** of **[int or str]**): Strings denoting chemical symbol of adsorbate atom.
* `ads_dist` (**list** of **[float or str]**, optional): Distance of adsorbate from surface site. If it is a string denoting the chemical symbol of an adsorbate atom, then distance is set by atomic radii. Defaults to covalent radii of atoms mentioned in `ads_id`.
* `print_movie` (**bool**, optional): Return a movie of all sites or one random site. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique sites. Defaults to `True`.
* `ads_rotate` (**bool**, optional): Rotate atoms such that they point in +z direction. Defaults to `True`.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.

#### `AddBi.get_modified_atoms(atoms: Atoms) → Atoms`

**Args:**

* `atoms` (**ase.Atoms**): The atoms object on which the adsorbate will be added.

**Returns:**

* `ase.Atoms` if `print_movie = True`
* `list[ase.Atoms]` if `print_movie = False`

---

## Remove Adsorbate

Remove an adsorbate from the surface.

```python
#Example for removing OH from the surface
from gg.modifiers import Remove
from gg.predefined_sites import FlexibleSites

FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
remove_OH = Remove(FS, to_del="OH", print_movie=True)

atoms = read("POSCAR_with_OH") #The atoms object that has OHs to be removed
modified_atoms = remove_OH.get_modified_atoms(atoms) #Atoms with the adsorbate

```

### `gg.modifiers.modifiers.Remove`

**Bases:** `ParentModifier`

**Args:**

* `surface_sites` (**gg.Sites**): Class that figures out surface sites.
* `to_del` (**str** or **ase.Atoms**): Atoms to delete. If a string is provided, it tries to make a molecule out of it.
* `max_bond_ratio` (**float**, optional): Max bond ratio to make graph of atoms to delete. Defaults to `1.2`.
* `max_bond` (**int**, optional): Max bond ratio to make graph of atoms to delete. Defaults to `0`.
* `print_movie` (**bool**, optional): Return a movie of all sites or one random site. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique sites.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.

#### `Remove.get_modified_atoms(atoms: Atoms) → Atoms`

**Args:**

* `atoms` (**ase.Atoms**): The atoms object on which the adsorbate will be added.

**Returns:**

* `ase.Atoms` if `print_movie = True`
* `list[ase.Atoms]` if `print_movie = False`

---

## Replace Atoms/Molecules

Remove an adsorbate from the surface and replace it.

```python
#Example for replacing surface OH with surface NO
from gg.modifiers import Replace
from gg.predefined_sites import FlexibleSites
from ase.io import read

FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface sites
remove_OH = Replace(FS, to_del="OH", with_replace="NO", print_movie=True)

atoms = read("POSCAR_with_OH") #The atoms object that has OHs to be replaced
modified_atoms = remove_OH.get_modified_atoms(atoms) #Atoms with the new adsorbate

```

### `gg.modifiers.modifiers.Replace`

**Bases:** `Remove`

**Args:**

* `surface_sites` (**gg.Sites**): Class that figures out surface sites.
* `to_del` (**str** or **ase.Atoms**): Atoms to delete. If a string is provided, it tries to make a molecule out of it.
* `with_replace` (**str** or **ase.Atoms**): Atoms to replace with. If a string is provided, it tries to make a molecule out of it.
* `max_bond_ratio` (**float**, optional): Max bond ratio to make graph of atoms to delete. Defaults to `1.2`.
* `max_bond` (**int**, optional): Max bond ratio to make graph of atoms to delete. Defaults to `0`.
* `print_movie` (**bool**, optional): Return a movie of all sites or one random site. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique sites.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.

#### `Replace.get_modified_atoms(atoms: Atoms) → Atoms`

**Args:**

* `atoms` (**ase.Atoms**): The atoms object on which the adsorbate will be added.

**Returns:**

* `ase.Atoms` if `print_movie = True`
* `list[ase.Atoms]` if `print_movie = False`

---

## Swap Atoms

Swap two atoms on the surface.

```python
from gg.modifiers import Swap
from gg.predefined_sites import FlexibleSites
from ase.io import read

FS = FlexibleSites(constraints=True,max_bond_ratio=1.2) #Define class to figure out surface
swap = Swap(FS, swap_sym=["Pt","Au"], print_movie=True)

atoms = read('POSCAR_PtAu') #The atoms object with Pt and Au that can be swapped
modified_atoms = swap.get_modified_atoms(atoms) #Atoms with the adsorbate

```

### `gg.modifiers.modifiers.Swap`

**Bases:** `ParentModifier`

**Args:**

* `surface_sites` (**gg.Sites**): Class which figures out surface sites.
* `swap_sym` (**list**): List of atom symbols that are allowed to swap.
* `swap_ind` (**list**): List of indices to swap. Default to `None`.
* `print_movie` (**bool**, optional): Return a movie of all sites or one random site. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique sites. Defaults to `True`.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.

---

## Cluster modifiers

Cluster modifiers operate on free clusters rather than slabs. They are useful for exploring alternative cluster orientations or positions.

### `gg.modifiers.cluster.ClusterRotate`

Rotates a cluster of atoms around a surface normal determined from tagged sites.
Make sure to tag the cluster atoms using atom.tag = -1

### Example of rotating Pt cluster on TiO2

```python
from gg.modifiers import ClusterRotate
from gg.predefined_sites import FlexibleSites
from ase.io import read

atoms = read("POSCAR_Pt_TiO2")
for a in atoms:
  if a.symbol=='Pt':
    a.tag=-1
FS = FlexibleSites(tag=-1)
rotate = ClusterRotate(FS,contact_error=0.2,rotate_vector=(0,0,1))
modified_atoms = rotate.get_modified_atoms(atoms)
```


**Args:**
* `surface_sites` (**gg.Sites**): Class that figures out surface sites.
* `max_angle` (**float**, optional): Maximum rotation angle allowed (degrees). Defaults to `180`.
* `rotate_vector` (**tuple[float, float, float] | None**, optional): Vector along which atoms are rotated. If `None`, the surface normal is used. Defaults to `None`.
* `contact_error` (**float**, optional): Allowable tolerance in atoms touching. Defaults to `0.2`.
* `nmovie` (**int**, optional): Number of rotation attempts to generate. Defaults to `1`.
* `print_movie` (**bool**, optional): Return a movie of all rotations or one random rotation. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique structures. Defaults to `True`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.

### `gg.modifiers.cluster.ClusterTranslate`

Translates a cluster by a random displacement vector subject to allowed directions.
Make sure to tag the cluster atoms using atom.tag = -1

### Example of translating Pt cluster on TiO2

```python
from gg.modifiers import ClusterTranslate
from gg.predefined_sites import FlexibleSites
from ase.io import read

atoms = read("POSCAR_Pt_TiO2")
for a in atoms:
  if a.symbol=='Pt':
    a.tag=-1
FS = FlexibleSites(tag=-1)
trans = ClusterTranslate(FS,contact_error=0.2)
modified_atoms = trans.get_modified_atoms(atoms)
```


**Args:**
* `surface_sites` (**gg.Sites**): Class that figures out surface sites.
* `max_displace` (**float**, optional): Maximum displacement allowed (angstrom). Defaults to `5`.
* `allowed_direction` (**tuple[bool, bool, bool]**, optional): Allow displacement in x/y/z. Defaults to `(True, True, False)`.
* `contact_error` (**float**, optional): Allowable tolerance in atoms touching. Defaults to `0.2`.
* `nmovie` (**int**, optional): Number of translation attempts to generate. Defaults to `1`.
* `print_movie` (**bool**, optional): Return a movie of all translations or one random translation. Defaults to `False`.
* `unique` (**bool**, optional): Return only unique structures. Defaults to `True`.
* `unique_method` (**str | list[str]**, optional): How uniqueness is computed (e.g., `"fullgraph"` or `["C"]`). Defaults to `"fullgraph"`.
* `unique_depth` (**int**, optional): Subgraph depth for uniqueness when `unique_method` is not `"fullgraph"`. Defaults to `3`.
* `weight` (**float**): Weight for gcbh. Defaults to `1`.

#### `Swap.get_modified_atoms(atoms: Atoms) → Atoms`

**Args:**

* `atoms` (**ase.Atoms**): The atoms object on which the adsorbate will be added.

**Returns:**

* `ase.Atoms` if `print_movie = True`
* `list[ase.Atoms]` if `print_movie = False`

## ModifierAdder

Combine multiple modifiers into a single **sequential** modifier. This is useful when you want a *single move* in GCBH to perform a multi-step surface event (e.g., dissociative adsorption: add OH **then** add H).

### Example: dissociative H₂O adsorption (OH* + H*)

```python
from gg.modifiers import Add, ModifierAdder
from gg.predefined_sites import FlexibleSites
from ase.build import fcc111

atoms = fcc111("Pt", size=(3, 3, 4), vacuum=10.0)

FS = FlexibleSites(max_bond_ratio=1.2, com=0.5)

add_H  = Add(FS, ads="H",  surf_coord=[1], ads_id=["H"], surf_sym=["Pt","O"], print_movie=True)
add_OH = Add(FS, ads="OH", surf_coord=[1], ads_id=["O"], surf_sym=["Pt"], print_movie=True)

# Apply add_OH first, then add_H
add_H2O = ModifierAdder([add_OH, add_H], print_movie=True, unique=True)

modified_atoms = add_H2O.get_modified_atoms(atoms)
```

---

### API

```python
class gg.modifiers.modifiers.ModifierAdder(
    modifier_instances: list[ParentModifier],
    max_bond_ratio: float = 1.2,
    max_bond: float = 0,
    print_movie: bool = False,
    unique: bool = True,
    unique_method: str = "fullgraph",
    unique_depth: int = 3,
    weight: float = 1,
)
```

**Bases:** `ParentModifier`

#### Args

* `modifier_instances` (`list[gg.ParentModifier]`): **Required.** List of modifier instances to apply sequentially.
* `max_bond_ratio` (`float`, optional): Max bond ratio used to build the atom graph (e.g., for deletion/neighbor logic). Default `1.2`.
* `max_bond` (`int`, optional): Max bond cutoff override. Default `0`.
* `print_movie` (`bool`, optional): Save/return a movie of all sites or one random site. Default `False`.
* `unique` (`bool`, optional): Return only unique modified structures. Default `True`.
* `unique_method` (`str | list[str]`): How uniqueness is computed. Can be `"fullgraph"` or a list like `["C"]` to build subgraphs around selected species. Default `"fullgraph"`.
* `unique_depth` (`int`): Subgraph depth for uniqueness (ignored if `unique_method="fullgraph"`). Default `3`.
* `weight` (`float`): Weight for GCBH move selection. Default `1`.

#### Method

```python
get_modified_atoms(atoms: Atoms) -> Atoms
```

**Returns:** `ase.Atoms`

## Available Modifiers for Adsorbates

### String → `ase.Atoms` conversion logic

```python
if isinstance(ads, str):
    if ads in adsorbates:
        _ads = adsorbates[ads]
    elif ads in ase.collections.g2.names:
        _ads = ase.build.molecule(ads)
    elif ads in chemical_symbols:
        _ads = ase.Atoms(ads, positions=[(0, 0, 0)])
    else:
        raise RuntimeError(f"Cannot convert string to Formula {ads}")
```

### Supported adsorbate names (`adsorbate_names`)

```python
adsorbate_names = {
    "hydrogen": adsorbates["H2"],
    "oxygen": adsorbates["O2"],
    "water": adsorbates["H2O"],
    "hydroxyl": adsorbates["OH"],
    "hydroperoxyl": adsorbates["OOH"],
    "carbon monoxide": adsorbates["CO"],
    "carbon dioxide": adsorbates["CO2"],
    "nitric oxide": adsorbates["NO"],
    "ammonia": adsorbates["NH3"],
    "amino": adsorbates["NH2"],
    "imidogen": adsorbates["NH2"],
    "methyl": adsorbates["CH3"],
    "methylene": adsorbates["CH2"],
    "methylidyne": adsorbates["CH"],
    "formate": adsorbates["HCOO"],
}
```
