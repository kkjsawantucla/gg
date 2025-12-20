# gg (graph-gcbh) RAG Reference

Use this file as retrieval context when prompting an LLM to work with the **gg** package. It summarizes what the library does, its key APIs, and where to look for examples.

## What gg is
- **gg** builds graph-based *grand canonical basin hopping* (GCBH) simulations for surface/cluster chemistry using **ASE** `Atoms` objects.
- The main simulation driver is `gg.gcbh.Gcbh`, which manages basin hopping, logging, and optimization steps.
- Structure modifications (adsorption, removal, swap, replace, cluster moves) are handled by **modifiers** in `gg/modifiers/`.
- Surface site identification is handled by **sites** classes in `gg/predefined_sites.py` and `gg/sites.py`.

## Installation (quick)
```bash
git clone https://github.com/kkjsawantucla/gg.git
cd gg
pip install .
```

## Documentation (GitHub Pages / ReadTheDocs)
- Documentation is hosted at **https://graph-gcbh.readthedocs.io/en/latest/**.
- The docs source lives in `docs/source/` and is built into `docs/build/`.
- For RAG context, prioritize `docs/source/*.rst` and the online docs above.

## Core entry points
### Simulation driver
- **`gg.gcbh.Gcbh`**: main basin hopping engine.
  - Requires an `ase.Atoms` with a calculator attached.
  - Reads YAML config via `config_file=...`.
  - Use `add_modifier(...)` to register structure modifications.
  - Run with `G.run(steps=...)`.
- **`gg.gcbh.GcbhFlexOpt`**: variant with flexible optimization (see `docs/source/GCBH.rst`).

### Sites (surface selection)
Common choices for identifying surface sites:
- **`gg.predefined_sites.FlexibleSites`**
  - Flexible selection (constraints, tags, COM fraction).
- **`gg.predefined_sites.SurfaceSites`**
  - Uses coordination numbers with `max_coord`.
- **`gg.sites.RuleSites`**
  - Combine multiple rules (see `docs/source/examples.rst`).

### Modifiers (structure changes)
Modifiers operate on `ase.Atoms` and return a list of modified structures:
- `gg.modifiers.add.Add` (monodentate adsorbate)
- `gg.modifiers.add.AddBi` (bidentate adsorbate)
- `gg.modifiers.modifiers.Remove`
- `gg.modifiers.modifiers.Replace`
- `gg.modifiers.modifiers.Swap`
- `gg.modifiers.cluster.ClusterRotate`
- `gg.modifiers.cluster.ClusterTranslate`
- `gg.modifiers.modifiers.ModifierAdder` (compose modifiers)

Each modifier exposes `get_modified_atoms(atoms)`.

## Typical workflow (Python)
```python
from ase.io import read
from gg.gcbh import Gcbh
from gg.predefined_sites import FlexibleSites
from gg.modifiers import Add, Remove

atoms = read("POSCAR")
atoms.calc = ...  # attach ASE calculator

sites = FlexibleSites(constraints=True, max_bond_ratio=1.2)
add_O = Add(sites, "O", surf_coord=[2,3], ads_id=["O"], surf_sym=["Cu", "Zn"])
rem_O = Remove(sites, "O", max_bond_ratio=1.2)

G = Gcbh(atoms, config_file="input.yaml")
G.add_modifier(add_O, "Add O")
G.add_modifier(rem_O, "Remove O")
G.run(steps=1000)
```

## End-to-end `run.py` template (for a new system)
Use this as a complete, minimal starting point. Replace placeholders with your system-specific inputs.

```python
# run.py
from ase.io import read
from gg.gcbh import Gcbh
from gg.predefined_sites import FlexibleSites, SurfaceSites
from gg.modifiers import Add, Remove, Replace, Swap, ModifierAdder, ClusterRotate, ClusterTranslate

# 1) Load structure and attach a calculator
atoms = read("POSCAR")  # or .xyz/.traj

# Example: user-provided ASE calculator
# from mace.calculators import mace_mp
# atoms.calc = mace_mp(model="medium-mpa-0", default_dtype="float64", device="cuda")
atoms.calc = ...  # REQUIRED: attach an ASE calculator

# 2) Define surface/site selection
# Option A: flexible sites (use constraints/tags/COM)
surface_sites = FlexibleSites(constraints=True, max_bond_ratio=1.2)

# Option B: coordination-based sites
# max_coord = {"Pt": 12, "O": 4, "H": 2}
# surface_sites = SurfaceSites(max_coord=max_coord, max_bond_ratio=1.2)

# 3) Define modifiers (examples below; choose what applies)
add_O = Add(surface_sites, "O", surf_coord=[2, 3], ads_id=["O"], surf_sym=["Pt"])
rem_O = Remove(surface_sites, "O", max_bond_ratio=1.2)

# Example: composite modifier (swap-like)
swap_O = ModifierAdder([add_O, rem_O])

# Example: cluster moves (tag cluster atoms with atom.tag = -1)
# for a in atoms:
#     if a.symbol == "Pt":
#         a.tag = -1
# cluster_sites = FlexibleSites(tag=-1)
# cl_rot = ClusterRotate(cluster_sites, contact_error=0.25, rotate_vector=(0, 0, 1))
# cl_trans = ClusterTranslate(cluster_sites, contact_error=0.2)

# 4) Initialize GCBH with YAML config
G = Gcbh(atoms, config_file="input.yaml")

# 5) Register modifiers (weights can be controlled in YAML via mod_weights)
G.add_modifier(add_O, "Add O")
G.add_modifier(rem_O, "Remove O")
G.add_modifier(swap_O, "Swap O")
# G.add_modifier(cl_rot, "Cluster Rotate")
# G.add_modifier(cl_trans, "Cluster Translate")

# 6) Optional: remove gas species generated during the run
# G.add_delete_gas(gas_species=["H2"])

# 7) Run
G.run(steps=1000)
```

## Minimal `input.yaml` for a new system
`chemical_potential` is required. These values control the grand canonical reference.

```yaml
temp: 550
stop_steps: 2000
stop_opt: 500

chemical_potential:
  Pt: -5.0
  O: -7.5
  H: -3.4

check_graphs: true
graph_method: fullgraph
max_bond_ratio: 1.2
initialize: true

# Optional tuning
# max_bond: 2.0
# max_history: 25
# area: false
# vib_correction: false
# detect_gas: null
# mod_weights:
#   Add O: 1.0
#   Remove O: 1.0
```

### Notes for creating `run.py` from scratch
- **Calculator is mandatory**: `Gcbh` raises if `atoms.calc` is `None`.
- **chemical_potential is mandatory** in the YAML config.
- If `initialize: true`, the initial structure is optimized before hopping.
- Use **`mod_weights`** to bias modifier selection.
- Use `restart=True` with `Gcbh(..., restart=True)` to resume from `current_status.pkl`.
- Modifier `unique_method` and `graph_method` control graph-based deduplication (see `gg/utils_graph.py`).

## Typical workflow (CLI)
The `bin/` scripts wrap simple Add modifiers:
- **`add_mono`**: monodentate adsorbate
- **`add_bi`**: bidentate adsorbate

Example:
```bash
add_mono -s POSCAR_Pt -a OH -sc 1 2 3 -aa O
```

## Configuration (YAML)
`Gcbh` reads a YAML config with fields like:
- `temp`, `stop_steps`, `stop_opt`
- `chemical_potential` (required)
- `check_graphs`, `max_bond_ratio`, `initialize`

Example in `docs/source/examples.rst` and `examples/`.

## Output files from a run
Typical outputs for `G.run(...)`:
- `local_minima.traj` (accepted structures)
- `gcbh.log` (run log)
- `gcbh.traj` (all structures)
- `current_status.pkl` (restart state)
- `opt_folder/opt_*` (optimization subfolders)

## Where to find examples
- `docs/source/examples.rst`
- `examples/`

## Pointers for LLMs
- **Use ASE `Atoms`** as the core object passed through gg APIs.
- **Always attach a calculator** before creating `Gcbh`.
- **Use `FlexibleSites` or `SurfaceSites`** to generate surface atoms rather than hardcoding indices.
- **Modifiers return lists** of candidate structures; gg runs accept/reject via basin hopping.
- For uniqueness and graph-based filtering, see `gg/utils_graph.py` and modifier options like `unique_method`.
