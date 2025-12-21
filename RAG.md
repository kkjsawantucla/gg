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

#### Defining surface sites by rule type
Use these patterns to tailor site selection to a new system:
- **Center-of-mass (COM) filters**: `FlexibleSites` supports COM fraction filtering to
  focus on atoms near the slab/cluster surface.
- **Atom-type filters**: pass element symbols (e.g., `surf_sym=["Cu", "Zn"]`) to
  modifiers and site classes to target specific surface species.
- **Coordination-based filters**: `SurfaceSites` selects atoms below a coordination
  threshold (via `max_coord`) to capture undercoordinated surface atoms.
- **Rule composition**: `RuleSites` can combine multiple rules (element, coordination,
  distance, tags) in a single selector; see `docs/source/examples.rst` for patterns.

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

## Minimal `run.py` scaffold (from scratch)
Use this as a complete starting point for a new system. Replace the calculator
block and any system-specific settings (surface symbols, adsorbate IDs, etc.).
```python
from ase.io import read
from gg.gcbh import Gcbh
from gg.predefined_sites import FlexibleSites
from gg.modifiers import Add, Remove

# 1) Load structure
atoms = read("POSCAR")  # or .xyz/.traj

# 2) Attach calculator (must be set before Gcbh)
# Example placeholder; replace with your calculator of choice
# from ase.calculators.emt import EMT
# atoms.calc = EMT()
atoms.calc = ...  # REQUIRED

# 3) Define surface sites
sites = FlexibleSites(constraints=True, max_bond_ratio=1.2)

# 4) Define modifiers (examples)
add_ads = Add(
    sites,
    "O",             # adsorbate element
    surf_coord=[2, 3],
    ads_id=["O"],
    surf_sym=["Cu"]  # surface element symbols
)
remove_ads = Remove(sites, "O", max_bond_ratio=1.2)

# 5) Create driver from YAML config and register modifiers
gcbh = Gcbh(atoms, config_file="input.yaml")
gcbh.add_modifier(add_ads, "Add O")
gcbh.add_modifier(remove_ads, "Remove O")

# 6) Run
gcbh.run(steps=1000)
```

**Required files and inputs**
- `input.yaml` must include `chemical_potential` and GCBH parameters.
- Structure file (`POSCAR`, `.xyz`, `.traj`) must match the calculator.
- Ensure `max_bond_ratio` and site selection parameters reflect your system.

**Common extensions**
- Use `gg.predefined_sites.SurfaceSites` for coordination-based selection.
- Use `gg.modifiers.AddBi` for bidentate adsorbates.
- Use `gg.gcbh.GcbhFlexOpt` for flexible optimization workflows.

## Calculator examples (MACE, FAIR-Chem, NequIP)
Attach any ASE-compatible calculator to `atoms.calc` before creating `Gcbh`.

**MACE**
```python
from mace.calculators import mace_mp

# MACE-MP default (downloads model weights if not present)
atoms.calc = mace_mp(model="large", dispersion=True, device="cuda")
```

**FAIR-Chem (Open Catalyst Project models)**
```python
from fairchem.core.models.model_registry import model_name_to_pretrained
from fairchem.core.models.ocp import OCPCalculator

model = model_name_to_pretrained["oc20"]
atoms.calc = OCPCalculator(checkpoint_path=model, cpu=False)
```

**NequIP**
```python
from nequip.ase import NequIPCalculator

atoms.calc = NequIPCalculator.from_deployed_model(
    "path/to/deployed_model.pth",
    device="cuda"
)
```

Notes:
- Exact import paths and model names can vary by release; consult each package’s docs.
- Ensure the calculator supports the elements in your system and any required device/cutoff settings.

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

## Updating this RAG on GitHub
- The RAG file lives at the repo root: `RAG.md`.
- To edit via GitHub’s web UI, open the file in the repository and use the
  “Edit” (pencil) icon, then commit the change in a branch/PR.
- Direct URL pattern: `https://github.com/<org>/<repo>/blob/<branch>/RAG.md`
  (replace org/repo/branch). For this repository, use the main branch.

## Pointers for LLMs
- **Use ASE `Atoms`** as the core object passed through gg APIs.
- **Always attach a calculator** before creating `Gcbh`.
- **Use `FlexibleSites` or `SurfaceSites`** to generate surface atoms rather than hardcoding indices.
- **Modifiers return lists** of candidate structures; gg runs accept/reject via basin hopping.
- For uniqueness and graph-based filtering, see `gg/utils_graph.py` and modifier options like `unique_method`.
