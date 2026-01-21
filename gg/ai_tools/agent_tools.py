import os
from pathlib import Path
from typing import Any, Dict
from collections import Counter

import numpy as np
from ase.io import read
from openai import OpenAI


SYSTEM_CODER = """
System: You are a Python code generator for the `gg` computational catalysis library. Produce a single, fully executable Python script that uses `gg` to build or modify atomic structures per the user's request.

## Core `gg` Concepts
- **Sites**: Use `FlexibleSites`, `RuleSites`, or `SurfaceSites` for adsorption or substitution sites as appropriate.
- **Modifiers**: Use modifiers like `Add`, `Remove`, `Replace`, `Swap`, `ClusterTranslate`, `ClusterRotate`, and combinators such as `ModifierAdder`.
- **Structures/Data**:
  - For surfaces: if no file is provided, construct with `ase.build` (e.g., `fcc111`, `fcc100`).
  - If no surface is specified, only infer a default if intent is clear; otherwise, select the minimal valid structure.
- Bidentate adsorbates: If the user requests bidentate binding (e.g., “bidentate,”, “bind through two atoms,” or provides two binding atoms), use AddBi (or the most specific bidentate-capable gg modifier) instead of Add.

## Script Requirements
- Output script must execute end-to-end as given.
- For adsorbates: if not found in `gg.data` or ASE, add a comment noting absence.
- Structure script in this order:
  1. Import statements
  2. Creation or loading of `atoms`
  3. Site definition
  4. Modifier(s) definition
  5. Modifier application via `.get_modified_atoms(atoms)`
  6. Assign output to a variable like `modified_atoms`
- For movies, uniqueness, or sampling: set `print_movie=True`, `unique=True`, or relevant flags as needed.
- Default to `FlexibleSites` unless user specifies otherwise.

## Strict Output Rules
1. Output only valid Python code—no Markdown or explanations.
2. Include all needed import statements (e.g., `from gg.modifiers import ...`, `from gg.predefined_sites import ...`, `from ase.build import ...`).
3. Match user order: define sites, then modifiers, then call `.get_modified_atoms()`.
4. Use exact parameter names/types as required by `gg` modifier signatures; do not invent them.
5. For structure files and modifications: ensure all referenced atom symbols exist in the structure.

## Import Information
1. `surf_coord` (with `Add`, `AddBi`): 1 = top, 2 = bridge, 3 = hollow site.
2. `gg` does not initially distinguish fcc/hcp hollow sites; uniqueness handled when checking for uniqueness.

## Example
User: "Dissociatively add H2O to Pt(111) surface atoms at top sites"

from gg.modifiers import Add, ModifierAdder
from gg.predefined_sites import FlexibleSites
from ase.build import fcc111

atoms = fcc111("Pt", size=(3,3,4), vacuum=10.0)

FS = FlexibleSites(max_bond_ratio=1.2, com=0.5)

add_H = Add(
    FS,
    ads="H",
    surf_coord=[1],
    ads_id=["H"],
    surf_sym=["Pt", "O"],
    print_movie=True
)

add_OH = Add(
    FS,
    ads="OH",
    surf_coord=[1],
    ads_id=["O"],
    surf_sym=["Pt"],
    print_movie=True
)

add_H2O = ModifierAdder(
    [add_OH, add_H],
    print_movie=True,
    unique=True
)

modified_atoms = add_H2O.get_modified_atoms(atoms)
"""


VECTOR_STORE_ID = "vs_696ff974f1a48191960a94c8b2ccef84"

# 1. Configuration
api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise ValueError("OPENAI_API_KEY must be set to use the OpenAI client.")
client = OpenAI(api_key=api_key)

def extract_final_text(response) -> str:
    """Python SDK returns typed objects, not dicts"""
    for item in reversed(response.output):
        if getattr(item, "type", None) == "message":
            for content in getattr(item, "content", []) or []:
                if getattr(content, "type", None) == "output_text":
                    return content.text
    return ""


def generate_gg_code(
    user_prompt: str,
    model: str = "gpt-4.1-mini",
    instructions: str = SYSTEM_CODER,
    structure_path: Path | None = None,
    include_retrieval: bool = True,
    use_search_tool: bool = True,
) -> str:
    """Generate gg-based Python code from a natural-language prompt."""
    if structure_path is not None:
        instructions = f"{instructions}\n\n{attach_structure_context(structure_path)}"
    kwargs: Dict[str, Any] = {
        "model": model,
        "instructions": instructions,
        "input": user_prompt,
    }
    if use_search_tool:
        kwargs["tools"] = [
            {
                "type": "file_search",
                "vector_store_ids": [VECTOR_STORE_ID],
            }
        ]
    if include_retrieval and use_search_tool:
        # Helpful for debugging what got retrieved
        kwargs["include"] = ["file_search_call.results"]

    response = client.responses.create(**kwargs)
    code = extract_final_text(response)
    assert code.strip(), "Model returned no code"
    return code


def _load_ai_tools_context() -> str:
    ai_tools_dir = Path(__file__).resolve().parent
    context_parts = []
    for doc_name in ("modifiers.md", "sites.md"):
        doc_path = ai_tools_dir / doc_name
        context_parts.append(f"## {doc_name}\n{doc_path.read_text(encoding='utf-8')}")
    return "\n\n".join(context_parts)


def generate_gg_code_with_local_docs(
    user_prompt: str,
    model: str = "gpt-4.1-mini",
    structure_path: Path | None = None,
) -> str:
    """Generate gg-based Python code using local docs for context."""
    local_context = _load_ai_tools_context()
    structure_context = ""
    if structure_path is not None:
        structure_context = attach_structure_context(structure_path)
    if structure_context:
        instructions = f"{SYSTEM_CODER}\n\nSYSTEM PLANNER CONTEXT:\n{local_context}\n\n{structure_context}"
    else:
        instructions = f"{SYSTEM_CODER}\n\nSYSTEM PLANNER CONTEXT:\n{local_context}"
    return generate_gg_code(
        user_prompt,
        model=model,
        instructions=instructions,
        include_retrieval=False,
        use_search_tool=False,
    )


def attach_structure_context(poscar_path: Path) -> str:
    """_summary_"""
    poscar_path = Path(poscar_path)
    if not poscar_path.exists():
        raise FileNotFoundError(f"Structure file not found: {poscar_path}")

    atoms = read(str(poscar_path))
    total_atoms = len(atoms)
    counts = Counter(atoms.get_chemical_symbols())
    composition = ", ".join(f"{element}:{counts[element]}" for element in sorted(counts))
    pbc = atoms.get_pbc()
    pbc_str = f"{bool(pbc[0])},{bool(pbc[1])},{bool(pbc[2])}"
    cell = atoms.get_cell().array
    cell_str = "; ".join(
        ",".join(f"{value:.3f}" for value in vector) for vector in cell
    )
    center_of_mass = atoms.get_center_of_mass()
    com_str = ",".join(f"{value:.3f}" for value in center_of_mass)
    z_positions = atoms.get_positions()[:, 2]
    z_min = float(np.min(z_positions))
    z_max = float(np.max(z_positions))
    z_range_str = f"{z_min:.3f},{z_max:.3f}"

    filename = poscar_path.name
    lines = [
        "Structure context (ASE-parsed):",
        f"Total atoms: {total_atoms}",
        f"Composition: {composition}",
        f"PBC (x,y,z): {pbc_str}",
        f"Cell vectors (A): {cell_str}",
        f"Center of mass (A): {com_str}",
        f"Z range (A): {z_range_str}",
        "Instructions:",
        "Structure already parsed with ASE.",
        f'Generated code MUST load structure using ase.io.read("{filename}").',
        "Do NOT build surfaces using ase.build.",
        "Apply all gg modifiers directly to the loaded atoms.",
        "If using SurfaceSites(max_coord=...), max_coord MUST be a dict covering ALL elements in atoms.",
    ]
    return "\n".join(lines)
