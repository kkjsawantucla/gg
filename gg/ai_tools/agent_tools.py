import os
from pathlib import Path
from typing import Any, Dict
from collections import Counter

import numpy as np
from ase.io import read
from openai import OpenAI

SYSTEM_CODER = """
You are a Python Code Generator for the 'gg' computational catalysis library.
Your goal is to translate user requests into executable Python scripts using the 'gg' library.

### Library Context (gg)
- **Sites**: Use `FlexibleSites`, `RuleSites`, or `SurfaceSites` to define where atoms go.
- **Modifiers**: Use `Add`, `Remove`, `Replace`, `Swap` to change the structure.
- **Data**: Use `ase.build` for surfaces (fcc111, etc.) unless a file is provided.

### Strict Output Rules
1. Return **ONLY** valid Python code. No markdown, no explanations, no `python` tags.
2. Always include necessary imports (e.g., `from gg.modifiers import ...`, `from ase.build import ...`).
3. Follow the user's style: define sites, define the modifier, and call `.get_modified_atoms()`.

### Few-Shot Examples

User: "Add OH to Pt(111) atoms at hollow sites"
Assistant:
from gg.modifiers import Add
from gg.predefined_sites import FlexibleSites

FS = FlexibleSites(constraints=True, max_bond_ratio=1.2) 
add_OH = Add(
    FS, 
    ads="OH",
    surf_coord=[3],
    ads_id=["O"],
    surf_sym=["Pt"],
    print_movie=True,
    unique=True
)

User: "Dissociatively add H2O to Pt(111) surface atoms at top sites"
Assistant:
from gg.modifiers import Add, ModifierAdder
from gg.predefined_sites import FlexibleSites

FS = FlexibleSites(max_bond_ratio=1.2,com=0.5)
add_H = Add(
    FS,
    ads="H",
    surf_coord=[1],
    ads_id=["H"],
    surf_sym=["Pt","O"]
    )
add_OH = Add(
    FS,
    ads="OH",
    surf_coord=[1],
    ads_id=["O"],
    surf_sym=["Pt"],print_movie=True)
add_H2O = ModifierAdder(
    [add_OH, add_H],
    print_movie=True,
    unique=True)
"""

SYSTEM_REPAIR = """
You are a Python Code Debugger for the 'gg' computational catalysis library.

You will be given:
- the user's original request
- a previous Python script you generated
- the runtime error (traceback / stderr) produced when running that script

Your job:
- Return ONLY corrected, executable Python code that satisfies the user's original request.
- Preserve the user's intent and use the 'gg' library idioms (sites, modifiers, etc.).
- Fix imports, missing variables, wrong API calls, and any runtime issues.
- Do NOT include markdown or explanations. Output must be only Python code.
"""

VECTOR_STORE_ID = "vs_6966fc692edc8191b14d5899006ca850"

# 1. Configuration
api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise ValueError("OPENAI_API_KEY must be set to use the OpenAI client.")
client = OpenAI(api_key=api_key)


def generate_gg_code(
    user_prompt: str,
    *,
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
    return response.output_text


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
    ]
    return "\n".join(lines)
