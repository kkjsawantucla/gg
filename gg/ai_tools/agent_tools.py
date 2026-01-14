import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Tuple

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

VECTOR_STORE_ID = "vs_6963443265b481919b82d51c2841850f"

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
    include_retrieval: bool = True,
) -> str:
    """Generate gg-based Python code from a natural-language prompt."""
    kwargs: Dict[str, Any] = {
        "model": model,
        "instructions": instructions,
        "input": user_prompt,
        "tools": [
            {
                "type": "file_search",
                "vector_store_ids": [VECTOR_STORE_ID],
            }
        ],
    }
    if include_retrieval:
        # Helpful for debugging what got retrieved
        kwargs["include"] = ["file_search_call.results"]

    response = client.responses.create(**kwargs)
    return response.output_text


def dry_run(python_code: str, *, timeout_s: int = 30) -> Dict[str, Any]:
    """
    Runs python_code in a fresh interpreter process.

    Returns a dict with:
      ok (bool), returncode (int), stdout (str), stderr (str)
    """
    with tempfile.TemporaryDirectory() as td:
        run_path = Path(td) / "run.py"
        run_path.write_text(python_code, encoding="utf-8")

        proc = subprocess.run(
            ["python", str(run_path)],
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )

        return {
            "ok": proc.returncode == 0,
            "returncode": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }


def _compose_repair_prompt(
    *,
    original_request: str,
    previous_prompt: str,
    previous_code: str,
    run_result: Dict[str, Any],
    iter_idx: int,
) -> str:
    stdout = (run_result.get("stdout") or "").strip()
    stderr = (run_result.get("stderr") or "").strip()
    rc = run_result.get("returncode")

    # Keep the prompt reasonably sized; still include enough context to fix the issue.
    max_code_chars = 10000
    code_snippet = previous_code if len(previous_code) <= max_code_chars else previous_code[:max_code_chars] + "\n# ... (truncated) ...\n"

    return (
        "ORIGINAL REQUEST:\n"
        f"{original_request}\n\n"
        f"ITERATION: {iter_idx}\n\n"
        "PREVIOUS PROMPT USED (for context):\n"
        f"{previous_prompt}\n\n"
        "THE PREVIOUSLY GENERATED CODE (that failed):\n"
        "```\n"
        f"{code_snippet}\n"
        "```\n\n"
        "RUNTIME RESULT WHEN EXECUTED (python run.py):\n"
        f"Return code: {rc}\n\n"
        "STDOUT:\n"
        f"{stdout if stdout else '(empty)'}\n\n"
        "STDERR / TRACEBACK:\n"
        f"{stderr if stderr else '(empty)'}\n\n"
        "TASK:\n"
        "Fix the code so it runs successfully and still satisfies the ORIGINAL REQUEST. "
        "Return ONLY corrected Python code."
    )


def iterative_generate_gg_code(
    user_prompt: str,
    *,
    max_iters: int = 3,
    model: str = "gpt-4.1-mini",
    timeout_s: int = 30,
    include_retrieval: bool = True,
    run_check: bool = True,
) -> Tuple[str, Dict[str, Any]]:
    """
    Loop through: prompt -> code -> error -> new_prompt -> new_code

    - max_iters is user-defined.
    - If run_check=False, the function will generate once and skip execution.

    Returns (final_code, metadata).
    """
    if max_iters < 1:
        raise ValueError("max_iters must be >= 1")

    prompt = user_prompt
    last_code = ""
    last_run: Dict[str, Any] = {}

    for i in range(1, max_iters + 1):
        instructions = SYSTEM_CODER if i == 1 else SYSTEM_REPAIR
        code_temp = generate_gg_code(
            prompt,
            model=model,
            instructions=instructions,
            include_retrieval=include_retrieval,
        )
        last_code = code_temp

        if not run_check:
            return code_temp, {
                "ok": True,
                "iterations": i,
                "note": "run_check=False (skipped execution).",
                "last_run": {},
            }

        last_run = dry_run(code_temp, timeout_s=timeout_s)
        if last_run.get("ok"):
            return code_temp, {
                "ok": True,
                "iterations": i,
                "last_run": last_run,
            }

        # Build the next prompt using the error
        prompt = _compose_repair_prompt(
            original_request=user_prompt,
            previous_prompt=prompt,
            previous_code=code_temp,
            run_result=last_run,
            iter_idx=i,
        )

    return last_code, {
        "ok": False,
        "iterations": max_iters,
        "last_run": last_run,
        "note": "Reached max_iters without a successful run. Consider increasing max_iters or improving the initial prompt.",
    }
