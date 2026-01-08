from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple
from openai import OpenAI
from plan_contract import plan_json_schema, validate_plan

# 1) Tool Specs (what you pass to the model)
TOOLS: List[Dict[str, Any]] = [
    {
        "type": "function",
        "function": {
            "name": "retrieve_gg_docs",
            "description": (
                "Search local gg docs/examples and return relevant API snippets "
                "(signatures/examples)."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query, e.g. 'AddBi formate bidentate'",
                    },
                    "top_k": {
                        "type": "integer",
                        "description": "How many snippets to return",
                        "default": 6,
                    },
                },
                "required": ["query"],
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "emit_plan",
            "description": (
                "Create a PlanJSON that conforms to the gg contract, given a "
                "natural-language request and doc snippets."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "nl_request": {"type": "string"},
                    "snippets": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "source": {"type": "string"},
                                "text": {"type": "string"},
                            },
                            "required": ["source", "text"],
                            "additionalProperties": False,
                        },
                    },
                },
                "required": ["nl_request", "snippets"],
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "compile_plan",
            "description": (
                "Compile a validated PlanJSON into deterministic python code "
                "that uses gg primitives."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "object",
                        "description": "PlanJSON object matching the contract",
                    },
                },
                "required": ["plan"],
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "dry_run",
            "description": (
                "Execute compiled python code in a subprocess and return pass/fail "
                "+ traceback."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "python_code": {"type": "string"},
                    "timeout_s": {"type": "integer", "default": 30},
                },
                "required": ["python_code"],
                "additionalProperties": False,
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "revise",
            "description": (
                "Revise a PlanJSON to fix validation/runtime errors, returning "
                "an updated PlanJSON."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {"type": "object"},
                    "errors": {"type": "string"},
                    "snippets": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "source": {"type": "string"},
                                "text": {"type": "string"},
                            },
                            "required": ["source", "text"],
                            "additionalProperties": False,
                        },
                    },
                },
                "required": ["plan", "errors", "snippets"],
                "additionalProperties": False,
            },
        },
    },
]

# 2) Retrieval (retrieve_gg_docs)
@dataclass
class Snippet:
    source: str
    text: str

def _load_text_files(root: Path, exts: Tuple[str, ...]) -> List[Tuple[str, str]]:
    items: List[Tuple[str, str]] = []
    for p in root.rglob("*"):
        if not (p.is_file() and p.suffix.lower() in exts):
            continue
        try:
            txt = p.read_text(encoding="utf-8", errors="ignore")
        except (OSError, UnicodeError):
            continue
        items.append((str(p), txt))
    return items


def _score(query: str, text: str) -> int:
    q = re.findall(r"[A-Za-z0-9_]+", query.lower())
    t = text.lower()
    return sum(t.count(w) for w in q)


def retrieve_gg_docs(query: str, top_k: int = 6) -> List[Dict[str, str]]:
    root = Path(os.environ.get("GG_DOC_ROOT", ".")).resolve()
    candidates: List[Tuple[str, str]] = []

    for sub, exts in [
        ("docs", (".md", ".rst", ".txt", ".py")),
        ("examples", (".py", ".md")),
        ("gg", (".py",)),
        ("tests", (".py",)),
        ("bin", (".py", ".sh", ".bash")),
    ]:
        d = root / sub
        if d.is_dir():
            candidates += _load_text_files(d, exts)

    for f in ("README.md", "RAG.md"):
        p = root / f
        if p.is_file():
            candidates.append((str(p), p.read_text(encoding="utf-8", errors="ignore")))

    scored = [(s, src, txt) for (src, txt) in candidates if (s := _score(query, txt)) > 0]
    scored.sort(key=lambda x: x[0], reverse=True)

    focus = (re.findall(r"[A-Za-z0-9_]+", query.lower()) or [""])[0]
    out: List[Snippet] = []
    for _, src, txt in scored[:top_k]:
        i = txt.lower().find(focus) if focus else 0
        i = 0 if i < 0 else i
        out.append(Snippet(src, txt[max(0, i - 600) : min(len(txt), i + 1000)].strip()))

    return [{"source": s.source, "text": s.text} for s in out]

# 3) emit_plan + revise (Structured Outputs via Responses API)
PLAN_JSON_SCHEMA: Dict[str, Any] = {
    "name": "gg_plan",
    "schema": plan_json_schema(),
}

SYSTEM_PLANNER = """You translate natural language catalyst-surface modification requests into a STRICT PlanJSON.
Rules:
- Output MUST conform exactly to the provided JSON schema.
- Only use supported gg primitives (Sites: FlexibleSites/SurfaceSites/RuleSites; Modifiers: Add/AddBi/Remove/Replace/Swap/Rattle/Translate/ClusterRotate/ClusterTranslate/ModifierAdder).
- Do not invent kwargs. Prefer kwargs visible in the snippets.
- Make step names valid python identifiers and unique.
"""

def emit_plan(
    client: OpenAI, model: str, nl_request: str, snippets: List[Dict[str, str]]
) -> Dict[str, Any]:
    user_payload = {
        "nl_request": nl_request,
        "snippets": snippets,
        "instructions": "Return a single PlanJSON that is valid and minimal (no extra steps).",
    }

    resp = client.responses.create(
        model=model,
        input=[
            {"role": "system", "content": SYSTEM_PLANNER},
            {"role": "user", "content": json.dumps(user_payload)},
        ],
        text={
            "format": {
                "type": "json_schema",
                "strict": True,
                "name": PLAN_JSON_SCHEMA.get("name", "gg_plan"),
                "schema": PLAN_JSON_SCHEMA["schema"],
            }
        },
    )
    out = resp.output_text
    return json.loads(out)


def revise(
    client: OpenAI,
    model: str,
    plan: Dict[str, Any],
    errors: str,
    snippets: List[Dict[str, str]],
) -> Dict[str, Any]:
    user_payload = {
        "current_plan": plan,
        "errors": errors,
        "snippets": snippets,
        "instructions": "Fix the plan so it validates and avoids the runtime error. Keep changes minimal.",
    }
    resp = client.responses.create(
        model=model,
        input=[
            {"role": "system", "content": SYSTEM_PLANNER},
            {"role": "user", "content": json.dumps(user_payload)},
        ],
        text={
            "format": {
                "type": "json_schema",
                "strict": True,
                "name": PLAN_JSON_SCHEMA.get("name", "gg_plan"),
                "schema": PLAN_JSON_SCHEMA["schema"],
            }
        },
    )
    return json.loads(resp.output_text)

# 4) compile_plan (deterministic codegen)
def _py(obj: Any) -> str:
    """Render python literal deterministically."""
    return (
        json.dumps(obj)
        .replace("true", "True")
        .replace("false", "False")
        .replace("null", "None")
    )

def _rule_parser_expr(rule: Dict[str, Any]) -> str:
    name = rule["name"]
    kwargs = rule.get("kwargs", {}) or {}
    kwargs_expr = _py(kwargs)
    return f"lambda atoms: {name}(atoms, **{kwargs_expr})"

def compile_plan(plan: Dict[str, Any]) -> str:
    """
    Deterministic codegen: plan -> python code using gg.
    Assumes plan matches your contract.
    """
    validate_plan(plan)

    sites = plan["sites"]
    steps = plan["steps"]

    lines: List[str] = []
    lines += [
        "from gg.predefined_sites import FlexibleSites, SurfaceSites",
        "from gg.sites import RuleSites, get_com_sites, get_surface_sites_by_coordination, get_tagged_sites, get_unconstrained_sites",
        "from gg.modifiers import Add, AddBi, Remove, Replace, Swap, Rattle, Translate, ClusterRotate, ClusterTranslate, ModifierAdder",
        "",
        "# NOTE: you must supply `atoms` (ASE Atoms) before running modifiers.",
        "# e.g., from ase.build import fcc111; atoms = fcc111('Pt', size=(3,3,4), vacuum=10.0)",
        "",
    ]

    stype = sites["type"]
    skw = sites.get("kwargs", {})
    if stype == "FlexibleSites":
        lines.append(f"FS = FlexibleSites(**{_py(skw)})")
        sites_var = "FS"
    elif stype == "SurfaceSites":
        lines.append(f"SS = SurfaceSites(**{_py(skw)})")
        sites_var = "SS"
    elif stype == "RuleSites":
        rule_parsers = [
            _rule_parser_expr(rule) for rule in skw.get("index_parsers", [])
        ]
        rule_parsers_expr = "[" + ", ".join(rule_parsers) + "]"
        rule_kwargs = {
            key: value for key, value in skw.items() if key != "index_parsers"
        }
        lines.append(f"RS = RuleSites(index_parsers={rule_parsers_expr}, **{_py(rule_kwargs)})")
        sites_var = "RS"
    else:
        raise ValueError(f"Unsupported sites type: {stype}")
    lines.append("")

    name_to_var: Dict[str, str] = {}

    for step in steps:
        op = step["op"]
        name = step["name"]
        comment = step.get("comment")
        if comment:
            lines.append(f"# {comment}")

        if op == "ModifierAdder":
            seq = step["sequence"]
            kwargs = step.get("kwargs", {})
            seq_vars = [name_to_var[s] for s in seq]
            lines.append(
                f"{name} = ModifierAdder([{', '.join(seq_vars)}], **{_py(kwargs)})"
            )
            name_to_var[name] = name
            lines.append("")
            continue

        kwargs = step.get("kwargs", {})
        ctor = op
        lines.append(f"{name} = {ctor}({sites_var}, **{_py(kwargs)})")
        name_to_var[name] = name
        lines.append("")

    lines += [
        "# Example execution:",
        "# modified = " + steps[-1]["name"] + ".get_modified_atoms(atoms)",
        "",
    ]

    return "\n".join(lines)


# 5) dry_run (subprocess exec, capture tracebacks)
def dry_run(python_code: str, timeout_s: int = 30) -> Dict[str, Any]:
    """
    Runs python_code in a fresh interpreter process.
    Expects the environment to have gg + ase installed, and your code to define atoms if needed.
    """
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "run.py"
        path.write_text(python_code, encoding="utf-8")

        try:
            proc = subprocess.run(
                ["python", str(path)],
                capture_output=True,
                text=True,
                timeout=timeout_s,
                check=False,
            )
            ok = proc.returncode == 0
            return {
                "ok": ok,
                "returncode": proc.returncode,
                "stdout": proc.stdout[-4000:],
                "stderr": proc.stderr[-8000:],
            }
        except subprocess.TimeoutExpired as exc:
            return {
                "ok": False,
                "returncode": None,
                "stdout": "",
                "stderr": f"Timeout after {timeout_s}s: {exc}",
            }

# 6) Simple orchestrator loop (retrieve -> emit -> compile -> dry_run -> revise)
def run_nl_to_code(
    nl_request: str,
    atoms_code: str,                 # <-- user-provided python that defines `atoms`
    model: str = "gpt-4o-mini",
    max_iters: int = 3,
    timeout_s: int = 30,
) -> Dict[str, Any]:
    """
    End-to-end NL -> Plan -> gg code. `atoms_code` must define a variable named `atoms`
    (ASE Atoms) that gg modifiers can operate on.
    """
    client = OpenAI()

    snippets = retrieve_gg_docs(query=nl_request, top_k=6)
    plan = emit_plan(client=client, model=model, nl_request=nl_request, snippets=snippets)

    code_with_atoms = ""
    result: Dict[str, Any] = {}

    for i in range(max_iters):
        validate_plan(plan)

        code = compile_plan(plan)

        code_with_atoms = (
            textwrap.dedent(atoms_code).rstrip()
            + "\n\n"
            + code
            + "\n\n"
            + f"{plan['steps'][-1]['name']}.get_modified_atoms(atoms)\n"
        )

        result = dry_run(code_with_atoms, timeout_s=timeout_s)

        if result["ok"]:
            return {"ok": True, "plan": plan, "python_code": code_with_atoms, "dry_run": result}

        errors = (
            f"dry_run failed (iter {i+1}/{max_iters}). stderr:\n{result['stderr']}\n"
            f"stdout:\n{result['stdout']}"
        )
        plan = revise(client=client, model=model, plan=plan, errors=errors, snippets=snippets)

    return {"ok": False, "plan": plan, "python_code": code_with_atoms, "dry_run": result}
