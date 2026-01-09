from __future__ import annotations

import json
import logging
import os
import re
import subprocess
import tempfile
import textwrap
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List
from urllib.request import urlopen
from openai import OpenAI
from .plan_contract import validate_plan

LOGGER = logging.getLogger(__name__)


def _log_llm_call(phase: str, model: str, payload: Dict[str, Any], output_text: str) -> None:
    LOGGER.info(
        "LLM %s call model=%s input=%s output=%s",
        phase,
        model,
        json.dumps(payload, ensure_ascii=False),
        output_text,
    )

# 1) Tool Specs (what you pass to the model)
TOOLS: List[Dict[str, Any]] = [
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

# 2) Retrieval (modifiers RAG)
MODIFIERS_RAG_URL = "https://graph-gcbh.readthedocs.io/en/latest/modifiers.html"


def _strip_html(value: str) -> str:
    no_script = re.sub(r"(?is)<(script|style).*?>.*?</\1>", " ", value)
    no_tags = re.sub(r"(?is)<[^>]+>", " ", no_script)
    cleaned = re.sub(r"\s+", " ", no_tags)
    return cleaned.strip()


@lru_cache(maxsize=1)
def _load_modifiers_rag() -> str:
    with urlopen(MODIFIERS_RAG_URL) as response:
        html = response.read().decode("utf-8", errors="ignore")
    return _strip_html(html)


def get_modifiers_rag_snippets() -> List[Dict[str, str]]:
    return [{"source": MODIFIERS_RAG_URL, "text": _load_modifiers_rag()}]

# 3) emit_plan + revise
SYSTEM_PLANNER = """You translate natural language catalyst-surface modification requests into a STRICT PlanJSON.
Rules:
- Output MUST be a JSON object with `sites` and `steps`.
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
    log_payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": SYSTEM_PLANNER},
            {"role": "user", "content": json.dumps(user_payload)},
        ],
    }

    resp = client.responses.create(
        model=model,
        input=[
            {"role": "system", "content": SYSTEM_PLANNER},
            {"role": "user", "content": json.dumps(user_payload)},
        ],
    )
    out = resp.output_text
    _log_llm_call("emit_plan", model, log_payload, out)
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
    log_payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": SYSTEM_PLANNER},
            {"role": "user", "content": json.dumps(user_payload)},
        ],
    }
    resp = client.responses.create(
        model=model,
        input=[
            {"role": "system", "content": SYSTEM_PLANNER},
            {"role": "user", "content": json.dumps(user_payload)},
        ],
    )
    out = resp.output_text
    _log_llm_call("revise", model, log_payload, out)
    return json.loads(out)

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
            seq = step["modifier_instances"]
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

# 6) Simple orchestrator loop (RAG -> emit -> compile -> dry_run -> revise)
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
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO)
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise ValueError("OPENAI_API_KEY must be set to use the OpenAI client.")
    client = OpenAI(api_key=api_key)

    snippets = get_modifiers_rag_snippets()
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
