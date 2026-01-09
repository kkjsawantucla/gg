"""
Lightweight validation for gg plans.

This module keeps validation intentionally simple so the planner can emit a
compact JSON object without a full JSON schema contract.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List

IDENTIFIER_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

ALLOWED_SITE_TYPES = {"FlexibleSites", "SurfaceSites", "RuleSites"}
ALLOWED_OPS = {
    "Add",
    "AddBi",
    "Remove",
    "Replace",
    "Swap",
    "Rattle",
    "Translate",
    "ClusterRotate",
    "ClusterTranslate",
    "ModifierAdder",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _validate_identifier(name: Any, label: str) -> str:
    _require(isinstance(name, str), f"{label} must be a string.")
    _require(bool(IDENTIFIER_RE.match(name)), f"{label} must be a valid identifier.")
    return name


def _validate_sites(sites: Any) -> Dict[str, Any]:
    _require(isinstance(sites, dict), "sites must be an object.")
    site_type = sites.get("type")
    _require(site_type in ALLOWED_SITE_TYPES, "sites.type must be a supported sites type.")
    kwargs = sites.get("kwargs", {})
    _require(isinstance(kwargs, dict), "sites.kwargs must be an object.")
    return {"type": site_type, "kwargs": kwargs}


def _validate_steps(steps: Any) -> List[Dict[str, Any]]:
    _require(isinstance(steps, list), "steps must be a list.")
    _require(bool(steps), "steps must contain at least one step.")

    seen_names: set[str] = set()
    validated_steps: List[Dict[str, Any]] = []

    for idx, step in enumerate(steps, start=1):
        _require(isinstance(step, dict), f"steps[{idx}] must be an object.")
        op = step.get("op")
        _require(op in ALLOWED_OPS, f"steps[{idx}].op must be a supported operation.")

        name = _validate_identifier(step.get("name"), f"steps[{idx}].name")
        _require(name not in seen_names, f"steps[{idx}].name must be unique.")
        seen_names.add(name)

        comment = step.get("comment")
        if comment is not None:
            _require(isinstance(comment, str), f"steps[{idx}].comment must be a string.")

        if op == "ModifierAdder":
            modifier_instances = step.get("modifier_instances")
            _require(
                isinstance(modifier_instances, list) and modifier_instances,
                f"steps[{idx}].modifier_instances must be a non-empty list.",
            )
            for ref in modifier_instances:
                _validate_identifier(ref, f"steps[{idx}].modifier_instances entry")
                _require(
                    ref in seen_names,
                    f"steps[{idx}].modifier_instances must reference earlier step names.",
                )
        else:
            modifier_instances = None

        kwargs = step.get("kwargs", {})
        _require(isinstance(kwargs, dict), f"steps[{idx}].kwargs must be an object.")

        normalized_step: Dict[str, Any] = {
            "op": op,
            "name": name,
            "kwargs": kwargs,
        }
        if comment is not None:
            normalized_step["comment"] = comment
        if modifier_instances is not None:
            normalized_step["modifier_instances"] = modifier_instances
        validated_steps.append(normalized_step)

    return validated_steps


def validate_plan(plan_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Validate a plan dictionary and return a normalized plan dict."""
    _require(isinstance(plan_dict, dict), "plan must be an object.")
    sites = _validate_sites(plan_dict.get("sites"))
    steps = _validate_steps(plan_dict.get("steps"))
    return {"sites": sites, "steps": steps}
