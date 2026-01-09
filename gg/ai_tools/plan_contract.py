"""
Contract for translating natural language instructions into gg plan steps.

This file defines a **JSON-serializable** plan schema (Pydantic v2) that an LLM
can output, and your orchestrator can validate, then compile into gg objects.

It is aligned to the public gg docs (graph-gcbh readthedocs) for:
- Sites: FlexibleSites, SurfaceSites
- Modifiers: Add, AddBi, Remove, Replace, Swap, ClusterRotate, ClusterTranslate, ModifierAdder
"""

from __future__ import annotations
from typing import Annotated, Literal, Union
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

#Helpers for inputs
IdentifierStr = Annotated[str, Field(pattern=r"^[A-Za-z_][A-Za-z0-9_]*$")]
SurfCoord = Union[int, list[int], list[list[int]]]
SymbolList = Union[str, list[str]]

# Base / shared config
class BaseKwargs(BaseModel):
    model_config = ConfigDict(extra="forbid")

# Sites specs (gg.predefined_sites)
class FlexibleSitesKwargs(BaseKwargs):
    constraints: bool | None = False
    index: list[int] | None = None
    tag: int | None = None
    opp_tag: bool | None = False
    com: float | None = None
    max_bond_ratio: float | None = 1.2
    max_bond: float | None = 0.0
    contact_error: float | None = 0.3


class SurfaceSitesKwargs(BaseKwargs):
    max_coord: dict[str, int]
    max_bond_ratio: float | None = 1.2
    max_bond: float | None = 0.0
    contact_error: float | None = 0.3
    com: float | None = 0.1

class FlexibleSitesSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["FlexibleSites"]
    kwargs: FlexibleSitesKwargs = Field(default_factory=FlexibleSitesKwargs)

class SurfaceSitesSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: Literal["SurfaceSites"]
    kwargs: SurfaceSitesKwargs

SitesSpec = Annotated[
    Union[FlexibleSitesSpec, SurfaceSitesSpec],
    Field(discriminator="type"),
]

# Modifier kwargs (gg.modifiers.*)
class CommonModifierKwargs(BaseKwargs):
    print_movie: bool | None = False
    unique: bool | None = True
    unique_method: str | None = "fullgraph"
    unique_depth: int | None = 3
    weight: float | None = 1.0

class AddKwargs(CommonModifierKwargs):
    ads: str
    surf_coord: SurfCoord
    surf_sym: SymbolList
    ads_id: SymbolList | None = None
    ads_dist: float | None = None
    ads_rotate: bool | None = True
    normal_method: str | None = "svd"
    tag: bool | None = True

    @model_validator(mode="after")
    def _validate_add(self) -> "AddKwargs":
        # Ensure lists are non-empty if provided as lists
        if isinstance(self.surf_sym, list) and len(self.surf_sym) == 0:
            raise ValueError(
                "Add.surf_sym must be a non-empty symbol list (or a single symbol string).")
        if isinstance(self.ads_id, list) and len(self.ads_id) == 0:
            raise ValueError(
                "Add.ads_id must be a non-empty symbol list (or a single symbol string).")
        return self

class AddBiKwargs(CommonModifierKwargs):
    ads: str
    surf_coord: SurfCoord
    surf_sym: SymbolList
    ads_id: SymbolList
    ads_dist: float | str | None = None
    ads_rotate: bool | None = True
    add_ads_error: float | None = 0.5
    normal_method: str | None = "mean"
    tag: bool | None = True

    @model_validator(mode="after")
    def _validate_addbi(self) -> "AddBiKwargs":
        if isinstance(self.surf_sym, list) and len(self.surf_sym) == 0:
            raise ValueError(
                "AddBi.surf_sym must be a non-empty symbol list (or a single symbol string).")
        if isinstance(self.ads_id, list) and len(self.ads_id) == 0:
            raise ValueError(
                "AddBi.ads_id must be a non-empty symbol list (or a single symbol string).")
        # If bidentate is expressed as a nested list, it should look like [[...],[...]]
        if isinstance(self.surf_coord, list) and self.surf_coord and isinstance(self.surf_coord[0], list):
            if len(self.surf_coord) != 2:
                raise ValueError(
                    "AddBi.surf_coord as list[list[int]] must have exactly two site entries.")
        return self

class RemoveKwargs(BaseKwargs):
    # gg docs: Remove(surface_sites, to_del: Atoms | str, ...)
    # Keep JSON-serializable: require str (e.g. "OH", "HCOO", "O")
    to_del: str
    max_bond_ratio: float | None = 1.2
    max_bond: float | None = 0.0
    allow_external_symbols: list[str] | None = None
    max_external_neighbors: int | None = None
    print_movie: bool | None = False
    unique: bool | None = False
    unique_method: str | None = "fullgraph"
    unique_depth: int | None = 3
    weight: float | None = 1.0

class ReplaceKwargs(BaseKwargs):
    # gg docs: Replace(surface_sites, to_del: Atoms|str, with_replace: Atoms|str, ...)
    to_del: str
    with_replace: str
    max_bond_ratio: float | None = 1.2
    max_bond: float | None = 0.0
    print_movie: bool | None = False
    unique: bool | None = True
    unique_method: str | None = "fullgraph"
    unique_depth: int | None = 3
    weight: float | None = 1.0

class SwapKwargs(CommonModifierKwargs):
    # gg docs: Swap(surface_sites, swap_sym: list, swap_ind: list=None, ...)
    swap_sym: list[str]
    swap_ind: list[int] | None = None

    @model_validator(mode="after")
    def _validate_swap(self) -> "SwapKwargs":
        if len(self.swap_sym) != 2:
            raise ValueError("Swap.swap_sym must contain exactly two element symbols (e.g. ['Ni','N']).")
        if self.swap_ind is not None and len(self.swap_ind) != 2:
            raise ValueError("Swap.swap_ind must contain exactly two indices if provided.")
        return self

class ClusterRotateKwargs(CommonModifierKwargs):
    max_angle: float | int | None = 180
    rotate_vector: tuple[float, float, float] | None = None
    contact_error: float | None = 0.2
    nmovie: int | None = 1

class ClusterTranslateKwargs(CommonModifierKwargs):
    max_displace: float | None = 5.0
    allowed_direction: tuple[bool, bool, bool] | list[bool] | None = (True, True, False)
    contact_error: float | None = 0.2
    nmovie: int | None = 1

    @model_validator(mode="after")
    def _validate_allowed_direction(self) -> "ClusterTranslateKwargs":
        if isinstance(self.allowed_direction, list):
            if len(self.allowed_direction) != 3:
                raise ValueError(
                    "ClusterTranslate.allowed_direction must have length 3 if provided as a list.")
            self.allowed_direction = (bool(self.allowed_direction[0]),
                                      bool(self.allowed_direction[1]),
                                      bool(self.allowed_direction[2]))
        return self


class ModifierAdderKwargs(CommonModifierKwargs):
    # gg docs: ModifierAdder(modifier_instances: list[ParentModifier], ...)
    # In-plan we reference earlier steps by name; compiler will materialize instances.
    max_bond_ratio: float | None = 1.2
    max_bond: float | None = 0.0

# Steps
class StepBase(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: IdentifierStr
    comment: str | None = None

class AddStep(StepBase):
    op: Literal["Add"]
    kwargs: AddKwargs

class AddBiStep(StepBase):
    op: Literal["AddBi"]
    kwargs: AddBiKwargs

class RemoveStep(StepBase):
    op: Literal["Remove"]
    kwargs: RemoveKwargs

class ReplaceStep(StepBase):
    op: Literal["Replace"]
    kwargs: ReplaceKwargs

class SwapStep(StepBase):
    op: Literal["Swap"]
    kwargs: SwapKwargs

class ClusterRotateStep(StepBase):
    op: Literal["ClusterRotate"]
    kwargs: ClusterRotateKwargs = Field(default_factory=ClusterRotateKwargs)

class ClusterTranslateStep(StepBase):
    op: Literal["ClusterTranslate"]
    kwargs: ClusterTranslateKwargs = Field(default_factory=ClusterTranslateKwargs)

class ModifierAdderStep(StepBase):
    op: Literal["ModifierAdder"]
    modifier_instances: list[IdentifierStr]
    kwargs: ModifierAdderKwargs = Field(default_factory=ModifierAdderKwargs)

    @model_validator(mode="after")
    def validate_modifier_instances(self) -> "ModifierAdderStep":
        if not self.modifier_instances:
            raise ValueError(
                "ModifierAdder.modifier_instances must be a non-empty list of previous step names.")
        return self

Step = Annotated[
    Union[
        AddStep,
        AddBiStep,
        RemoveStep,
        ReplaceStep,
        SwapStep,
        ClusterRotateStep,
        ClusterTranslateStep,
        ModifierAdderStep,
    ],
    Field(discriminator="op"),
]

# Plan
class Plan(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sites: SitesSpec
    steps: list[Step]

    @model_validator(mode="after")
    def validate_steps(self) -> "Plan":
        if not self.steps:
            raise ValueError("Plan requires at least one step.")
        # Ensure unique step names (so a compiler can reference them)
        names = [s.name for s in self.steps]
        if len(names) != len(set(names)):
            raise ValueError("Step names must be unique.")
        return self

def validate_plan(plan_dict: dict) -> Plan:
    """Validate a plan dictionary and return the parsed Plan."""
    try:
        return Plan.model_validate(plan_dict)
    except ValidationError as exc:
        raise ValueError(f"Invalid plan: {exc}") from exc

def _enforce_required_properties(schema: object) -> object:
    if isinstance(schema, dict):
        if schema.get("type") == "object" and "properties" in schema:
            properties = schema.get("properties", {})
            required = schema.get("required", [])
            if isinstance(properties, dict):
                required_keys = set(properties.keys())
                existing = set(required) if isinstance(required, list) else set()
                schema["required"] = sorted(required_keys | existing)
        for key in ("allOf", "anyOf", "oneOf"):
            items = schema.get(key)
            if isinstance(items, list):
                for item in items:
                    _enforce_required_properties(item)
        items = schema.get("items")
        if items is not None:
            _enforce_required_properties(items)
        prefix_items = schema.get("prefixItems")
        if isinstance(prefix_items, list):
            for item in prefix_items:
                _enforce_required_properties(item)
        for defs_key in ("$defs", "definitions"):
            defs = schema.get(defs_key)
            if isinstance(defs, dict):
                for item in defs.values():
                    _enforce_required_properties(item)
    elif isinstance(schema, list):
        for item in schema:
            _enforce_required_properties(item)
    return schema


def plan_json_schema() -> dict:
    """Return the JSON Schema for the Plan contract."""
    schema = Plan.model_json_schema(mode="validation")
    return _enforce_required_properties(schema)

# Minimal tests (run as script)
def _tests() -> None:
    valid_add = {
        "sites": {"type": "FlexibleSites", "kwargs": {"constraints": True, "max_bond_ratio": 1.2}},
        "steps": [
            {
                "op": "Add",
                "name": "add_oh",
                "kwargs": {
                    "ads": "OH",
                    "surf_coord": [1, 2, 3],
                    "ads_id": ["O"],
                    "surf_sym": ["Pt"],
                    "unique": True,
                    "unique_method": "fullgraph",
                },
            }
        ],
    }
    assert validate_plan(valid_add).steps[0].op == "Add"

    valid_add_bi = {
        "sites": {"type": "FlexibleSites", "kwargs": {"constraints": True}},
        "steps": [
            {
                "op": "AddBi",
                "name": "add_formate",
                "kwargs": {
                    "ads": "HCOO",
                    "surf_coord": [[1], [2]],
                    "ads_id": ["O", "O"],
                    "surf_sym": ["Pt"],
                },
            }
        ],
    }
    assert validate_plan(valid_add_bi).steps[0].op == "AddBi"

    valid_remove = {
        "sites": {"type": "SurfaceSites", "kwargs": {"max_coord": {"Pt": 12, "O": 4, "H": 2}}},
        "steps": [
            {
                "op": "Remove",
                "name": "remove_oh",
                "kwargs": {"to_del": "OH"},
            }
        ],
    }
    assert validate_plan(valid_remove).steps[0].op == "Remove"

    valid_swap = {
        "sites": {"type": "SurfaceSites", "kwargs": {"max_coord": {"Pt": 12, "Ni": 12, "N": 3}}},
        "steps": [
            {"op": "Swap", "name": "swap_ni_n", "kwargs": {"swap_sym": ["Ni", "N"]}},
        ],
    }
    assert validate_plan(valid_swap).steps[0].op == "Swap"

    valid_modifier_adder = {
        "sites": {"type": "FlexibleSites", "kwargs": {"constraints": True}},
        "steps": [
            {
                "op": "Add",
                "name": "add_o",
                "kwargs": {"ads": "O", "surf_coord": [1], "surf_sym": ["Pt"], "ads_id": ["O"]},
            },
            {
                "op": "Add",
                "name": "add_h",
                "kwargs": {"ads": "H", "surf_coord": [1], "surf_sym": ["Pt"], "ads_id": ["H"]},
            },
            {
                "op": "ModifierAdder",
                "name": "add_oh",
                "modifier_instances": ["add_o", "add_h"],
                "kwargs": {"unique": True, "print_movie": False},
            },
        ],
    }
    assert validate_plan(valid_modifier_adder).steps[-1].op == "ModifierAdder"

    invalid_unknown_op = {
        "sites": {"type": "FlexibleSites", "kwargs": {}},
        "steps": [{"op": "Unknown", "name": "bad", "kwargs": {}}],
    }
    try:
        validate_plan(invalid_unknown_op)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected unknown op to fail")

    invalid_swap_len = {
        "sites": {"type": "FlexibleSites", "kwargs": {}},
        "steps": [{"op": "Swap", "name": "bad_swap", "kwargs": {"swap_sym": ["Ni"]}}],
    }
    try:
        validate_plan(invalid_swap_len)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected bad swap_sym length to fail")

    invalid_allowed_direction = {
        "sites": {"type": "FlexibleSites", "kwargs": {}},
        "steps": [{"op": "ClusterTranslate",
                   "name": "bad_ct",
                   "kwargs": {"allowed_direction": [True, False]}}],
    }
    try:
        validate_plan(invalid_allowed_direction)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected bad allowed_direction length to fail")

if __name__ == "__main__":
    _tests()
