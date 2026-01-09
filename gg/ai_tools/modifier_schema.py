from __future__ import annotations

MODIFIER_SCHEMA = {
    "title": "gg.modifiers (lightweight)",
    "type": "object",
    "description": "Minimal modifier kwargs distilled from docs/source/modifiers.rst.",
    "properties": {
        "Add": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "ads": {"type": "string"},
                "surf_coord": {"type": "array", "items": {"type": "integer"}},
                "ads_id": {"type": "array", "items": {"type": "string"}},
                "surf_sym": {"type": "array", "items": {"type": "string"}},
                "print_movie": {"type": "boolean"},
                "unique": {"type": "boolean"},
            },
        },
        "AddBi": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "ads": {"type": "string"},
                "surf_coord": {"type": "array", "items": {"type": "integer"}},
                "ads_id": {"type": "array", "items": {"type": "string"}},
                "surf_sym": {"type": "array", "items": {"type": "string"}},
                "print_movie": {"type": "boolean"},
                "unique": {"type": "boolean"},
            },
        },
        "Remove": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "to_del": {"type": "string"},
                "print_movie": {"type": "boolean"},
            },
        },
        "Replace": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "to_del": {"type": "string"},
                "with_replace": {"type": "string"},
                "print_movie": {"type": "boolean"},
            },
        },
        "Swap": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "swap_sym": {"type": "array", "items": {"type": "string"}},
                "print_movie": {"type": "boolean"},
            },
        },
        "ClusterRotate": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "contact_error": {"type": "number"},
                "rotate_vector": {"type": "array", "items": {"type": "number"}},
            },
        },
        "ClusterTranslate": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "contact_error": {"type": "number"},
            },
        },
        "ModifierAdder": {
            "type": "object",
            "properties": {
                "sites": {"type": "string", "description": "Site class instance name."},
                "modifier_instances": {"type": "array", "items": {"type": "string"}},
                "print_movie": {"type": "boolean"},
                "unique": {"type": "boolean"},
            },
        },
    },
}
