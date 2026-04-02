"""
Disease profile loader for the tolerogenic epitope pipeline.

A disease profile is a JSON file in ``data/diseases/{disease_id}.json``
that specifies everything the pipeline needs to run for a particular
autoimmune disease: target antigens, IEDB disease filter, gold standard
path, and calibration specification.

Adding a new disease requires only creating a new JSON profile and the
corresponding gold standard file — no code changes.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

_PROFILES_DIR = Path("data/diseases")

_REQUIRED_FIELDS = [
    "disease_id",
    "disease_name",
    "iedb_disease_filter",
    "primary_antigen",
    "antigens",
    "gold_standard_path",
    "calibration",
]


def load_disease_profile(
    disease_id: str,
    profiles_dir: str | os.PathLike[str] = _PROFILES_DIR,
) -> dict[str, Any]:
    """Load and validate a disease profile JSON.

    Parameters
    ----------
    disease_id :
        Short identifier (e.g. ``"itp"``, ``"ms"``).
    profiles_dir :
        Directory containing ``{disease_id}.json`` files.

    Returns
    -------
    The parsed profile dict.

    Raises
    ------
    FileNotFoundError
        If the profile JSON does not exist.
    ValueError
        If required fields are missing.
    """
    path = Path(profiles_dir) / f"{disease_id}.json"
    if not path.exists():
        available = [p.stem for p in Path(profiles_dir).glob("*.json")]
        raise FileNotFoundError(
            f"Disease profile '{disease_id}' not found at {path}. "
            f"Available profiles: {available}"
        )

    profile = json.loads(path.read_text())

    missing = [f for f in _REQUIRED_FIELDS if f not in profile]
    if missing:
        raise ValueError(
            f"Disease profile '{disease_id}' is missing required fields: {missing}"
        )

    return profile


def get_antigens(profile: dict[str, Any]) -> list[dict[str, str]]:
    """Return the list of antigen dicts from a disease profile."""
    return profile["antigens"]


def get_antigen_dict(profile: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Build a ``{uniprot_id: metadata}`` dict matching ``ITP_ANTIGENS`` format."""
    return {
        ag["uniprot_id"]: {k: v for k, v in ag.items() if k != "uniprot_id"}
        for ag in profile["antigens"]
    }


def get_primary_antigen(profile: dict[str, Any]) -> str:
    """Return the primary antigen UniProt ID."""
    return profile["primary_antigen"]


def get_iedb_filter(profile: dict[str, Any]) -> str:
    """Return the IEDB disease filter string."""
    return profile["iedb_disease_filter"]


def get_gold_standard_path(profile: dict[str, Any]) -> str:
    """Return the path to the gold standard JSON."""
    return profile["gold_standard_path"]


def get_calibration_spec(profile: dict[str, Any]) -> dict[str, Any]:
    """Return the calibration specification dict."""
    return profile["calibration"]
