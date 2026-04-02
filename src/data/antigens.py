"""
Antigen definitions for the tolerogenic epitope pipeline.

``ITP_ANTIGENS`` is the original hardcoded dict for backward compatibility.
New code should use ``load_antigens_from_profile()`` with a disease profile
from ``src.data.disease_profile``.
"""

from __future__ import annotations

from typing import Any

# Legacy ITP antigen dict — kept for backward compatibility with existing
# notebooks and imports.  Do not delete.
ITP_ANTIGENS: dict[str, dict[str, str]] = {
    "P08514": {
        "gene": "ITGA2B",
        "name": "Integrin subunit alpha-IIb",
        "complex": "GPIIb/IIIa",
    },
    "P05106": {
        "gene": "ITGB3",
        "name": "Integrin subunit beta-3",
        "complex": "GPIIb/IIIa",
    },
    "P07359": {
        "gene": "GP1BA",
        "name": "Platelet glycoprotein Ib alpha chain",
        "complex": "GPIb-IX-V",
    },
    "P13224": {
        "gene": "GP1BB",
        "name": "Platelet glycoprotein Ib beta chain",
        "complex": "GPIb-IX-V",
    },
    "P14770": {
        "gene": "GP9",
        "name": "Platelet glycoprotein IX",
        "complex": "GPIb-IX-V",
    },
    "P40197": {
        "gene": "GP5",
        "name": "Platelet glycoprotein V",
        "complex": "GPIb-IX-V",
    },
}


def load_antigens_from_profile(profile: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Build an antigen dict from a disease profile.

    Returns a dict matching the ``ITP_ANTIGENS`` format:
    ``{uniprot_id: {"gene": ..., "name": ..., ...}}``.
    """
    from src.data.disease_profile import get_antigen_dict
    return get_antigen_dict(profile)
