"""
ITP target antigen definitions.

Maps UniProt accession numbers to metadata for the six platelet surface
glycoproteins targeted by autoantibodies in Immune Thrombocytopenic
Purpura (ITP).  Two protein complexes are represented:

  GPIIb/IIIa  — integrin αIIbβ3, the primary autoantibody target
  GPIb-IX-V   — the von Willebrand factor receptor complex
"""

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
