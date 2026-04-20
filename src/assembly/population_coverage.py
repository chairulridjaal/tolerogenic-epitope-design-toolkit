"""
HLA population coverage computation.

Computes what fraction of a given population can present at least one
(or at least k) epitopes from a multi-epitope construct, based on
published HLA allele frequencies and the peptide-allele binding data
from Phase 2 MHC-II predictions.

Formula: Bui et al. 2006, BMC Bioinformatics 7:153.
Allele frequencies: Allele Frequency Net Database (AFND) / published sources.
"""

from __future__ import annotations

from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# HLA allele frequencies by population (gene/allele frequency, not phenotype)
# Sources: AFND (allelefrequencies.net), Gonzalez-Galarza et al. 2020
#
# For DQ and DP heterodimers, we use the beta-chain frequency as an
# approximation (alpha-beta linkage disequilibrium is strong).
# ---------------------------------------------------------------------------

HLA_FREQUENCIES: dict[str, dict[str, float]] = {
    # HLA-DRB1 alleles
    "HLA-DRB1*01:01": {
        "European": 0.095, "East_Asian": 0.020, "African": 0.030, "South_Asian": 0.045,
    },
    "HLA-DRB1*03:01": {
        "European": 0.110, "East_Asian": 0.020, "African": 0.080, "South_Asian": 0.065,
    },
    "HLA-DRB1*04:01": {
        "European": 0.090, "East_Asian": 0.030, "African": 0.015, "South_Asian": 0.040,
    },
    "HLA-DRB1*04:05": {
        "European": 0.010, "East_Asian": 0.055, "African": 0.005, "South_Asian": 0.035,
    },
    "HLA-DRB1*07:01": {
        "European": 0.125, "East_Asian": 0.050, "African": 0.060, "South_Asian": 0.065,
    },
    "HLA-DRB1*09:01": {
        "European": 0.015, "East_Asian": 0.130, "African": 0.020, "South_Asian": 0.030,
    },
    "HLA-DRB1*11:01": {
        "European": 0.065, "East_Asian": 0.030, "African": 0.070, "South_Asian": 0.080,
    },
    "HLA-DRB1*13:01": {
        "European": 0.055, "East_Asian": 0.025, "African": 0.045, "South_Asian": 0.050,
    },
    "HLA-DRB1*15:01": {
        "European": 0.130, "East_Asian": 0.080, "African": 0.050, "South_Asian": 0.100,
    },
    # DQ alleles (beta-chain frequency as proxy)
    "HLA-DQA1*01:01/DQB1*05:01": {
        "European": 0.100, "East_Asian": 0.050, "African": 0.080, "South_Asian": 0.070,
    },
    "HLA-DQA1*01:02/DQB1*06:02": {
        "European": 0.100, "East_Asian": 0.060, "African": 0.040, "South_Asian": 0.080,
    },
    # DP allele
    "HLA-DPA1*01:03/DPB1*04:01": {
        "European": 0.200, "East_Asian": 0.100, "African": 0.100, "South_Asian": 0.150,
    },
}

POPULATIONS = ["European", "East_Asian", "African", "South_Asian"]


def allele_to_phenotype_freq(allele_freq: float) -> float:
    """Convert allele frequency to phenotype frequency via Hardy-Weinberg.

    An individual is positive for an allele if they carry at least one copy.
    P(phenotype) = 1 - (1 - p_allele)^2
    """
    return 1.0 - (1.0 - allele_freq) ** 2


def compute_coverage(
    epitope_alleles: dict[str, list[str]],
    population: str,
) -> float:
    """Compute fraction of population that can present at least one epitope.

    Parameters
    ----------
    epitope_alleles :
        Dict mapping peptide → list of alleles it binds (rank ≤ 10%).
    population :
        Population name (key in HLA_FREQUENCIES).

    Returns
    -------
    Coverage fraction in [0, 1].

    Reference: Bui et al. 2006, BMC Bioinformatics 7:153.
    Formula: Coverage = 1 - prod_over_alleles(1 - p_phenotype_i)
    where the product runs over all UNIQUE alleles that bind at least
    one epitope in the set.
    """
    # Collect unique alleles that bind any epitope
    covered_alleles: set[str] = set()
    for alleles in epitope_alleles.values():
        covered_alleles.update(alleles)

    # Compute coverage
    prob_not_covered = 1.0
    for allele in covered_alleles:
        freq_data = HLA_FREQUENCIES.get(allele, {})
        allele_freq = freq_data.get(population, 0.0)
        pheno_freq = allele_to_phenotype_freq(allele_freq)
        prob_not_covered *= (1.0 - pheno_freq)

    return 1.0 - prob_not_covered


def compute_coverage_table(
    epitope_alleles: dict[str, list[str]],
) -> dict[str, float]:
    """Compute coverage across all populations.

    Returns dict mapping population → coverage fraction.
    """
    return {
        pop: compute_coverage(epitope_alleles, pop)
        for pop in POPULATIONS
    }


def get_epitope_alleles(
    peptides: list[str],
    predictions_df: "pd.DataFrame",
    threshold: float = 10.0,
) -> dict[str, list[str]]:
    """Extract which alleles each peptide binds from predictions.

    Parameters
    ----------
    peptides :
        List of peptide sequences.
    predictions_df :
        DataFrame with columns peptide, allele, percentile_rank.
    threshold :
        Maximum percentile rank to count as a binder.

    Returns
    -------
    Dict mapping peptide → list of allele names.
    """
    import pandas as pd

    result: dict[str, list[str]] = {}
    for pep in peptides:
        rows = predictions_df[
            (predictions_df["peptide"] == pep)
            & (predictions_df["percentile_rank"] <= threshold)
        ]
        result[pep] = rows["allele"].unique().tolist()

    return result


def format_coverage_report(
    epitope_alleles: dict[str, list[str]],
    construct_name: str = "",
) -> str:
    """Generate a formatted coverage report string."""
    coverage = compute_coverage_table(epitope_alleles)

    lines = []
    if construct_name:
        lines.append(f"Population Coverage — {construct_name}")
    else:
        lines.append("Population Coverage")
    lines.append("-" * 45)

    n_epitopes = len(epitope_alleles)
    n_alleles = len(set(a for alleles in epitope_alleles.values() for a in alleles))
    lines.append(f"Epitopes: {n_epitopes}, covering {n_alleles} unique alleles")
    lines.append("")

    for pop in POPULATIONS:
        cov = coverage[pop]
        lines.append(f"  {pop:<15}: {cov:.1%}")

    return "\n".join(lines)
