"""
Peptide scanning utilities for the tolerogenic epitope pipeline.

Generates overlapping peptide windows from a protein sequence and
provides benchmark-recovery scoring against known epitopes.
Standard library only — no external dependencies.
"""

from __future__ import annotations


def scan_sequence(
    sequence: str,
    min_len: int = 9,
    max_len: int = 25,
) -> list[dict[str, str | int]]:
    """Slide windows of every length from *min_len* to *max_len* over *sequence*.

    Parameters
    ----------
    sequence:
        Full amino-acid sequence string.
    min_len:
        Shortest peptide to generate (default 9).
    max_len:
        Longest peptide to generate (default 25).

    Returns
    -------
    List of dicts, each with keys ``peptide``, ``start`` (0-based),
    ``end`` (exclusive), and ``length``.
    """
    peptides: list[dict[str, str | int]] = []
    seq_len = len(sequence)
    for length in range(min_len, max_len + 1):
        for start in range(seq_len - length + 1):
            end = start + length
            peptides.append({
                "peptide": sequence[start:end],
                "start": start,
                "end": end,
                "length": length,
            })
    return peptides


def filter_peptides(
    peptides: list[dict[str, str | int]],
    min_len: int = 15,
    max_len: int = 15,
) -> list[dict[str, str | int]]:
    """Keep only peptides whose length is within [*min_len*, *max_len*].

    Parameters
    ----------
    peptides:
        Output of :func:`scan_sequence`.
    min_len:
        Minimum peptide length to keep (inclusive, default 15).
    max_len:
        Maximum peptide length to keep (inclusive, default 15).

    Returns
    -------
    Filtered list of peptide dicts.
    """
    return [p for p in peptides if min_len <= p["length"] <= max_len]


def benchmark_recovery(
    predicted_peptides: list[str],
    known_epitopes: list[str],
    top_n_percent: float = 10.0,
) -> dict[str, float | int]:
    """Measure what fraction of *known_epitopes* are recovered in the top predictions.

    A known epitope is considered *recovered* if any predicted peptide in
    the top-N% set either contains the epitope as a substring or is itself
    a substring of the epitope.  This accommodates length mismatches
    between fixed-window predictions and variable-length experimental
    epitopes.

    Parameters
    ----------
    predicted_peptides:
        Peptide sequences sorted by binding score (best/strongest first).
    known_epitopes:
        Experimentally characterised epitope sequences.
    top_n_percent:
        Percentage of the prediction list to consider (default 10).

    Returns
    -------
    Dict with keys:

    - ``recovery_rate`` — fraction of known epitopes found (0.0–1.0)
    - ``recovered`` — count of known epitopes found
    - ``total_known`` — total known epitopes evaluated
    - ``top_n_count`` — how many predictions were in the top-N% set
    """
    if not predicted_peptides or not known_epitopes:
        return {
            "recovery_rate": 0.0,
            "recovered": 0,
            "total_known": len(known_epitopes),
            "top_n_count": 0,
        }

    n = max(1, int(len(predicted_peptides) * top_n_percent / 100.0))
    top_set = predicted_peptides[:n]

    # Substring matching in both directions: a known epitope is "recovered"
    # if any top-N% prediction contains it OR is contained in it.  This is
    # intentionally generous — it accommodates the length mismatch between
    # fixed-window predictions (e.g. 15-mers) and variable-length IEDB
    # epitopes (6–26 aa).  A short known epitope (say 9-mer) will match
    # any 15-mer that happens to contain it, even if the binding is driven
    # by a different part of that 15-mer.  Keep this in mind when
    # interpreting recovery rates — they are an upper-bound estimate.
    recovered = 0
    for epitope in known_epitopes:
        for pred in top_set:
            if epitope in pred or pred in epitope:
                recovered += 1
                break

    return {
        "recovery_rate": recovered / len(known_epitopes),
        "recovered": recovered,
        "total_known": len(known_epitopes),
        "top_n_count": n,
    }
