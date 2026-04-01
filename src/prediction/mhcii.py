"""
MHC Class II binding prediction via the IEDB Tools API.

Wraps the IEDB ``netmhciipan_el`` prediction endpoint to score
pre-scanned peptides against a panel of HLA Class II alleles.
Results are cached as TSV files to avoid redundant API calls.

Dependencies: ``requests``, ``pandas``.
"""

from __future__ import annotations

import hashlib
import io
import os
import time
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_IEDB_PREDICT_URL = (
    "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
)

# Maximum FASTA entries per API request (IEDB limit is 200; stay under).
_BATCH_SIZE = 100

# Reference panel of common HLA Class II alleles covering broad global
# population diversity.  DRB1 alleles are the most polymorphic and
# best-studied; DQ and DP add coverage breadth.
#
# DQ and DP alleles require the full alpha/beta chain specification
# for the IEDB prediction API (beta-chain-only names return HTTP 500).
HLA_PANEL: list[str] = [
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*04:05",
    "HLA-DRB1*07:01",
    "HLA-DRB1*09:01",
    "HLA-DRB1*11:01",
    "HLA-DRB1*13:01",
    "HLA-DRB1*15:01",
    "HLA-DQA1*01:01/DQB1*05:01",
    "HLA-DQA1*01:02/DQB1*06:02",
    "HLA-DPA1*01:03/DPB1*04:01",
]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _peptides_hash(peptides: list[str]) -> str:
    """First 8 hex chars of the MD5 of the newline-joined peptide string."""
    joined = "\n".join(peptides)
    return hashlib.md5(joined.encode()).hexdigest()[:8]


def _safe_allele_filename(allele: str) -> str:
    """Sanitise an allele name for use in a file path."""
    return allele.replace("*", "_").replace(":", "").replace("/", "_")


def _build_fasta(peptides: list[str]) -> str:
    """Format peptide strings as a FASTA block."""
    lines: list[str] = []
    for i, pep in enumerate(peptides):
        lines.append(f">pep{i}")
        lines.append(pep)
    return "\n".join(lines)


def _post_iedb(fasta: str, allele: str, method: str) -> str:
    """POST one batch to the IEDB prediction API and return raw TSV text.

    Raises :class:`ValueError` on HTTP errors or non-tabular responses.
    """
    response = requests.post(
        _IEDB_PREDICT_URL,
        data={
            "method": method,
            "sequence_text": fasta,
            "allele": allele,
            "length": "asis",
        },
        timeout=300,
    )

    if response.status_code != 200:
        raise ValueError(
            f"IEDB prediction failed for {allele}: "
            f"HTTP {response.status_code} — {response.text[:300]}"
        )

    text = response.text.strip()
    if not text or "\t" not in text:
        raise ValueError(
            f"IEDB returned an unexpected response for {allele}: "
            f"{text[:300]}"
        )

    return text


def _parse_tsv(text: str) -> pd.DataFrame:
    """Parse IEDB's tab-separated response into a DataFrame.

    Renames the ``rank`` column to ``percentile_rank`` for clarity.
    """
    df = pd.read_csv(io.StringIO(text), sep="\t")
    # Normalise column names — IEDB uses "peptide" for most methods but
    # some return "sequence" instead.
    if "peptide" not in df.columns and "sequence" in df.columns:
        df = df.rename(columns={"sequence": "peptide"})
    if "rank" in df.columns:
        df = df.rename(columns={"rank": "percentile_rank"})
    return df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def predict_binding(
    peptides: list[str],
    allele: str = "HLA-DRB1*04:01",
    cache_dir: str | os.PathLike[str] = "data/processed",
    method: str = "netmhciipan_el",
) -> pd.DataFrame:
    """Predict MHC-II binding for *peptides* against a single *allele*.

    Sends peptides to the IEDB prediction API in FASTA format with
    ``length=asis`` so each peptide is scored at its actual length.
    Requests are batched into groups of 100 to stay within API limits.

    Parameters
    ----------
    peptides:
        List of amino-acid strings to score.
    allele:
        HLA Class II allele name (e.g. ``"HLA-DRB1*04:01"``).
    cache_dir:
        Directory for cached TSV responses.
    method:
        IEDB prediction method (default ``"netmhciipan_el"``).

    Returns
    -------
    DataFrame with columns ``peptide``, ``allele``, ``percentile_rank``,
    sorted ascending by percentile_rank (lower = stronger binder).

    Raises
    ------
    ValueError
        If the IEDB API returns an HTTP error or unexpected response body.
    """
    if not peptides:
        return pd.DataFrame(columns=["peptide", "allele", "percentile_rank"])

    # --- Cache check -------------------------------------------------------
    pep_hash = _peptides_hash(peptides)
    safe_allele = _safe_allele_filename(allele)
    cache_path = Path(cache_dir) / f"{safe_allele}_{method}_{pep_hash}.tsv"

    if cache_path.exists():
        df = pd.read_csv(cache_path, sep="\t")
        return (
            df[["peptide", "allele", "percentile_rank"]]
            .sort_values("percentile_rank")
            .reset_index(drop=True)
        )

    # --- Batched API calls --------------------------------------------------
    all_frames: list[pd.DataFrame] = []

    for batch_start in range(0, len(peptides), _BATCH_SIZE):
        batch = peptides[batch_start : batch_start + _BATCH_SIZE]
        fasta = _build_fasta(batch)
        raw_text = _post_iedb(fasta, allele, method)
        batch_df = _parse_tsv(raw_text)
        all_frames.append(batch_df)

        # Polite rate-limiting between batches
        if batch_start + _BATCH_SIZE < len(peptides):
            time.sleep(1)

    combined = pd.concat(all_frames, ignore_index=True)

    # --- Cache the full response (all columns) -----------------------------
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(cache_path, sep="\t", index=False)

    return (
        combined[["peptide", "allele", "percentile_rank"]]
        .sort_values("percentile_rank")
        .reset_index(drop=True)
    )


def predict_all_alleles(
    peptides: list[str],
    allele_panel: list[str] | None = None,
    cache_dir: str | os.PathLike[str] = "data/processed",
) -> pd.DataFrame:
    """Run :func:`predict_binding` for every allele in the panel.

    Parameters
    ----------
    peptides:
        List of amino-acid strings to score.
    allele_panel:
        List of HLA allele names.  Defaults to :data:`HLA_PANEL`.
    cache_dir:
        Directory for cached TSV responses.

    Returns
    -------
    Single DataFrame with all allele predictions, columns ``peptide``,
    ``allele``, ``percentile_rank``, sorted ascending by percentile_rank.
    """
    if allele_panel is None:
        allele_panel = HLA_PANEL

    frames: list[pd.DataFrame] = []

    for i, allele in enumerate(allele_panel):
        df = predict_binding(peptides, allele=allele, cache_dir=cache_dir)
        frames.append(df)

        # Rate-limit between alleles
        if i < len(allele_panel) - 1:
            time.sleep(1)

    if not frames:
        return pd.DataFrame(columns=["peptide", "allele", "percentile_rank"])

    return (
        pd.concat(frames, ignore_index=True)
        .sort_values("percentile_rank")
        .reset_index(drop=True)
    )


def top_binders(
    df: pd.DataFrame,
    percentile_threshold: float = 10.0,
) -> pd.DataFrame:
    """Filter to rows where ``percentile_rank`` <= *percentile_threshold*.

    Parameters
    ----------
    df:
        Output of :func:`predict_binding` or :func:`predict_all_alleles`.
    percentile_threshold:
        Maximum percentile rank to keep (default 10.0, the conventional
        weak-binder cutoff).  Use 2.0 for strong binders only.

    Returns
    -------
    Filtered DataFrame sorted by percentile_rank ascending.
    """
    return (
        df[df["percentile_rank"] <= percentile_threshold]
        .sort_values("percentile_rank")
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from src.data.uniprot import fetch_sequence
    from src.prediction.scanner import filter_peptides, scan_sequence

    print("Smoke test: ITGB3 (P05106) → 15-mers → HLA-DRB1*04:01\n")

    # 1. Fetch the sequence
    entry = fetch_sequence("P05106")
    seq = entry["sequence"]
    print(f"Sequence: {entry['name']} ({len(seq)} aa)")

    # 2. Scan for 15-mers, take first 50
    all_peptides = filter_peptides(scan_sequence(seq), min_len=15, max_len=15)
    first_50 = [p["peptide"] for p in all_peptides[:50]]
    print(f"Peptides: {len(first_50)} x 15-mers (first 50 of {len(all_peptides)})")

    # 3. Predict binding
    print("Predicting binding for HLA-DRB1*04:01 ...")
    results = predict_binding(first_50, allele="HLA-DRB1*04:01")
    print(f"Results: {len(results)} predictions\n")

    # 4. Print top 10
    top10 = results.head(10)
    print(f"{'Rank':<6} {'Peptide':<20} {'Percentile':>10}")
    print("-" * 38)
    for rank, (_, row) in enumerate(top10.iterrows(), 1):
        print(f"{rank:<6} {row['peptide']:<20} {row['percentile_rank']:>10.1f}")

    # 5. Summary
    binders = top_binders(results, percentile_threshold=10.0)
    print(f"\nBinders (percentile <= 10): {len(binders)} / {len(results)}")
