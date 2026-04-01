"""
IEDB epitope retrieval for the ITP tolerogenic epitope pipeline.

Queries the IEDB Query API (PostgREST-based) for experimentally
characterised T-cell and B-cell epitopes mapped to ITP platelet
antigens.  Results are cached as raw JSON to avoid redundant API calls.

No external dependencies beyond ``requests``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import requests

from src.data.antigens import ITP_ANTIGENS

_IEDB_BASE = "https://query-api.iedb.org"

# PostgREST filter values ---------------------------------------------------
# parent_source_antigen_iri  = eq.UNIPROT:<accession>   (exact match)
# host_organism_iri_search   = cs.{NCBITaxon:9606}      (array contains)
_HUMAN_HOST_FILTER = "cs.{NCBITaxon:9606}"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_params(uniprot_id: str) -> dict[str, str]:
    """Build the shared query-string parameters for an IEDB search."""
    return {
        "parent_source_antigen_iri": f"eq.UNIPROT:{uniprot_id}",
        "host_organism_iri_search": _HUMAN_HOST_FILTER,
    }


def _fetch_iedb(
    endpoint: str,
    uniprot_id: str,
    cache_path: Path,
) -> list[dict[str, Any]]:
    """Fetch from an IEDB search endpoint, with disk caching.

    Returns the parsed JSON array (may be empty).
    """
    if cache_path.exists():
        return json.loads(cache_path.read_text())

    url = f"{_IEDB_BASE}/{endpoint}"
    params = _build_params(uniprot_id)
    response = requests.get(url, params=params, timeout=60)

    if response.status_code != 200:
        raise ValueError(
            f"IEDB request failed for {uniprot_id} ({endpoint}): "
            f"HTTP {response.status_code} — {response.reason}"
        )

    data: list[dict[str, Any]] = response.json()

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(data, indent=2))

    return data


def _extract_record(raw: dict[str, Any]) -> dict[str, str | None]:
    """Pull the fields we care about from one IEDB result row.

    IEDB field names (PostgREST view):
        linear_sequence   — peptide string, null for discontinuous epitopes
        mhc_allele_name   — e.g. "HLA-A*02:01", "human", or null
        assay_names        — pipe-delimited, e.g. "proliferation|3H-thymidine"
        disease_names      — JSON array or null, e.g. ["autoimmune thrombocytopenic purpura"]
        pubmed_id          — string, e.g. "17272505"
    """
    # disease_names is an array or null; flatten to a comma-separated string
    diseases = raw.get("disease_names")
    if isinstance(diseases, list):
        disease_str = ", ".join(str(d) for d in diseases)
    else:
        disease_str = None

    return {
        "epitope_sequence": raw.get("linear_sequence"),
        "mhc_allele": raw.get("mhc_allele_name"),
        "assay_type": raw.get("assay_names"),
        "disease": disease_str,
        "pubmed_id": raw.get("pubmed_id"),
    }


def _apply_disease_filter(
    records: list[dict[str, str | None]],
    disease_filter: str | None,
) -> list[dict[str, str | None]]:
    """Keep only records whose disease field contains *disease_filter* (case-insensitive).

    Common ITP filter values and their record counts (ITGB3 / P05106):

        "thrombocytopenic purpura"   → 88 T-cell records (matches the IEDB
                                       annotation "autoimmune thrombocytopenic
                                       purpura")
        "thrombocytopenia"           →  0 T-cell records (different IEDB term,
                                       only matches a handful of B-cell records
                                       on other antigens)

    Use ``"thrombocytopenic purpura"`` to select ITP-specific records.
    """
    if disease_filter is None:
        return records
    needle = disease_filter.lower()
    return [
        r for r in records
        if r["disease"] is not None and needle in r["disease"].lower()
    ]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fetch_tcell_epitopes(
    uniprot_id: str,
    cache_dir: str | os.PathLike[str] = "data/raw",
    disease_filter: str | None = None,
) -> list[dict[str, str | None]]:
    """Fetch T-cell epitopes for *uniprot_id* from IEDB.

    Parameters
    ----------
    uniprot_id:
        UniProt accession (e.g. ``"P05106"``).
    cache_dir:
        Directory for cached JSON responses.
    disease_filter:
        If provided, only return records whose disease field contains
        this substring (case-insensitive).  ``None`` returns all records.

    Returns
    -------
    List of dicts with keys ``epitope_sequence``, ``mhc_allele``,
    ``assay_type``, ``disease``, ``pubmed_id``.
    Empty list if IEDB has no results for this antigen.
    """
    cache_path = Path(cache_dir) / f"{uniprot_id}_tcell.json"
    raw_results = _fetch_iedb("tcell_search", uniprot_id, cache_path)
    records = [_extract_record(r) for r in raw_results]
    return _apply_disease_filter(records, disease_filter)


def fetch_bcell_epitopes(
    uniprot_id: str,
    cache_dir: str | os.PathLike[str] = "data/raw",
    disease_filter: str | None = None,
) -> list[dict[str, str | None]]:
    """Fetch B-cell epitopes for *uniprot_id* from IEDB.

    Parameters
    ----------
    uniprot_id:
        UniProt accession (e.g. ``"P05106"``).
    cache_dir:
        Directory for cached JSON responses.
    disease_filter:
        If provided, only return records whose disease field contains
        this substring (case-insensitive).  ``None`` returns all records.

    Returns
    -------
    List of dicts with keys ``epitope_sequence``, ``mhc_allele``,
    ``assay_type``, ``disease``, ``pubmed_id``.
    Empty list if IEDB has no results for this antigen.
    """
    cache_path = Path(cache_dir) / f"{uniprot_id}_bcell.json"
    raw_results = _fetch_iedb("bcell_search", uniprot_id, cache_path)
    records = [_extract_record(r) for r in raw_results]
    return _apply_disease_filter(records, disease_filter)


def fetch_all_itp_epitopes(
    cache_dir: str | os.PathLike[str] = "data/raw",
    disease_filter: str | None = None,
) -> list[dict[str, str | None]]:
    """Fetch T-cell and B-cell epitopes for every antigen in ``ITP_ANTIGENS``.

    Parameters
    ----------
    cache_dir:
        Directory for cached JSON responses.
    disease_filter:
        Passed through to the individual fetch functions.

    Returns
    -------
    Combined list of dicts.  Each dict has the standard five epitope
    fields plus ``epitope_type`` (``"T-cell"`` or ``"B-cell"``) and
    ``uniprot_id``.
    """
    combined: list[dict[str, str | None]] = []

    for accession in ITP_ANTIGENS:
        for fetch_fn, etype in [
            (fetch_tcell_epitopes, "T-cell"),
            (fetch_bcell_epitopes, "B-cell"),
        ]:
            records = fetch_fn(accession, cache_dir=cache_dir, disease_filter=disease_filter)
            for rec in records:
                rec["epitope_type"] = etype
                rec["uniprot_id"] = accession
            combined.extend(records)

    return combined


# ---------------------------------------------------------------------------
# CLI summary
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    cache = "data/raw"

    # The IEDB disease annotation for ITP is "autoimmune thrombocytopenic
    # purpura".  Filter on "thrombocytopenic purpura" (not just
    # "thrombocytopenia", which is a different IEDB term and matches
    # almost nothing for these antigens).
    itp_filter = "thrombocytopenic purpura"

    print(f"Fetching IEDB epitopes for all ITP antigens (filter: {itp_filter!r})...\n")

    # Gather counts per antigen
    rows: list[tuple[str, str, str, int, int]] = []
    for accession, meta in ITP_ANTIGENS.items():
        t = fetch_tcell_epitopes(accession, cache_dir=cache, disease_filter=itp_filter)
        b = fetch_bcell_epitopes(accession, cache_dir=cache, disease_filter=itp_filter)
        rows.append((accession, meta["gene"], meta["complex"], len(t), len(b)))

    # Print table
    gene_w = max(len(r[1]) for r in rows)
    cx_w = max(len(r[2]) for r in rows)

    header = (
        f"{'Accession':<10}  "
        f"{'Gene':<{gene_w}}  "
        f"{'Complex':<{cx_w}}  "
        f"{'T-cell':>7}  "
        f"{'B-cell':>7}"
    )
    print(header)
    print("-" * len(header))

    total_t = total_b = 0
    for acc, gene, cx, nt, nb in rows:
        total_t += nt
        total_b += nb
        print(
            f"{acc:<10}  "
            f"{gene:<{gene_w}}  "
            f"{cx:<{cx_w}}  "
            f"{nt:>7}  "
            f"{nb:>7}"
        )

    print("-" * len(header))
    print(
        f"{'TOTAL':<10}  "
        f"{'':<{gene_w}}  "
        f"{'':<{cx_w}}  "
        f"{total_t:>7}  "
        f"{total_b:>7}"
    )
