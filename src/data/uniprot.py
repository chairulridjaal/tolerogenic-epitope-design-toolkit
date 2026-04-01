"""
UniProt sequence retrieval for the ITP tolerogenic epitope pipeline.

Fetches protein sequences from the UniProt REST API in FASTA format,
caches them to disk, and parses them into plain sequence strings.
No external dependencies beyond ``requests``.
"""

from __future__ import annotations

import os
from pathlib import Path

import requests

from src.data.antigens import ITP_ANTIGENS

_UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def parse_fasta(fasta_text: str) -> dict[str, str]:
    """Parse raw FASTA text into its components.

    Parameters
    ----------
    fasta_text:
        A single-entry FASTA string (header line starting with ``>``
        followed by sequence lines).

    Returns
    -------
    dict with keys ``accession``, ``name``, and ``sequence``.

    Raises
    ------
    ValueError
        If *fasta_text* is empty or lacks a valid header line.
    """
    lines = fasta_text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("FASTA text is empty or missing a header line")

    header = lines[0][1:]  # drop the leading '>'

    # UniProt header format:  sp|P08514|ITA2B_HUMAN Description OS=...
    # We want the accession (second pipe-delimited field) and everything
    # after the third pipe up to " OS=" as the protein name.
    parts = header.split("|")
    if len(parts) >= 3:
        accession = parts[1].strip()
        # Name sits between the third pipe and the first " OS=" tag.
        rest = parts[2]
        os_idx = rest.find(" OS=")
        name = rest[rest.index(" ") + 1 : os_idx] if os_idx != -1 else rest.strip()
    else:
        # Non-standard header — use the whole line as the name.
        accession = header.split()[0]
        name = header

    sequence = "".join(line.strip() for line in lines[1:])
    return {"accession": accession, "name": name, "sequence": sequence}


def fetch_sequence(
    accession: str,
    cache_dir: str | os.PathLike[str] = "data/raw",
) -> dict[str, str]:
    """Fetch a single protein sequence from UniProt.

    If a cached FASTA file already exists at
    ``<cache_dir>/<accession>.fasta``, the network call is skipped and the
    cached file is read instead.

    Parameters
    ----------
    accession:
        UniProt accession number (e.g. ``"P08514"``).
    cache_dir:
        Directory where raw FASTA files are stored.

    Returns
    -------
    dict with keys ``accession``, ``name``, and ``sequence``.

    Raises
    ------
    ValueError
        If the HTTP response status is not 200.
    """
    cache_path = Path(cache_dir) / f"{accession}.fasta"

    if cache_path.exists():
        fasta_text = cache_path.read_text()
    else:
        url = _UNIPROT_FASTA_URL.format(accession=accession)
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            raise ValueError(
                f"UniProt request failed for {accession}: "
                f"HTTP {response.status_code} — {response.reason}"
            )

        fasta_text = response.text
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(fasta_text)

    return parse_fasta(fasta_text)


def fetch_all_itp_antigens(
    cache_dir: str | os.PathLike[str] = "data/raw",
) -> list[dict[str, str]]:
    """Fetch sequences for every antigen defined in ``ITP_ANTIGENS``.

    Parameters
    ----------
    cache_dir:
        Directory where raw FASTA files are stored.

    Returns
    -------
    List of dicts, each containing ``accession``, ``name``, ``sequence``,
    and ``complex`` (the platelet surface complex the protein belongs to).
    """
    results: list[dict[str, str]] = []
    for accession, meta in ITP_ANTIGENS.items():
        entry = fetch_sequence(accession, cache_dir=cache_dir)
        entry["complex"] = meta["complex"]
        entry["gene"] = meta["gene"]
        results.append(entry)
    return results


if __name__ == "__main__":
    antigens = fetch_all_itp_antigens()

    # Column widths
    acc_w = max(len(a["accession"]) for a in antigens)
    name_w = max(len(a["name"]) for a in antigens)
    cx_w = max(len(a["complex"]) for a in antigens)

    header = (
        f"{'Accession':<{acc_w}}  "
        f"{'Name':<{name_w}}  "
        f"{'Complex':<{cx_w}}  "
        f"{'Length':>6}"
    )
    print(header)
    print("-" * len(header))

    for a in antigens:
        print(
            f"{a['accession']:<{acc_w}}  "
            f"{a['name']:<{name_w}}  "
            f"{a['complex']:<{cx_w}}  "
            f"{len(a['sequence']):>6}"
        )
