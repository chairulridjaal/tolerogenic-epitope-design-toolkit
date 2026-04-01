"""
Build a local 9-mer index from the human reviewed proteome.

Downloads the Swiss-Prot (reviewed) human proteome from UniProt,
extracts every unique 9-mer, and saves them as a compressed pickle.
This index powers the JMX proxy (Criterion 7) — a fast check for
whether a peptide's 9-mer windows appear in human self-proteins.

Usage (one-time):
    python -m src.assembly.build_jmx_index

Output:
    data/models/human_9mers.pkl.gz
"""

from __future__ import annotations

import gzip
import pickle
import re
import sys
from pathlib import Path

import requests

_UNIPROT_STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=reviewed:true+AND+organism_id:9606"
    "&format=fasta"
)

_OUTPUT_PATH = Path("data/models/human_9mers.pkl.gz")
_FASTA_CACHE = Path("data/raw/human_proteome_reviewed.fasta")


def download_human_proteome(cache_path: Path = _FASTA_CACHE) -> str:
    """Download the human reviewed proteome FASTA from UniProt."""
    if cache_path.exists():
        print(f"Using cached proteome: {cache_path} ({cache_path.stat().st_size / 1e6:.1f} MB)")
        return cache_path.read_text()

    print("Downloading human reviewed proteome from UniProt ...")
    resp = requests.get(_UNIPROT_STREAM_URL, timeout=300)
    if resp.status_code != 200:
        raise ValueError(f"UniProt download failed: HTTP {resp.status_code}")

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(resp.text)
    print(f"Saved to {cache_path} ({len(resp.text) / 1e6:.1f} MB)")
    return resp.text


def extract_sequences(fasta_text: str) -> list[str]:
    """Extract protein sequences from FASTA text."""
    sequences: list[str] = []
    current: list[str] = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            # Keep only standard amino acid characters
            cleaned = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", line.upper())
            current.append(cleaned)

    if current:
        sequences.append("".join(current))

    return sequences


def build_9mer_set(sequences: list[str]) -> set[str]:
    """Extract all unique 9-mers from a list of protein sequences."""
    ninemer_set: set[str] = set()
    for seq in sequences:
        for i in range(len(seq) - 8):
            ninemer_set.add(seq[i : i + 9])
    return ninemer_set


def save_index(ninemer_set: set[str], output_path: Path = _OUTPUT_PATH) -> None:
    """Save the 9-mer set as a gzip-compressed pickle."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wb") as f:
        pickle.dump(ninemer_set, f, protocol=pickle.HIGHEST_PROTOCOL)
    size_mb = output_path.stat().st_size / 1e6
    print(f"Saved {len(ninemer_set):,} unique 9-mers to {output_path} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    print("JMX Proxy — Human Proteome 9-mer Index Builder")
    print("=" * 50)

    fasta = download_human_proteome()
    sequences = extract_sequences(fasta)
    print(f"Parsed {len(sequences):,} protein sequences")

    total_residues = sum(len(s) for s in sequences)
    print(f"Total residues: {total_residues:,}")

    print("Extracting 9-mers ...")
    ninemers = build_9mer_set(sequences)
    print(f"Unique 9-mers: {len(ninemers):,}")

    save_index(ninemers)
    print("\nDone. Run the scorer to use JMX proxy scores.")
