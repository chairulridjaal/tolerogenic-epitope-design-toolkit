"""
Multi-epitope mRNA construct assembly for the ITP tolerogenic pipeline.

Takes top-ranked tolerogenic peptides from Phase 3 scoring and assembles
them into candidate mRNA vaccine constructs with:
  - Flexible linker joining (GPGPG or AAY)
  - JMX proxy scoring (Criterion 7, replaces neutral 0.5)
  - B-cell epitope safety filtering (Parker hydrophilicity)
  - Construct-level scoring with bonuses and penalties
  - Human codon-optimized mRNA sequence generation
  - Experimental priority tiering

All operations are 100% offline.

Dependencies: ``pandas``, ``joblib`` (for cached 9-mer index), standard library.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import pickle
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_JMX_INDEX_PATH = Path("data/models/human_9mers.pkl.gz")

# Linker sequences — standard choices from vaccine literature.
LINKERS: dict[str, str] = {
    "GPGPG": "GPGPG",   # flexible, favours coil; most common in multi-epitope designs
    "AAY":   "AAY",      # rigid, favours processing by cathepsins
}

# Parker hydrophilicity scale (Parker et al., 1986).
# Higher values = more hydrophilic = more surface-exposed = higher B-cell risk.
PARKER_HYDROPHILICITY: dict[str, float] = {
    "A":  2.1, "R":  4.2, "N":  7.0, "D": 10.0, "C":  1.4,
    "E":  7.8, "Q":  6.0, "G":  5.7, "H":  2.1, "I": -8.0,
    "L": -9.2, "K":  5.7, "M": -4.2, "F": -9.2, "P":  2.1,
    "S":  6.5, "T":  5.2, "W": -10.0, "Y": -1.9, "V": -3.7,
}

# Human codon usage frequencies (per thousand) from the Kazusa Codon Usage
# Database for Homo sapiens (CDS: 93487 sequences).
# Each amino acid maps to a dict of {codon: frequency_per_thousand}.
HUMAN_CODON_TABLE: dict[str, dict[str, float]] = {
    "A": {"GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG":  7.4},
    "R": {"CGT":  4.5, "CGC": 10.4, "CGA":  6.2, "CGG": 11.4, "AGA": 12.2, "AGG": 12.0},
    "N": {"AAT": 17.0, "AAC": 19.1},
    "D": {"GAT": 21.8, "GAC": 25.1},
    "C": {"TGT": 10.6, "TGC": 12.6},
    "E": {"GAA": 29.0, "GAG": 39.6},
    "Q": {"CAA": 12.3, "CAG": 34.2},
    "G": {"GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5},
    "H": {"CAT": 10.9, "CAC": 15.1},
    "I": {"ATT": 16.0, "ATC": 20.8, "ATA":  7.5},
    "L": {"TTA":  7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6, "CTA":  7.2, "CTG": 39.6},
    "K": {"AAA": 24.4, "AAG": 31.9},
    "M": {"ATG": 22.0},
    "F": {"TTT": 17.6, "TTC": 20.3},
    "P": {"CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG":  6.9},
    "S": {"TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG":  4.4, "AGT": 12.1, "AGC": 19.5},
    "T": {"ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG":  6.1},
    "W": {"TGG": 13.2},
    "Y": {"TAT": 12.2, "TAC": 15.3},
    "V": {"GTT": 11.0, "GTC": 14.5, "GTA":  7.1, "GTG": 28.1},
}

# mRNA structural elements.
_KOZAK_5UTR = "GCCACCATG"  # Kozak consensus (includes start AUG in DNA form)
_STOP_CODON = "TGA"
# Minimal human beta-globin 3' UTR (first 50 nt of HBB 3'UTR, widely used).
_3UTR = (
    "GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAG"
)
_POLYA = "A" * 120


# ---------------------------------------------------------------------------
# JMX Proxy — Criterion 7
# ---------------------------------------------------------------------------

# Module-level cache so the 9-mer set is loaded only once per process.
_jmx_cache: set[str] | None = None


def _load_jmx_index(path: Path = _JMX_INDEX_PATH) -> set[str]:
    """Load the human 9-mer index, caching in module memory."""
    global _jmx_cache
    if _jmx_cache is not None:
        return _jmx_cache
    if not path.exists():
        raise FileNotFoundError(
            f"JMX index not found at {path}. "
            f"Run: python -m src.assembly.build_jmx_index"
        )
    with gzip.open(path, "rb") as f:
        _jmx_cache = pickle.load(f)
    return _jmx_cache


def score_jmx_proxy(peptide: str) -> float:
    """Score a peptide's self-similarity to the human proteome.

    Generates all 9-mer windows from the peptide and checks how many
    appear in the human reviewed proteome.  Higher fraction = more
    self-like = more likely to engage natural Tregs.

    Returns a score in [0, 1].  Returns 0.5 (neutral) if the JMX
    index has not been built.

    **Criterion 7** in ``docs/tolerogenic_criteria.md``.
    """
    try:
        index = _load_jmx_index()
    except FileNotFoundError:
        return 0.5

    pep = peptide.upper()
    if len(pep) < 9:
        return 0.0

    ninemers = [pep[i : i + 9] for i in range(len(pep) - 8)]
    hits = sum(1 for nm in ninemers if nm in index)
    return hits / len(ninemers)


# ---------------------------------------------------------------------------
# B-Cell Epitope Safety Filter
# ---------------------------------------------------------------------------

def score_bcell_risk(peptide: str, window: int = 7, threshold: float = 4.0) -> tuple[bool, float]:
    """Check if a peptide contains a likely linear B-cell epitope.

    Uses the Parker hydrophilicity scale (Parker et al., 1986) with a
    sliding window.  If any window's average hydrophilicity exceeds
    *threshold*, the peptide is flagged as a B-cell risk.

    A peptide that triggers antibody production against itself would
    undermine the tolerogenic goal.

    Returns
    -------
    (is_risky, penalty)
        ``is_risky`` is True if any window exceeds the threshold.
        ``penalty`` is -0.15 if risky, 0.0 otherwise.
    """
    pep = peptide.upper()
    if len(pep) < window:
        return False, 0.0

    for i in range(len(pep) - window + 1):
        win = pep[i : i + window]
        values = [PARKER_HYDROPHILICITY.get(aa, 0.0) for aa in win]
        if sum(values) / len(values) > threshold:
            return True, -0.15

    return False, 0.0


# ---------------------------------------------------------------------------
# Construct Assembly
# ---------------------------------------------------------------------------

def assemble_construct(peptides: list[str], linker: str = "GPGPG") -> str:
    """Join peptides with a linker sequence into a multi-epitope construct.

    Parameters
    ----------
    peptides :
        Ordered list of epitope amino-acid strings.
    linker :
        Linker sequence.  Default ``GPGPG`` (flexible).
        Alternative: ``AAY`` (rigid, favours cathepsin processing).

    Returns
    -------
    The full multi-epitope amino-acid string.
    """
    return linker.join(peptides)


def detect_junction_epitopes(
    construct: str,
    linker: str,
    known_binders: set[str] | None = None,
) -> list[str]:
    """Scan linker junctions for unintended strong MHC-II binders.

    Extracts 15-mers that span each junction (overlap with the linker
    region) and checks whether any appear in *known_binders* (the set
    of peptides with percentile_rank < 2 from Phase 2 predictions).

    Returns a list of flagged junction peptides.
    """
    if known_binders is None:
        return []

    flagged: list[str] = []

    # Find each linker occurrence and extract junction 15-mers
    linker_len = len(linker)
    start = 0
    while True:
        idx = construct.find(linker, start)
        if idx == -1:
            break

        # Junction region: from 7 residues before linker to 7 after
        region_start = max(0, idx - 7)
        region_end = min(len(construct), idx + linker_len + 7)
        region = construct[region_start:region_end]

        # Scan 15-mers in this region
        for j in range(len(region) - 14):
            kmer = region[j : j + 15]
            if kmer in known_binders:
                flagged.append(kmer)

        start = idx + 1

    return flagged


# ---------------------------------------------------------------------------
# Construct-Level Scoring
# ---------------------------------------------------------------------------

def score_full_construct(
    peptides: list[str],
    scores_df: pd.DataFrame,
    gold_standard: list[dict[str, Any]],
    linker: str = "GPGPG",
    known_strong_binders: set[str] | None = None,
) -> dict[str, Any]:
    """Score an assembled multi-epitope construct.

    Computes a construct-level score based on:
      - Average composite score of component peptides
      - +0.2 bonus if ≥ 3 epitopes from the same antigen
      - +0.1 bonus for spatial clustering (gold-standard positions
        within 50 residues of each other)
      - -0.15 penalty per junction with new strong binders

    Returns a dict with ``construct_score``, ``junction_flags``, and
    ``bonuses_applied``.
    """
    # Average composite score of components
    pep_scores = []
    for pep in peptides:
        row = scores_df[scores_df["peptide"] == pep]
        if not row.empty:
            pep_scores.append(row.iloc[0]["composite_score"])
    avg = sum(pep_scores) / max(len(pep_scores), 1)

    bonuses: list[str] = []

    # Bonus: ≥3 epitopes from same antigen (always true for single-antigen ITP)
    if len(peptides) >= 3:
        avg += 0.2
        bonuses.append("+0.2 multi-epitope (≥3 from same antigen)")

    # Bonus: spatial clustering — check if any component peptides overlap
    # gold-standard positions within 50 residues of each other
    gs_positions = [
        (gs.get("position_start", 0), gs.get("position_end", 0))
        for gs in gold_standard
    ]
    gs_peptide_seqs = {gs["sequence"] for gs in gold_standard}
    matching_gs = [p for p in peptides if p in gs_peptide_seqs]
    if len(matching_gs) >= 2:
        avg += 0.1
        bonuses.append("+0.1 spatial clustering (multiple gold-standard regions)")

    # Penalty: junction epitopes
    construct_seq = assemble_construct(peptides, linker)
    junction_flags = detect_junction_epitopes(construct_seq, linker, known_strong_binders)
    if junction_flags:
        penalty = 0.15 * len(junction_flags)
        avg -= penalty
        bonuses.append(f"-{penalty:.2f} junction epitope penalty ({len(junction_flags)} flagged)")

    return {
        "construct_score": max(0.0, min(1.0, avg)),
        "junction_flags": junction_flags,
        "bonuses_applied": bonuses,
    }


# ---------------------------------------------------------------------------
# Codon Optimization & mRNA Generation
# ---------------------------------------------------------------------------

def optimize_codons(
    aa_sequence: str,
    seed: int = 42,
    gc_max_window: float = 0.62,
    window_size: int = 18,
) -> str:
    """Convert an amino-acid sequence to a GC-balanced codon-optimized DNA string.

    Samples codons proportionally to their human usage frequency (Kazusa
    database) with a sliding-window GC constraint.  After selecting each
    codon, checks the last *window_size* codons; if GC content exceeds
    *gc_max_window*, the next codon is forced to the lowest-GC synonym.

    Deterministic for a given *seed*.  Targets overall GC of 50-60%.

    Parameters
    ----------
    aa_sequence : str
        Amino acid string.
    seed : int
        Random seed for reproducibility (default 42).
    gc_max_window : float
        Maximum GC fraction in the sliding window before forcing low-GC
        codons (default 0.62).
    window_size : int
        Number of codons in the sliding GC window (default 18, = 54 nt).
    """
    import random
    rng = random.Random(seed)

    dna: list[str] = []

    for aa in aa_sequence.upper():
        if aa not in HUMAN_CODON_TABLE:
            continue

        codons_freqs = HUMAN_CODON_TABLE[aa]
        codons = list(codons_freqs.keys())
        weights = list(codons_freqs.values())

        # Check sliding window GC — if too high, prefer low-GC codons
        if len(dna) >= window_size:
            recent = "".join(dna[-window_size:])
            gc_frac = sum(1 for nt in recent if nt in "GC") / len(recent)
            if gc_frac > gc_max_window:
                # Sort codons by GC count (ascending) and pick the lowest
                gc_counts = [sum(1 for nt in c if nt in "GC") for c in codons]
                min_gc = min(gc_counts)
                low_gc = [c for c, g in zip(codons, gc_counts) if g == min_gc]
                # Among lowest-GC codons, pick by frequency weight
                low_weights = [codons_freqs[c] for c in low_gc]
                chosen = rng.choices(low_gc, weights=low_weights, k=1)[0]
                dna.append(chosen)
                continue

        # Normal weighted sampling
        chosen = rng.choices(codons, weights=weights, k=1)[0]
        dna.append(chosen)

    return "".join(dna)


def build_mrna(aa_sequence: str) -> dict[str, Any]:
    """Generate a full mRNA sequence from an amino-acid construct.

    Structure: 5'UTR (Kozak + ATG) | CDS | Stop | 3'UTR (β-globin) | polyA(120)

    Returns
    -------
    Dict with keys:
        ``mrna_sequence`` — full DNA template (T not U; synthesizers use DNA)
        ``cds_start`` — 0-based start of the CDS
        ``cds_end`` — 0-based end of the CDS (exclusive)
        ``total_length`` — nucleotides including polyA
        ``gc_content`` — fraction of G+C in the CDS
        ``notes`` — manufacturing notes
    """
    cds = optimize_codons(aa_sequence)

    # Assemble: 5'UTR already contains ATG, CDS starts after it
    # Actually, Kozak ATG IS the start codon, so CDS = codons after ATG
    # We include ATG in Kozak, then the rest of the CDS (skipping the
    # first Met if the construct doesn't already start with M).
    full_dna = _KOZAK_5UTR + cds + _STOP_CODON + _3UTR + _POLYA

    cds_start = len(_KOZAK_5UTR) - 3  # ATG position within Kozak
    cds_end = cds_start + 3 + len(cds) + 3  # ATG + CDS + stop

    gc_bases = sum(1 for nt in cds if nt in "GC")
    gc_content = gc_bases / max(len(cds), 1)

    return {
        "mrna_sequence": full_dna,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "total_length": len(full_dna),
        "gc_content": round(gc_content, 3),
        "notes": (
            "All U residues → 1-methylpseudouridine (m1Ψ) during synthesis; "
            "formulate in dsRNA-depleted LNP for spleen/liver targeting; "
            "5' cap: CleanCap AG or m7GpppAm"
        ),
    }


# ---------------------------------------------------------------------------
# Experimental Tier Assignment
# ---------------------------------------------------------------------------

def assign_experimental_tier(
    peptide: str,
    composite_score: float,
    itp_proximity: float,
    gold_standard: list[dict[str, Any]],
) -> int:
    """Assign a testing priority tier.

    Tier 1 (top priority): composite ≥ 0.65 AND exact match to a
        tolerogenic-validated gold-standard peptide.
    Tier 2: composite ≥ 0.55 AND itp_proximity > 0.
    Tier 3: everything else.
    """
    validated_seqs = {
        gs["sequence"] for gs in gold_standard
        if gs.get("tolerogenic_validated")
    }
    if composite_score >= 0.65 and peptide in validated_seqs:
        return 1
    if composite_score >= 0.55 and itp_proximity > 0:
        return 2
    return 3


def select_diverse_peptides(
    scores_df: pd.DataFrame,
    antigen_sequence: str,
    top_n: int = 10,
    min_distance: int = 10,
) -> list[str]:
    """Select top peptides with positional diversity.

    Iterates through peptides ranked by composite_score.  After selecting
    each peptide, excludes all remaining candidates whose start position
    in *antigen_sequence* is within *min_distance* residues.  This
    prevents overlapping 15-mers from the same region from dominating
    the construct.

    Peptides not found in the antigen sequence (e.g. gold-standard
    peptides from a different numbering) are always included — they
    can't overlap anything.
    """
    selected: list[str] = []
    used_positions: list[int] = []

    for _, row in scores_df.iterrows():
        if len(selected) >= top_n:
            break

        pep = row["peptide"]
        pos = antigen_sequence.find(pep)

        # If peptide not found in antigen, always include it
        if pos == -1:
            selected.append(pep)
            continue

        # Check distance from all already-selected positions
        if any(abs(pos - used) < min_distance for used in used_positions):
            continue

        selected.append(pep)
        used_positions.append(pos)

    return selected


# ---------------------------------------------------------------------------
# Main Pipeline — Generate mRNA Constructs
# ---------------------------------------------------------------------------

def generate_mrna_constructs(
    top_n: int = 5,
    scores_path: str | os.PathLike[str] = "data/processed/itgb3_tolerogenic_scores.csv",
    predictions_path: str | os.PathLike[str] = "data/processed/itgb3_top_binders.csv",
    gold_standard_path: str | os.PathLike[str] = "data/processed/itp_gold_standard.json",
    output_path: str | os.PathLike[str] = "data/processed/itgb3_tolerogenic_mrna_constructs.csv",
    linker: str = "GPGPG",
    antigen_sequence: str | None = None,
) -> pd.DataFrame:
    """Generate multi-epitope mRNA vaccine constructs from scored peptides.

    Takes the top *top_n* peptides by composite score with positional
    diversity filtering (overlapping 15-mers from the same region are
    excluded), assembles them into a construct, and generates the
    codon-optimized mRNA sequence.

    Returns a DataFrame and saves to *output_path*.
    """
    # Load data
    scores_df = pd.read_csv(scores_path)
    predictions_df = pd.read_csv(predictions_path)

    gs_data = json.loads(Path(gold_standard_path).read_text())
    gold_standard = gs_data["peptides"]

    # Get strong binders for junction checking
    strong_binders = set(
        predictions_df[predictions_df["percentile_rank"] < 2.0]["peptide"].unique()
    )

    # Select top N peptides with positional diversity
    if antigen_sequence is not None:
        top_peptides = select_diverse_peptides(
            scores_df, antigen_sequence, top_n=top_n, min_distance=10,
        )
    else:
        # Fallback: no diversity filter
        top_peptides = scores_df.head(top_n)["peptide"].tolist()

    # Build a single construct from all top_n peptides
    construct_seq = assemble_construct(top_peptides, linker=linker)

    # Score the construct
    construct_result = score_full_construct(
        top_peptides, scores_df, gold_standard,
        linker=linker, known_strong_binders=strong_binders,
    )

    # Per-peptide analysis
    peptide_details: list[dict[str, Any]] = []
    for pep in top_peptides:
        row = scores_df[scores_df["peptide"] == pep]
        composite = row.iloc[0]["composite_score"] if not row.empty else 0.0
        itp_prox = row.iloc[0]["itp_proximity"] if not row.empty else 0.0

        jmx = score_jmx_proxy(pep)
        bcell_risk, bcell_penalty = score_bcell_risk(pep)
        tier = assign_experimental_tier(pep, composite, itp_prox, gold_standard)
        final_adj = composite + bcell_penalty

        # Notes
        notes_parts: list[str] = []
        for gs in gold_standard:
            if pep == gs["sequence"]:
                notes_parts.append(f"matches gold-standard {gs['id']}")
                if gs.get("tolerogenic_validated"):
                    notes_parts.append("(tolerogenic-validated)")
        if bcell_risk:
            notes_parts.append("B-cell risk flagged")

        peptide_details.append({
            "peptide": pep,
            "composite_score": round(composite, 4),
            "jmx_proxy": round(jmx, 4),
            "bcell_risk": bcell_risk,
            "bcell_penalty": round(bcell_penalty, 4),
            "final_adjusted_score": round(final_adj, 4),
            "experimental_tier": tier,
            "notes": "; ".join(notes_parts) if notes_parts else "",
        })

    # Build the mRNA
    mrna = build_mrna(construct_seq)

    # Create construct-level output row
    construct_rows = []
    construct_rows.append({
        "construct_id": f"ITP-ITGB3-{linker}-{top_n}ep",
        "peptide_list": " | ".join(top_peptides),
        "n_epitopes": len(top_peptides),
        "linker": linker,
        "full_construct_seq": construct_seq,
        "construct_length_aa": len(construct_seq),
        "mrna_length_nt": mrna["total_length"],
        "gc_content": mrna["gc_content"],
        "construct_score": round(construct_result["construct_score"], 4),
        "junction_flags": len(construct_result["junction_flags"]),
        "bonuses": "; ".join(construct_result["bonuses_applied"]),
        "mrna_sequence": mrna["mrna_sequence"],
        "manufacturing_notes": mrna["notes"],
    })

    # Also generate a variant with AAY linker
    if linker == "GPGPG":
        alt_construct_seq = assemble_construct(top_peptides, linker="AAY")
        alt_result = score_full_construct(
            top_peptides, scores_df, gold_standard,
            linker="AAY", known_strong_binders=strong_binders,
        )
        alt_mrna = build_mrna(alt_construct_seq)
        construct_rows.append({
            "construct_id": f"ITP-ITGB3-AAY-{top_n}ep",
            "peptide_list": " | ".join(top_peptides),
            "n_epitopes": len(top_peptides),
            "linker": "AAY",
            "full_construct_seq": alt_construct_seq,
            "construct_length_aa": len(alt_construct_seq),
            "mrna_length_nt": alt_mrna["total_length"],
            "gc_content": alt_mrna["gc_content"],
            "construct_score": round(alt_result["construct_score"], 4),
            "junction_flags": len(alt_result["junction_flags"]),
            "bonuses": "; ".join(alt_result["bonuses_applied"]),
            "mrna_sequence": alt_mrna["mrna_sequence"],
            "manufacturing_notes": alt_mrna["notes"],
        })

    construct_df = pd.DataFrame(construct_rows)
    peptide_df = pd.DataFrame(peptide_details)

    # Save both
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    construct_df.to_csv(out, index=False)

    peptide_out = out.parent / out.name.replace(".csv", "_peptide_detail.csv")
    peptide_df.to_csv(peptide_out, index=False)

    print(f"Saved {len(construct_df)} constructs to {out}")
    print(f"Saved {len(peptide_df)} peptide details to {peptide_out}")

    return construct_df, peptide_df


# ---------------------------------------------------------------------------
# CLI Smoke Test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    print("Phase 4 — Multi-Epitope mRNA Construct Assembly")
    print("=" * 55)

    # Load antigen sequence for positional diversity filtering
    from src.data.uniprot import fetch_sequence
    itgb3_seq = fetch_sequence("P05106")["sequence"]

    construct_df, peptide_df = generate_mrna_constructs(
        top_n=5, antigen_sequence=itgb3_seq,
    )

    print("\n--- Construct Summary ---\n")
    for _, row in construct_df.iterrows():
        print(f"Construct: {row['construct_id']}")
        print(f"  Epitopes: {row['n_epitopes']}")
        print(f"  Linker: {row['linker']}")
        print(f"  AA length: {row['construct_length_aa']}")
        print(f"  mRNA length: {row['mrna_length_nt']} nt")
        print(f"  GC content: {row['gc_content']:.1%}")
        print(f"  Construct score: {row['construct_score']:.4f}")
        print(f"  Junction flags: {row['junction_flags']}")
        print(f"  Bonuses: {row['bonuses']}")
        print()

    print("--- Peptide Details ---\n")
    header = f"{'Peptide':<20} {'Score':>6} {'JMX':>5} {'Bcell':>5} {'Final':>6} {'Tier':>4} Notes"
    print(header)
    print("-" * len(header))
    for _, row in peptide_df.iterrows():
        print(
            f"{row['peptide']:<20} "
            f"{row['composite_score']:>6.3f} "
            f"{row['jmx_proxy']:>5.2f} "
            f"{'YES' if row['bcell_risk'] else 'no':>5} "
            f"{row['final_adjusted_score']:>6.3f} "
            f"{row['experimental_tier']:>4} "
            f"{row['notes']}"
        )

    print("\n--- mRNA Preview (first construct, first 200 nt) ---\n")
    mrna = construct_df.iloc[0]["mrna_sequence"]
    for i in range(0, min(200, len(mrna)), 60):
        print(f"  {i+1:>4}  {mrna[i:i+60]}")
    print(f"  ...  ({len(mrna)} nt total)")
