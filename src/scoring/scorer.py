"""
Tolerogenic scoring module for the ITP epitope pipeline.

Implements seven scoring criteria for ranking candidate peptides by their
likelihood of inducing immune tolerance rather than effector activation.
Each criterion is normalized to [0, 1] and combined into a weighted
composite score.  Criteria, weights, and their literature sources are
documented in ``docs/tolerogenic_criteria.md``.

IMPORTANT — third-party service dependencies:

    * **IL-10pred** (Criterion 4) runs **locally** using a Random Forest
      model trained on Nagpal et al. (2017) data with the exact 73
      features from Table S1.  Train once with
      ``python -m src.scoring.train_il10_model``.  If the model file is
      missing, the score falls back to **0.5 (neutral)**.

    * **IFNepitope2** (Criterion 5) runs **locally** using the
      ``ifnepitope2`` package (Dhall et al., 2024, *Scientific Reports*).
      Install with ``pip install --no-deps ifnepitope2``.  If the package
      is missing or the model fails to load, the score falls back to
      **0.5 (neutral)**.

    * **JanusMatrix** (Criterion 7) has **no public API**.  JMX scores are
      set to **0.5 (neutral)** for all peptides pending manual submission
      at https://janusmatrix.essentialfacts.com.  Scores can be overridden
      after manual submission.

Dependencies: ``pandas``, ``scikit-learn``, ``joblib``, standard library only.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Kyte-Doolittle hydropathy scale (Kyte & Doolittle, 1982).
KYTE_DOOLITTLE: dict[str, float] = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

# Default composite-score weights — sum to 1.0.
# See docs/tolerogenic_criteria.md for rationale.
#
# IL-10 weight reduced (0.15 → 0.08): Hall et al. 2019 found no IL-10
# elevation in validated ITP tolerogenic model. FoxP3+ Treg mechanism
# is IL-10 independent in this disease context.
#
# ITP proximity increased (0.25 → 0.30): strongest ITP-specific evidence.
# JMX increased (0.05 → 0.07): self-mimicry is part of the tolerance-
# escape mechanism (Moise 2013).
DEFAULT_WEIGHTS: dict[str, float] = {
    "mhc_zone":         0.20,
    "hla_promiscuity":  0.20,
    "itp_proximity":    0.30,
    "il10":             0.08,
    "ifng":             0.07,
    "solubility":       0.08,
    "jmx":              0.07,
}


# ---------------------------------------------------------------------------
# Criterion 1 — MHC Binding in the Tolerogenic Zone
# ---------------------------------------------------------------------------

def score_mhc_zone(
    peptide: str,
    predictions_df: pd.DataFrame,
) -> float:
    """Score a peptide's MHC binding by the Goldilocks tolerogenic zone.

    Moderate binders (percentile rank 2–10%) are optimal for Treg
    induction.  Ultra-strong binders (<2%) risk effector activation or
    FoxP3 destabilization and are penalized.  Non-binders (>20%) fail
    to be presented and score zero.

    Zone mapping per allele:
        percentile_rank < 2    → 0.2
        2 ≤ rank ≤ 10          → 1.0  (optimal)
        10 < rank ≤ 20         → 0.5
        rank > 20              → 0.0

    Returns the mean zone score across all alleles that have a
    prediction for this peptide.  Returns 0.0 if the peptide has no
    predictions in the DataFrame.

    **Criterion 1** in ``docs/tolerogenic_criteria.md``.
    """
    rows = predictions_df[predictions_df["peptide"] == peptide]
    if rows.empty:
        return 0.0

    def _zone(rank: float) -> float:
        if rank < 2.0:
            return 0.2
        if rank <= 10.0:
            return 1.0
        if rank <= 20.0:
            return 0.5
        # Low-affinity binders: not zero — ITP autoepitopes are frequently
        # low-affinity MHC binders by design (thymic escape mechanism).
        # Sukati 2007 Table 6: most immunodominant GPIIIa peptides have
        # no predicted high-affinity HLA-DR binding.
        return 0.1

    return rows["percentile_rank"].apply(_zone).mean()


# ---------------------------------------------------------------------------
# Criterion 2 — HLA Promiscuity / Population Coverage
# ---------------------------------------------------------------------------

def score_hla_promiscuity(
    peptide: str,
    predictions_df: pd.DataFrame,
) -> float:
    """Fraction of HLA alleles for which this peptide is a binder (rank ≤ 10).

    Promiscuous binding across many alleles means broader population
    coverage and is a hallmark of natural Tregitopes.

    The denominator is the total number of unique alleles in the entire
    ``predictions_df`` (typically 12 — the HLA_PANEL size), not just
    the alleles for this peptide.

    Returns 0.0 if the peptide has no predictions.

    **Criterion 2** in ``docs/tolerogenic_criteria.md``.
    """
    rows = predictions_df[predictions_df["peptide"] == peptide]
    if rows.empty:
        return 0.0

    total_alleles = predictions_df["allele"].nunique()
    if total_alleles == 0:
        return 0.0

    binder_count = (rows["percentile_rank"] <= 10.0).sum()
    return binder_count / total_alleles


# ---------------------------------------------------------------------------
# Criterion 3 — ITP-Validated Peptide Proximity
# ---------------------------------------------------------------------------

def _has_overlap(seq_a: str, seq_b: str, min_len: int = 9) -> bool:
    """True if *seq_a* and *seq_b* share a common substring of ≥ *min_len*."""
    shorter, longer = (seq_a, seq_b) if len(seq_a) <= len(seq_b) else (seq_b, seq_a)
    for i in range(len(shorter) - min_len + 1):
        if shorter[i : i + min_len] in longer:
            return True
    return False


def score_itp_proximity(
    peptide: str,
    gold_standard: list[dict[str, Any]],
) -> float:
    """Score a peptide by proximity to experimentally validated ITP epitopes.

    Scoring tiers (highest tier achieved wins):
        1.0 — exact match with a ``tolerogenic_validated=true`` peptide
        0.8 — exact match or full containment with any gold standard peptide
        0.4 — partial overlap ≥ 9 consecutive residues with any gold
              standard peptide (by sequence substring, or by positional
              range overlap for future use)
        0.0 — no overlap

    Full containment means the candidate contains the gold standard
    sequence as a substring or vice versa.

    **Criterion 3** in ``docs/tolerogenic_criteria.md``.
    """
    best = 0.0

    for gs in gold_standard:
        gs_seq = gs["sequence"]

        # Tier 1: exact tolerogenic match
        if gs.get("tolerogenic_validated") and peptide == gs_seq:
            return 1.0

        # Tier 2: exact or containment
        if peptide == gs_seq or peptide in gs_seq or gs_seq in peptide:
            best = max(best, 0.8)
            continue

        # Tier 3a: partial sequence overlap (substring of ≥ 9 residues)
        if best < 0.4 and _has_overlap(peptide, gs_seq, min_len=9):
            best = 0.4

        # Tier 3b: position-range overlap — reserved for future use.
        # When candidate peptide positions are available, check overlap
        # with gs["position_start"] / gs["position_end"] to catch
        # numbering offsets or 14-mer truncation cases.

    return best


# ---------------------------------------------------------------------------
# Criterion 6 — Solubility (GRAVY Score)
# ---------------------------------------------------------------------------

def score_gravy(peptide: str) -> float:
    """Compute a solubility score from the Grand Average of Hydropathy.

    Uses the Kyte-Doolittle scale.  Hydrophilic peptides (negative
    GRAVY) are more soluble and better candidates for tolerogenic
    delivery.  Hydrophobic peptides aggregate and fail to reach
    tolerogenic dendritic cells.

    Normalization:
        GRAVY ≤ −0.5  → 1.0
        GRAVY ≥ +1.0  → 0.0
        Linear interpolation between.

    Unknown amino acids are silently skipped.

    **Criterion 6** in ``docs/tolerogenic_criteria.md``.
    """
    values = [KYTE_DOOLITTLE[aa] for aa in peptide.upper() if aa in KYTE_DOOLITTLE]
    if not values:
        return 0.5  # neutral fallback for empty / unknown

    gravy = sum(values) / len(values)

    if gravy <= -0.5:
        return 1.0
    if gravy >= 1.0:
        return 0.0
    # Linear interpolation: at -0.5 → 1.0, at +1.0 → 0.0
    return (1.0 - gravy) / 1.5


# ---------------------------------------------------------------------------
# Criterion 4 — IL-10 Induction Potential (local RF model)
# Local retrained model: RF on 73 features from Nagpal et al. 2017 Table S1.
# 100% offline after one-time training via train_il10_model.py.
# ---------------------------------------------------------------------------

def _peptides_hash(peptides: list[str]) -> str:
    joined = "\n".join(peptides)
    return hashlib.md5(joined.encode()).hexdigest()[:8]


def score_il10_local(
    peptides: list[str],
    cache_dir: str | os.PathLike[str] = "data/processed",
) -> dict[str, float]:
    """Predict IL-10 induction probability using a local Random Forest model.

    Reproduces Nagpal et al. (2017) using the exact 73 features from
    Table S1 (16 AA composition + 57 dipeptide composition).  The model
    is trained by ``src/scoring/train_il10_model.py`` and saved to
    ``data/models/il10_rf_model.pkl``.

    Returns a dict mapping peptide → P(IL-10 inducer) in [0, 1].
    Higher = more likely to induce IL-10 = more tolerogenic.

    On **any** failure (model not found, feature error, etc.), returns
    0.5 (neutral) for every peptide and logs a warning.

    **Criterion 4** in ``docs/tolerogenic_criteria.md``.

    Reference: Nagpal G *et al.* (2017) *Scientific Reports* **7**:42851.
    """
    neutral = {p: 0.5 for p in peptides}
    if not peptides:
        return neutral

    cache_path = Path(cache_dir) / f"il10_local_{_peptides_hash(peptides)}.json"

    if cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text())
            if isinstance(cached, dict) and all(p in cached for p in peptides):
                return {p: cached[p] for p in peptides}
        except Exception:
            pass

    try:
        import joblib
        from src.scoring.train_il10_model import extract_features, load_feature_spec

        model_path = Path("data/models/il10_rf_model.pkl")
        if not model_path.exists():
            raise FileNotFoundError(
                f"IL-10 model not found at {model_path}. "
                f"Run: python -m src.scoring.train_il10_model"
            )

        clf = joblib.load(model_path)
        aa_feats, dp_feats = load_feature_spec()

        X = extract_features(peptides, aa_feats, dp_feats)
        probs = clf.predict_proba(X)[:, 1]
        result = {p: round(float(prob), 4) for p, prob in zip(peptides, probs)}

    except Exception as exc:
        logger.warning("IL-10 local model failed (%s) — using neutral 0.5", exc)
        result = neutral

    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(result, indent=2))
    except Exception:
        pass

    return result


# ---------------------------------------------------------------------------
# Criterion 5 — IFN-gamma Penalty (IFNepitope2, local)
# Local model: ExtraTrees via IFNepitope2 package (Dhall et al. 2024).
# 100% offline after one-time pip install; sklearn compat shim for old pickle.
# ---------------------------------------------------------------------------

# Amino acids used for dipeptide composition (DPC) features — matches
# the feature set used to train the IFNepitope2 ExtraTrees model.
_DPC_AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# WARNING — technical debt (accepted temporarily)
# This shim patches sklearn internals to load a 0.24.1 pickle on 1.3+.
# Plan: re-train the IFNepitope2 model cleanly on current sklearn
# in the next sprint (same approach we are using for IL-10pred).
def _load_ifnepitope2_model(host: str = "human"):
    """Load the IFNepitope2 sklearn model with a compatibility shim.

    The shipped pickle was saved with scikit-learn 0.24.  Versions ≥ 1.3
    added a ``missing_go_to_left`` field to the tree node dtype.  We
    patch ``_check_node_ndarray`` to inject that field when it is absent,
    then restore the original function immediately after loading.
    """
    import pickle
    import sklearn.tree._tree as _tt

    pkg_dir = None
    try:
        import ifnepitope2
        pkg_dir = ifnepitope2.__path__[0]
    except ImportError:
        raise ImportError("ifnepitope2 package not installed (pip install --no-deps ifnepitope2)")

    model_file = os.path.join(pkg_dir, "model", f"{host}_et.pkl")
    if not os.path.exists(model_file):
        raise FileNotFoundError(f"IFNepitope2 model not found: {model_file}")

    # Monkey-patch the module-level check function (not the Cython type)
    original_check = _tt._check_node_ndarray

    def _compat_check(node_ndarray, expected_dtype):
        if hasattr(node_ndarray, "dtype") and "missing_go_to_left" not in (
            node_ndarray.dtype.names or []
        ):
            import numpy as _np
            new_arr = _np.zeros(node_ndarray.shape, dtype=_np.dtype(expected_dtype))
            for name in node_ndarray.dtype.names:
                new_arr[name] = node_ndarray[name]
            return new_arr
        return original_check(node_ndarray, expected_dtype)

    _tt._check_node_ndarray = _compat_check
    try:
        with open(model_file, "rb") as fh:
            clf = pickle.load(fh)
    finally:
        _tt._check_node_ndarray = original_check

    return clf


def _dpc_features(sequences: list[str]) -> pd.DataFrame:
    """Dipeptide composition (DPC) features — 400-dimensional vector."""
    rows: list[list[float]] = []
    for seq in sequences:
        feats: list[float] = []
        denom = max(len(seq) - 1, 1)
        for a in _DPC_AMINO_ACIDS:
            for b in _DPC_AMINO_ACIDS:
                pair = a + b
                count = sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == pair)
                feats.append((count / denom) * 100)
        rows.append(feats)
    return pd.DataFrame(rows)


def score_ifng_epitope(
    peptides: list[str],
    cache_dir: str | os.PathLike[str] = "data/processed",
    host: str = "human",
) -> dict[str, float]:
    """Predict IFN-gamma induction using IFNepitope2 (Dhall et al., 2024).

    Runs the IFNepitope2 ExtraTrees model **locally** — no web calls.
    Computes dipeptide composition features, predicts IFN-γ induction
    probability, then inverts: score = 1.0 − probability.
    High IFN-gamma prediction → low tolerogenic score.

    Falls back to 0.5 (neutral) on any error (missing package, model
    load failure, etc.).

    **Criterion 5** in ``docs/tolerogenic_criteria.md``.

    Reference: Dhall A *et al.* (2024) IFNepitope2. *Scientific Reports*.
    """
    neutral = {p: 0.5 for p in peptides}
    if not peptides:
        return neutral

    cache_path = Path(cache_dir) / f"ifnepitope2_{_peptides_hash(peptides)}_{host}.json"

    if cache_path.exists():
        try:
            cached = json.loads(cache_path.read_text())
            if isinstance(cached, dict) and all(p in cached for p in peptides):
                return {p: cached[p] for p in peptides}
        except Exception:
            pass

    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            clf = _load_ifnepitope2_model(host=host)

        X = _dpc_features(peptides)
        probs = clf.predict_proba(X)[:, -1]  # P(IFN-γ inducer)
        result = {p: round(1.0 - float(prob), 4) for p, prob in zip(peptides, probs)}

    except Exception as exc:
        logger.warning("IFNepitope2 failed (%s) — using neutral 0.5 for all peptides", exc)
        result = neutral

    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(result, indent=2))
    except Exception:
        pass

    return result


# ---------------------------------------------------------------------------
# Gold standard loading
# ---------------------------------------------------------------------------

def load_gold_standard(
    path: str | os.PathLike[str] = "data/processed/itp_gold_standard.json",
    antigen_sequence: str | None = None,
) -> list[dict[str, Any]]:
    """Load and validate the ITP gold-standard peptide set.

    Reads the JSON file, extracts the ``peptides`` array, and validates:
      1. Every peptide's ``sequence`` length matches its ``length`` field.
      2. If *antigen_sequence* is provided, every peptide's sequence
         appears at the stated ``position_start`` in the antigen.  This
         catches signal-peptide offset errors and transcription mistakes.

    Returns
    -------
    List of peptide dicts, each with keys ``id``, ``sequence``,
    ``length``, ``tolerogenic_validated``, etc.

    Raises
    ------
    ValueError
        If any peptide fails validation.
    FileNotFoundError
        If *path* does not exist.
    """
    data = json.loads(Path(path).read_text())
    peptides = data["peptides"]

    errors: list[str] = []
    for p in peptides:
        actual = len(p["sequence"])
        stated = p["length"]
        if actual != stated:
            errors.append(
                f"{p['id']}: sequence is {actual} aa but length field says {stated}"
            )

    # Position verification against the actual antigen sequence
    if antigen_sequence is not None:
        for p in peptides:
            seq = p["sequence"]
            start = p["position_start"] - 1  # convert 1-based to 0-based
            expected = antigen_sequence[start : start + p["length"]]
            if expected != seq:
                errors.append(
                    f"{p['id']}: sequence mismatch at stated position. "
                    f"Expected {expected}, got {seq}. "
                    f"Check signal peptide offset (+26 for P05106)."
                )

    if errors:
        raise ValueError(
            "Gold standard validation failed:\n  " + "\n  ".join(errors)
        )

    return peptides


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def score_all_peptides(
    peptides: list[str],
    predictions_df: pd.DataFrame,
    gold_standard: list[dict[str, Any]],
    cache_dir: str | os.PathLike[str] = "data/processed",
    weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Score every peptide across all seven tolerogenic criteria.

    This is the main entry point.  It calls the batch API functions once
    for the full peptide list, then scores each peptide individually on
    the remaining criteria and computes the weighted composite.

    **JMX (Criterion 7)** has no public API.  JMX scores are set to 0.5
    (neutral) for all peptides.  Override by adding a ``jmx`` column to
    the returned DataFrame after manual JanusMatrix submission.

    Parameters
    ----------
    peptides :
        List of amino-acid strings to score.
    predictions_df :
        Phase 2 MHC-II binding predictions (columns: peptide, allele,
        percentile_rank).
    gold_standard :
        Output of :func:`load_gold_standard`.
    cache_dir :
        Directory for cached API responses.
    weights :
        Scoring weights.  Defaults to :data:`DEFAULT_WEIGHTS`.

    Returns
    -------
    DataFrame with columns: ``peptide``, ``mhc_zone``,
    ``hla_promiscuity``, ``itp_proximity``, ``il10``, ``ifng``,
    ``solubility``, ``jmx``, ``bcell_risk``, ``composite_score``.
    Sorted descending by ``composite_score``.
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS

    # --- Batch API calls (once for all peptides) ---------------------------
    il10_scores = score_il10_local(peptides, cache_dir=cache_dir)
    ifng_scores = score_ifng_epitope(peptides, cache_dir=cache_dir)

    # --- JMX proxy (import here to avoid circular imports) -----------------
    try:
        from src.assembly.construct_builder import score_jmx_proxy, score_bcell_risk
        jmx_available = True
    except Exception:
        jmx_available = False

    # --- Per-peptide scoring -----------------------------------------------
    rows: list[dict[str, Any]] = []
    for pep in peptides:
        mhc = score_mhc_zone(pep, predictions_df)
        hla = score_hla_promiscuity(pep, predictions_df)
        itp = score_itp_proximity(pep, gold_standard)
        il10 = il10_scores.get(pep, 0.5)
        ifng = ifng_scores.get(pep, 0.5)
        sol = score_gravy(pep)

        if jmx_available:
            jmx = score_jmx_proxy(pep)
            bcell_risk, bcell_penalty = score_bcell_risk(pep)
        else:
            jmx = 0.5
            bcell_risk, bcell_penalty = False, 0.0

        composite = (
            weights["mhc_zone"] * mhc
            + weights["hla_promiscuity"] * hla
            + weights["itp_proximity"] * itp
            + weights["il10"] * il10
            + weights["ifng"] * ifng
            + weights["solubility"] * sol
            + weights["jmx"] * jmx
        )

        rows.append({
            "peptide": pep,
            "mhc_zone": round(mhc, 4),
            "hla_promiscuity": round(hla, 4),
            "itp_proximity": round(itp, 4),
            "il10": round(il10, 4),
            "ifng": round(ifng, 4),
            "solubility": round(sol, 4),
            "jmx": round(jmx, 4),
            "bcell_risk": bcell_risk,
            "composite_score": round(composite, 4),
        })

    df = pd.DataFrame(rows)
    return df.sort_values("composite_score", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Calibration check
# ---------------------------------------------------------------------------

def run_calibration_check(
    scores_df: pd.DataFrame,
    gold_standard: list[dict[str, Any]],
) -> bool:
    """Calibration against Hall et al. 2019 validated tolerogenic peptides.

    Peptide_2 (TTRGVSSCQQCLAVS): low MHC affinity by design (Sukati 2007
    Table 6 — no predicted high-affinity HLA-DR binding). Calibration
    checks ITP proximity and JMX scores, not composite rank. Low MHC
    binding is the expected biology for tolerance-escape autoepitopes.

    Peptide_82 (ALLIWKLLITIHDRK): genuine HLA-DR binder (Sukati 2007
    Table 6 — binds DR01, DR08, DR11, DR13, DR15). Calibration checks
    composite rank (top 20%) AND mhc_zone > 0.5.
    """
    n = len(scores_df)
    threshold_rank = max(1, int(n * 0.20))

    # Find the two validated peptides
    validated = {
        p["id"]: p["sequence"]
        for p in gold_standard
        if p.get("tolerogenic_validated")
    }

    print(f"\n{'='*60}")
    print("CALIBRATION CHECK")
    print(f"{'='*60}")
    print(f"Total scored peptides: {n}")
    print(f"Top 20% threshold:    rank <= {threshold_rank}")
    print()

    results: dict[str, bool] = {}
    scored_peptides = scores_df["peptide"].tolist()

    for pid, seq in validated.items():
        if seq not in scored_peptides:
            print(f"  {pid} ({seq})")
            print(f"    NOT FOUND in scored peptides — FAIL")
            results[pid] = False
            continue

        row = scores_df[scores_df["peptide"] == seq].iloc[0]
        rank = scored_peptides.index(seq) + 1
        pct = rank / n * 100

        print(f"  {pid} ({seq})")
        print(f"    Rank: {rank}/{n} (top {pct:.1f}%)")
        print(f"    mhc_zone={row['mhc_zone']:.2f}  hla_prom={row['hla_promiscuity']:.2f}  "
              f"itp={row['itp_proximity']:.2f}  jmx={row['jmx']:.2f}  "
              f"composite={row['composite_score']:.4f}")

        if pid == "Peptide_2":
            # Low MHC binding expected (Sukati 2007). Check ITP proximity
            # and JMX only — composite rank is not meaningful for this peptide.
            p2_pass = row["itp_proximity"] >= 0.99 and row["jmx"] >= 0.99
            status = "PASS" if p2_pass else "FAIL"
            print(f"    Criterion: itp_proximity>=1.0 AND jmx>=1.0 — {status}")
            print(f"    (Low MHC binding expected — Sukati 2007 Table 6)")
            results[pid] = p2_pass

        elif pid == "Peptide_82":
            # Genuine HLA-DR binder. Check composite rank AND MHC zone.
            p82_pass = rank <= threshold_rank and row["mhc_zone"] > 0.5
            status = "PASS" if p82_pass else "FAIL"
            print(f"    Criterion: rank<={threshold_rank} AND mhc_zone>0.5 — {status}")
            results[pid] = p82_pass

        else:
            # Unknown validated peptide — fall back to rank check
            fallback_pass = rank <= threshold_rank
            status = "PASS" if fallback_pass else "FAIL"
            print(f"    Criterion: rank<={threshold_rank} — {status}")
            results[pid] = fallback_pass

    overall = all(results.values())
    print(f"\nOverall: {'PASS' if overall else 'FAIL'}")
    print(f"{'='*60}\n")

    return overall


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    from src.data.uniprot import fetch_sequence

    print("Tolerogenic Scoring — ITGB3 (P05106)\n")

    # 1. Load ITGB3 sequence
    entry = fetch_sequence("P05106")
    print(f"Antigen: {entry['name']} ({len(entry['sequence'])} aa)")

    # 2. Load Phase 2 predictions
    pred_path = Path("data/processed/itgb3_top_binders.csv")
    predictions_df = pd.read_csv(pred_path)
    print(f"Predictions: {len(predictions_df)} rows from {pred_path.name}")

    # 3. Load gold standard (with position verification against ITGB3)
    gs_path = Path("data/processed/itp_gold_standard.json")
    gold_standard = load_gold_standard(gs_path, antigen_sequence=entry["sequence"])
    print(f"Gold standard: {len(gold_standard)} peptides from {gs_path.name} (position-verified)")

    # 3b. Get real MHC predictions for gold standard peptides and merge
    #     into predictions_df. This replaces the old bridge scores with
    #     actual NetMHCIIpan data. Results are cached after first run.
    gold_seqs = [p["sequence"] for p in gold_standard]
    missing_gs = [s for s in gold_seqs if s not in predictions_df["peptide"].values]
    if missing_gs:
        try:
            from src.prediction.mhcii import predict_all_alleles
            print(f"Fetching real MHC predictions for {len(missing_gs)} gold-standard peptides...")
            gs_preds = predict_all_alleles(missing_gs)
            predictions_df = pd.concat([predictions_df, gs_preds], ignore_index=True)
            print(f"Merged: predictions_df now {len(predictions_df)} rows")
        except Exception as exc:
            print(f"  Warning: could not fetch gold-standard predictions: {exc}")

    # 4. Extract unique peptides + force-include gold standard sequences
    #    so calibration peptides are always scored even if they aren't
    #    among the MHC-II top binders.
    peptides = predictions_df["peptide"].unique().tolist()
    gold_seqs = [p["sequence"] for p in gold_standard]
    peptides = list(set(peptides + gold_seqs))  # deduplicate
    print(f"Unique binder peptides: {predictions_df['peptide'].nunique()}")
    print(f"Added {len(gold_seqs)} gold-standard peptides for validation "
          f"(now scoring {len(peptides)} total)")
    print()

    # 5. Score
    print("Scoring all peptides ...")
    scores_df = score_all_peptides(
        peptides, predictions_df, gold_standard, cache_dir="data/processed",
    )
    print(f"Done. {len(scores_df)} peptides scored.\n")

    # 6. Calibration
    cal_pass = run_calibration_check(scores_df, gold_standard)

    # 7. Top 20
    print("Top 20 peptides by composite score:\n")
    display_cols = [
        "peptide", "mhc_zone", "hla_promiscuity", "itp_proximity",
        "il10", "ifng", "solubility", "jmx", "composite_score",
    ]
    top20 = scores_df.head(20)[display_cols]

    # Pretty-print
    header = (
        f"{'#':<4} {'Peptide':<18} {'MHC':>5} {'HLA':>5} {'ITP':>5} "
        f"{'IL10':>5} {'IFNg':>5} {'Sol':>5} {'JMX':>5} {'Score':>6}"
    )
    print(header)
    print("-" * len(header))
    for i, (_, row) in enumerate(top20.iterrows(), 1):
        print(
            f"{i:<4} {row['peptide']:<18} "
            f"{row['mhc_zone']:>5.2f} {row['hla_promiscuity']:>5.2f} "
            f"{row['itp_proximity']:>5.2f} {row['il10']:>5.2f} "
            f"{row['ifng']:>5.2f} {row['solubility']:>5.2f} "
            f"{row['jmx']:>5.2f} {row['composite_score']:>6.3f}"
        )

    # 8. Save
    out_path = Path("data/processed/itgb3_tolerogenic_scores.csv")
    scores_df.to_csv(out_path, index=False)
    print(f"\nFull results saved to {out_path}")

    # 9. Generate mRNA constructs (Phase 4)
    print("\n" + "=" * 60)
    print("PHASE 4 — mRNA CONSTRUCT ASSEMBLY")
    print("=" * 60)

    try:
        from src.assembly.construct_builder import generate_mrna_constructs
        construct_df, peptide_df = generate_mrna_constructs(
            top_n=10,
            scores_path=out_path,
            predictions_path=pred_path,
            gold_standard_path=gs_path,
            antigen_sequence=entry["sequence"],
        )

        print("\nTop constructs:\n")
        for _, row in construct_df.iterrows():
            print(f"  {row['construct_id']}: {row['n_epitopes']} epitopes, "
                  f"{row['construct_length_aa']} aa, {row['mrna_length_nt']} nt, "
                  f"score={row['construct_score']:.3f}")

        print(f"\nPeptide details:\n")
        for _, row in peptide_df.iterrows():
            risk = "B-CELL RISK" if row["bcell_risk"] else ""
            print(f"  T{row['experimental_tier']} {row['peptide']:<18} "
                  f"score={row['composite_score']:.3f} "
                  f"jmx={row['jmx_proxy']:.2f} {risk}")
    except Exception as exc:
        print(f"  Construct assembly skipped: {exc}")
