"""
Train a local IL-10 induction prediction model (Random Forest).

Reproduces the best model from Nagpal et al. (2017) "Computer-aided
designing of immunosuppressive peptides based on IL-10 inducing
potential" using the exact 73 features from Table S1:
  - 16 amino acid composition features
  - 57 dipeptide composition features

Reads:
    data/raw/S1.csv   — feature specification (73 features)
    data/raw/S4.csv   — training set (394 positive + 848 negative)
    data/raw/S5.csv   — validation set (461 positive + 621 negative)

Writes:
    data/models/il10_rf_model.pkl   — trained RandomForestClassifier
    data/models/il10_features.json  — the 73 feature names for reproducibility

Usage:
    python -m src.scoring.train_il10_model
"""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score


# ---------------------------------------------------------------------------
# Feature specification (from Table S1 of Nagpal et al. 2017)
# ---------------------------------------------------------------------------

def load_feature_spec(path: str | Path = "data/raw/S1.csv") -> tuple[list[str], list[str]]:
    """Load the 73 feature names: 16 AA + 57 dipeptide compositions."""
    s1 = pd.read_csv(path)
    aa_feats = s1["amino_acid_feature"].dropna().tolist()
    dp_feats = s1["dipeptide_feature"].dropna().tolist()
    return aa_feats, dp_feats


# ---------------------------------------------------------------------------
# Data loading (handles the split-format CSVs)
# ---------------------------------------------------------------------------

def _parse_split_csv(path: str | Path) -> pd.DataFrame:
    """Parse IL-10pred's two-section CSV format.

    The CSVs from the supplementary material have:
      Section 1: ``sequence,1`` lines (positive peptides)
      Section 2: bare sequence lines or ``sequence,0`` / ``,0`` lines
                 (negative peptides), separated by a blank line or
                 format change.
    """
    with open(path) as f:
        lines = [l.strip() for l in f.readlines()]

    positives: list[str] = []
    negatives: list[str] = []

    i = 1  # skip header
    # Phase 1: positive sequences (seq,1 format)
    while i < len(lines):
        line = lines[i]
        if not line:
            i += 1
            break
        if "," in line:
            parts = line.split(",")
            if len(parts) == 2 and parts[1] == "1":
                positives.append(parts[0])
                i += 1
                continue
        # Bare sequence or format change → negatives
        break
    else:
        # All lines consumed (no blank separator)
        pass

    # Phase 2: negative sequences
    while i < len(lines):
        line = lines[i]
        i += 1
        if not line or line == ",0":
            continue
        if "," in line:
            parts = line.split(",")
            negatives.append(parts[0])
        else:
            negatives.append(line)

    rows = [(s, 1) for s in positives] + [(s, 0) for s in negatives]
    return pd.DataFrame(rows, columns=["sequence", "label"])


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def extract_features(
    sequences: list[str],
    aa_feats: list[str],
    dp_feats: list[str],
) -> pd.DataFrame:
    """Compute the 73-feature vector for each peptide.

    16 amino acid composition features: for each AA in *aa_feats*,
    count occurrences / peptide length × 100.

    57 dipeptide composition features: for each dipeptide in *dp_feats*,
    count occurrences / (peptide length − 1) × 100.
    """
    rows: list[list[float]] = []

    for seq in sequences:
        seq_upper = seq.upper()
        length = len(seq_upper)
        feats: list[float] = []

        # AA composition (%)
        for aa in aa_feats:
            count = seq_upper.count(aa)
            feats.append((count / max(length, 1)) * 100)

        # Dipeptide composition (%)
        denom = max(length - 1, 1)
        for dp in dp_feats:
            count = sum(
                1 for j in range(length - 1) if seq_upper[j : j + 2] == dp
            )
            feats.append((count / denom) * 100)

        rows.append(feats)

    col_names = [f"AA_{a}" for a in aa_feats] + [f"DP_{d}" for d in dp_feats]
    return pd.DataFrame(rows, columns=col_names)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("IL-10pred local model training")
    print("=" * 50)

    # 1. Load feature spec
    aa_feats, dp_feats = load_feature_spec()
    print(f"Features: {len(aa_feats)} AA + {len(dp_feats)} dipeptide = {len(aa_feats) + len(dp_feats)} total")

    # 2. Load training data
    train_df = _parse_split_csv("data/raw/S4.csv")
    print(f"Training (S4): {len(train_df)} peptides ({(train_df['label']==1).sum()} pos, {(train_df['label']==0).sum()} neg)")

    # 3. Extract features
    X_train = extract_features(train_df["sequence"].tolist(), aa_feats, dp_feats)
    y_train = train_df["label"].values

    # 4. Train
    print("\nTraining RandomForestClassifier(n_estimators=500) ...")
    clf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)
    clf.fit(X_train, y_train)
    print("Done.")

    # 5. Save model + feature spec
    model_dir = Path("data/models")
    model_dir.mkdir(parents=True, exist_ok=True)

    model_path = model_dir / "il10_rf_model.pkl"
    joblib.dump(clf, model_path)
    print(f"Model saved to {model_path}")

    features_path = model_dir / "il10_features.json"
    features_path.write_text(json.dumps({"aa": aa_feats, "dp": dp_feats}, indent=2))
    print(f"Feature spec saved to {features_path}")

    # 6. Evaluate on S5
    print("\n" + "=" * 50)
    print("Validation on S5")
    print("=" * 50)

    val_df = _parse_split_csv("data/raw/S5.csv")
    print(f"Validation (S5): {len(val_df)} peptides ({(val_df['label']==1).sum()} pos, {(val_df['label']==0).sum()} neg)")

    X_val = extract_features(val_df["sequence"].tolist(), aa_feats, dp_feats)
    y_val = val_df["label"].values

    y_prob = clf.predict_proba(X_val)[:, 1]
    y_pred = (y_prob >= 0.5).astype(int)

    auc = roc_auc_score(y_val, y_prob)
    mcc = matthews_corrcoef(y_val, y_pred)
    acc = accuracy_score(y_val, y_pred)

    print(f"  AUC:      {auc:.4f}")
    print(f"  MCC:      {mcc:.4f}")
    print(f"  Accuracy: {acc:.4f}")
    print(f"  (Paper reports ~0.81 accuracy)")
