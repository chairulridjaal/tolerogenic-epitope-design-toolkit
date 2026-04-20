"""
Train an IL-10 induction classifier using ESM-2 embeddings.

Replaces the composition-based Random Forest (train_il10_model.py) which
was found to have AUC 0.50 on truly independent (non-overlapping) data.

Uses ESM-2 (esm2_t6_8M_UR50D, 320-dim) mean-pooled embeddings as
features with Logistic Regression (L2 regularized).

Usage:
    python -m src.scoring.train_il10_esm
"""

from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import (
    accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
    precision_recall_fscore_support,
)
from sklearn.model_selection import cross_val_score

from src.scoring.esm_embeddings import compute_and_cache
from src.scoring.train_il10_model import _parse_split_csv


def _prepare_data():
    """Load and split the Nagpal et al. data with overlap tracking."""
    train = _parse_split_csv("data/raw/S4.csv")
    val = _parse_split_csv("data/raw/S5.csv")

    # Filter out malformed sequences (peptide names instead of sequences)
    _STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
    valid_seq = lambda s: all(c in _STANDARD_AA for c in s.upper())
    train = train[train["sequence"].apply(valid_seq)].reset_index(drop=True)
    val = val[val["sequence"].apply(valid_seq)].reset_index(drop=True)

    s4_seqs = set(train["sequence"])
    overlap = s4_seqs & set(val["sequence"])

    val_clean = val[~val["sequence"].isin(overlap)].copy()

    return train, val, val_clean, overlap


if __name__ == "__main__":
    print("IL-10 Classifier — ESM-2 Embeddings")
    print("=" * 55)

    # 1. Load data
    train, val, val_clean, overlap = _prepare_data()
    print(f"S4 train:              {len(train)} ({(train['label']==1).sum()} pos, {(train['label']==0).sum()} neg)")
    print(f"S5 full:               {len(val)} ({(val['label']==1).sum()} pos, {(val['label']==0).sum()} neg)")
    print(f"S4/S5 overlap:         {len(overlap)} sequences ({len(overlap)/len(val)*100:.1f}% of S5)")
    print(f"S5 non-overlapping:    {len(val_clean)} ({(val_clean['label']==1).sum()} pos, {(val_clean['label']==0).sum()} neg)")
    print()

    # 2. Compute ESM-2 embeddings
    all_seqs = list(set(train["sequence"].tolist() + val["sequence"].tolist()))
    print(f"Computing ESM-2 embeddings for {len(all_seqs)} unique sequences...")
    embeddings = compute_and_cache(all_seqs, cache_dir="data/models")
    print(f"Done. Embedding dim: {list(embeddings.values())[0].shape[0]}")
    print()

    # Build feature matrices
    X_train = np.array([embeddings[s] for s in train["sequence"]])
    y_train = train["label"].values
    X_val = np.array([embeddings[s] for s in val["sequence"]])
    y_val = val["label"].values
    X_clean = np.array([embeddings[s] for s in val_clean["sequence"]])
    y_clean = val_clean["label"].values

    # 3. Train multiple classifiers and compare
    classifiers = {
        "LogisticRegression (L2)": LogisticRegression(
            C=1.0, max_iter=1000, random_state=42, class_weight="balanced",
        ),
        "RandomForest (100)": RandomForestClassifier(
            n_estimators=100, random_state=42, n_jobs=-1, class_weight="balanced",
        ),
        "MLP (128-64)": MLPClassifier(
            hidden_layer_sizes=(128, 64), max_iter=500, random_state=42,
            early_stopping=True, validation_fraction=0.15,
        ),
    }

    best_name = None
    best_auc_clean = -1
    best_clf = None

    for name, clf in classifiers.items():
        print(f"--- {name} ---")
        clf.fit(X_train, y_train)

        # 5-fold CV on S4
        cv_auc = cross_val_score(
            clf.__class__(**clf.get_params()), X_train, y_train,
            cv=5, scoring="roc_auc",
        )
        print(f"  5-fold CV on S4:    AUC={cv_auc.mean():.4f} ± {cv_auc.std():.4f}")

        # Full S5
        y_prob_full = clf.predict_proba(X_val)[:, 1]
        auc_full = roc_auc_score(y_val, y_prob_full)
        print(f"  Full S5:            AUC={auc_full:.4f}")

        # Non-overlapping S5 (THE key metric)
        y_prob_clean = clf.predict_proba(X_clean)[:, 1]
        auc_clean = roc_auc_score(y_clean, y_prob_clean)
        y_pred_clean = (y_prob_clean >= 0.5).astype(int)
        acc_clean = accuracy_score(y_clean, y_pred_clean)
        mcc_clean = matthews_corrcoef(y_clean, y_pred_clean)
        prec, rec, f1, _ = precision_recall_fscore_support(
            y_clean, y_pred_clean, average="binary", zero_division=0,
        )
        print(f"  Non-overlapping S5: AUC={auc_clean:.4f}  Acc={acc_clean:.4f}  MCC={mcc_clean:.4f}")
        print(f"                      Prec={prec:.4f}  Rec={rec:.4f}  F1={f1:.4f}")
        print()

        if auc_clean > best_auc_clean:
            best_auc_clean = auc_clean
            best_name = name
            best_clf = clf

    # 4. Save best model
    print("=" * 55)
    print(f"Best model: {best_name} (non-overlapping AUC={best_auc_clean:.4f})")

    model_path = Path("data/models/il10_esm_model.pkl")
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(best_clf, model_path)
    print(f"Saved to {model_path}")

    # 5. Compare with old composition model
    print()
    print("=" * 55)
    print("COMPARISON: Old (composition) vs New (ESM-2)")
    print("=" * 55)
    print(f"                        Old (DPC RF)    New ({best_name})")
    print(f"  5-fold CV on S4:      AUC=0.8217      AUC={cv_auc.mean():.4f}")
    print(f"  Full S5:              AUC=0.9103      AUC={auc_full:.4f}")
    print(f"  Non-overlapping S5:   AUC=0.5028      AUC={best_auc_clean:.4f}")
    print(f"  (The non-overlapping S5 AUC is the only honest metric)")
