"""
ESM-2 protein language model embedding infrastructure.

Provides shared embedding computation for all scoring criteria that
benefit from sequence-level representations.  Uses the smallest ESM-2
model (esm2_t6_8M_UR50D, 320-dim) — sufficient for short peptides
(15-mers) and appropriate for the training data sizes in this pipeline.

The model is loaded once per process and cached in module memory.
Embeddings are cached to disk after first computation.
"""

from __future__ import annotations

import hashlib
import logging
import os
import pickle
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

_MODEL_NAME = "esm2_t6_8M_UR50D"
_REPR_LAYER = 6  # final layer for t6 model
_EMBED_DIM = 320

# Module-level cache
_model = None
_alphabet = None
_batch_converter = None


def load_esm_model():
    """Load ESM-2 model and alphabet, caching in module memory."""
    global _model, _alphabet, _batch_converter

    if _model is not None:
        return _model, _alphabet, _batch_converter

    try:
        import esm
        import torch

        _model, _alphabet = esm.pretrained.esm2_t6_8M_UR50D()
        _batch_converter = _alphabet.get_batch_converter()
        _model.eval()

        logger.info("Loaded ESM-2 model: %s (%d-dim)", _MODEL_NAME, _EMBED_DIM)
        return _model, _alphabet, _batch_converter

    except ImportError as exc:
        raise ImportError(
            "ESM-2 requires fair-esm and torch. "
            "Install: pip install --break-system-packages fair-esm"
        ) from exc


def compute_embeddings(
    peptides: list[str],
    batch_size: int = 32,
) -> dict[str, np.ndarray]:
    """Compute mean-pooled ESM-2 embeddings for a list of peptides.

    Parameters
    ----------
    peptides :
        List of amino acid strings.
    batch_size :
        Number of peptides per forward pass.

    Returns
    -------
    Dict mapping peptide string → numpy array of shape (320,).
    Peptides with non-standard amino acids are silently skipped.
    """
    import torch

    model, alphabet, batch_converter = load_esm_model()

    _STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

    results: dict[str, np.ndarray] = {}

    # Deduplicate and filter non-standard sequences
    unique_peps = [p for p in set(peptides) if all(c in _STANDARD_AA for c in p.upper())]
    skipped = len(set(peptides)) - len(unique_peps)
    if skipped:
        logger.warning("Skipped %d peptides with non-standard characters", skipped)

    for i in range(0, len(unique_peps), batch_size):
        batch_peps = unique_peps[i : i + batch_size]
        data = [(f"pep{j}", pep) for j, pep in enumerate(batch_peps)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)

        with torch.no_grad():
            output = model(batch_tokens, repr_layers=[_REPR_LAYER])

        representations = output["representations"][_REPR_LAYER]

        for j, pep in enumerate(batch_peps):
            # Mean pool over residue positions, skipping BOS (0) and EOS (-1)
            seq_len = len(pep)
            token_repr = representations[j, 1 : seq_len + 1]  # skip BOS, take seq_len tokens
            mean_repr = token_repr.mean(dim=0).numpy()
            results[pep] = mean_repr

    return results


def compute_and_cache(
    peptides: list[str],
    cache_dir: str | os.PathLike[str] = "data/models",
    batch_size: int = 32,
) -> dict[str, np.ndarray]:
    """Compute embeddings with disk caching.

    Cache file: ``{cache_dir}/esm2_embeddings_{hash}.pkl``
    where hash is the first 8 hex chars of the MD5 of the sorted
    peptide set.
    """
    cache_path = Path(cache_dir) / f"esm2_embeddings_{_peptides_hash(peptides)}.pkl"

    if cache_path.exists():
        try:
            cached = pickle.loads(cache_path.read_bytes())
            # Verify all requested peptides are in cache
            if all(p in cached for p in peptides):
                logger.info("Loaded %d embeddings from cache: %s", len(peptides), cache_path.name)
                return {p: cached[p] for p in peptides}
        except Exception:
            pass  # corrupt cache

    embeddings = compute_embeddings(peptides, batch_size=batch_size)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(pickle.dumps(embeddings))
    logger.info("Cached %d embeddings to %s", len(embeddings), cache_path.name)

    return embeddings


def _peptides_hash(peptides: list[str]) -> str:
    joined = "\n".join(sorted(set(peptides)))
    return hashlib.md5(joined.encode()).hexdigest()[:8]
