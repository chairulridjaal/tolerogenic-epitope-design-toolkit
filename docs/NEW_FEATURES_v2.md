# Tolerogenic Epitope Design Toolkit - New Features Since v1

## Overview
This document comprehensively catalogs ALL new features, analysis tools, and architectural changes added to the tolerogenic-epitope-design-toolkit since the original v1 pipeline. The toolkit has evolved from a single-disease (ITP) framework to a disease-agnostic architecture supporting Multiple Sclerosis (MS), with significant improvements to scoring criteria, investigation of previous models, and infrastructure for mRNA construct assembly.

---

## 1. CRITERION 4: Treg TCR-Contact Self-Similarity Scoring (Replaces IL-10)

### Location
- **File:** `src/scoring/scorer.py`
- **Function:** `score_treg_tcr_contact()` (lines 307-355)
- **Supporting:** `_build_tcr_motif_index()` (lines 266-304)

### What It Does
Replaces the original IL-10 induction model (which had AUC 0.50 on independent data) with a literature-grounded, no-ML-model approach based on the Tregitope concept (De Groot 2008) and JanusMatrix framework (Moise 2013).

**Biological Rationale:** In MHC-II presentation, the 9-mer binding core has:
- 4 anchor residues (P1, P4, P6, P9) buried in the MHC groove
- 5 TCR-contact residues (P2, P3, P5, P7, P8) facing the T-cell receptor

TCR-contact residues that are "self-like" (highly conserved in the human proteome) engage natural regulatory T cells (nTregs) selected in the thymus on self-peptides.

### Key Function Signatures
```python
def score_treg_tcr_contact(
    peptide: str,
    predictions_df: pd.DataFrame,
) -> float:
    """Returns percentile [0, 1] of TCR-motif frequency in human proteome."""
```

### Implementation Details

**TCR Position Constants** (line 259):
```python
_TCR_POSITIONS = [1, 2, 4, 6, 7]  # 0-indexed positions in 9-mer core
```

**Data Source:**
1. Loads human proteome 9-mer index from `data/models/human_9mers.pkl.gz`
2. Extracts 5 TCR-facing residues from each 9-mer
3. Builds frequency counter: 2,356,207 unique TCR motifs (73.6% of 20^5 = 3.2M possible)
4. For each candidate peptide:
   - Gets 9-mer core(s) from NetMHCIIpan predictions (`core_peptide` column)
   - Extracts TCR motif
   - Looks up frequency in human proteome counter
   - Converts to percentile (higher frequency → higher score)

**Return Value:** Percentile score in [0, 1]; returns 0.5 (neutral) if no core data available

### Validation Results

**Key Peptide Scores:**
- KEFAKFEEERARAKW (ITP novel #1): 0.95 (very common TCR contacts → strong Treg)
- RNLGELSRTTSEDNE (MS MBP30-44, validated): 0.90
- ALLIWKLLITIHDRK (ITP Peptide_82, validated): 0.69
- ENPVVHFFKNIVTPR (MS novel): 0.56 (was IL-10 highest at 0.88 — completely opposite!)

**Critical Discovery:** ENPVVHFFKNIVTPR had the HIGHEST IL-10 score (0.88) but the LOWEST Treg TCR-contact score (0.56). The TCR analysis says it has the most "foreign-like" TCR contacts. These are opposite conclusions from the same peptide. Since the IL-10 model doesn't generalize, the TCR analysis is more trustworthy.

### Comparison with JMX (Criterion 7)

**Problem with JMX:** For human self-proteins, JMX binary 9-mer matching gives 1.00 for every peptide (zero variance). Tested alternatives:
- Protein-count frequency
- Sequence similarity (≥7/9 identity)
- Neighbor counting
Result: All give zero or near-zero variance.

**Why TCR Score Works:** The 5-mer TCR motif space (20^5 = 3.2M) is denser than 9-mer space (20^9 = 512B), with 73.6% coverage. Frequency variation provides genuine discrimination.

### Weight Assignment

Originally (v1):
- IL-10: 0.08
- JMX: 0.07

After upgrade:
- **Treg TCR: 0.15** (absorbs former IL-10 + part of JMX weights)
- **JMX: 0.00** (zero discrimination, weight redistributed)

---

## 2. IL-10 DATA LEAKAGE INVESTIGATION

### Location
- **Investigation Log:** `docs/esm2_upgrade_log.md` (lines 14-45, 99-142)
- **Code Location:** `src/scoring/train_il10_model.py` (lines 1-150+)
- **Training Data:** `data/raw/S4.csv`, `data/raw/S5.csv` (Nagpal et al. 2017 supplementary data)

### What Was Investigated
The original IL-10 random forest model (from Nagpal et al. 2017) was achieving AUC 0.91 on the independent validation set (S5) but only AUC 0.50 on a non-overlapping subset of S5 (same as random chance).

### Discovery Process

**Testing three ESM-2-based classifiers** to see if better embeddings could rescue the model:

| Classifier | CV on S4 | Full S5 | Non-overlapping S5 |
|------------|----------|---------|-------------------|
| Old (DPC + RF) | AUC 0.82 | AUC 0.91 | **AUC 0.50** ← ISSUE |
| ESM-2 + LogReg | AUC 0.66 | AUC 0.65 | AUC 0.51 |
| ESM-2 + RF | AUC 0.65 | AUC 0.92 | AUC 0.46 |
| ESM-2 + MLP | AUC 0.63 | AUC 0.60 | AUC 0.51 |

### Key Finding

**ESM-2 does NOT rescue the IL-10 model.** All three classifiers achieve AUC ~0.50 on non-overlapping data. The problem is NOT the feature representation (ESM-2's 320-dim embeddings have vastly more information than 73 dipeptide features), but the **training data itself**:
- Only 394 positive examples (IL-10 inducers)
- Possibly insufficient diversity in the positive set
- IL-10 induction may not be purely sequence-determined — depends on MHC context, TCR repertoire, cytokine environment

### Decision

**IL-10 weight set to 0.00.** The criterion is retained in codebase as a placeholder for future models trained on larger, better-curated datasets.

### Reference Evidence

From Hall et al. (2019) ITP tolerogenic model validation:
- **No significant IL-10 elevation** in their validated Peptide_2 and Peptide_82 model
- Suppression mediated by CD4+CD25+FoxP3+ **Tregs (FoxP3+ mechanism is IL-10 independent)**
- This makes IL-10 a low-confidence criterion specifically for ITP

### Implications

The IL-10 model is not inherently broken — it's a genuine limitation of the training data. This finding is scientifically valuable because:
1. It reveals that IL-10 induction is likely multifactorial (sequence alone insufficient)
2. It explains why MHC-II binding + TCR contact biology (Criterion 4) is more predictive
3. It validates focusing on FoxP3+ Treg mechanisms for ITP

---

## 3. ESM-2 EMBEDDING INVESTIGATION

### Location
- **Code:** `src/scoring/esm_embeddings.py` (full module, lines 1-147)
- **Training:** `src/scoring/train_il10_esm.py` (referenced in docs)
- **Log:** `docs/esm2_upgrade_log.md` (lines 1-46, 99-142)

### What It Does

Provides infrastructure for computing protein language model embeddings using the smallest ESM-2 model (`esm2_t6_8M_UR50D`, 320-dim).

### Key Function Signatures

```python
def load_esm_model() -> tuple[Any, Any, Any]:
    """Load ESM-2 model and alphabet, cached in module memory."""

def compute_embeddings(
    peptides: list[str],
    batch_size: int = 32,
) -> dict[str, np.ndarray]:
    """Compute mean-pooled embeddings for peptide list.
    Returns dict mapping peptide → (320,) numpy array."""

def compute_and_cache(
    peptides: list[str],
    cache_dir: str | os.PathLike = "data/models",
    batch_size: int = 32,
) -> dict[str, np.ndarray]:
    """Compute embeddings with disk caching."""
```

### Implementation Details

**Model Choice:**
- Name: `esm2_t6_8M_UR50D`
- Dimension: 320-dim (smallest model, appropriate for 15-mers)
- Representation layer: 6 (final layer for t6 model)

**Processing:**
- Filters standard amino acids only (ACDEFGHIKLMNPQRSTVWY)
- Mean-pools over residue positions
- Skips BOS (beginning-of-sequence) and EOS (end-of-sequence) tokens
- Batch processing (default 32 peptides per forward pass)

**Caching:**
- Cache file: `data/models/esm2_embeddings_{hash}.pkl`
- Hash: First 8 hex chars of MD5 of sorted peptide set

### Important Note

While ESM-2 was investigated as a potential replacement for the IL-10 model, it did NOT improve performance. The embeddings are highly informative but cannot compensate for insufficient or biased training data. The module is retained for potential future use with better-curated IL-10 training datasets.

---

## 4. CATHEPSIN S PROCESSING LIKELIHOOD SCORING

### Location
- **File:** `src/scoring/processing.py` (full module, lines 1-132)
- **Integration:** Called in `scorer.py` line 671-706

### What It Does

Scores the likelihood that a 15-mer peptide is generated by cathepsin S digestion of the source protein in the endosomal/lysosomal compartment. Cathepsin S is the **dominant protease in MHC-II antigen processing** in professional APCs (dendritic cells, B cells).

### Key Function Signatures

```python
def score_cleavage_site(
    antigen: str,
    position: int,  # 0-based cut position
) -> float:
    """Score log-odds cleavage at this position."""

def score_processing(
    peptide: str,
    antigen_sequence: str,
) -> float | None:
    """Score likelihood peptide generated by cathepsin S.
    Returns sum of N-terminal + C-terminal cleavage scores."""

def score_processing_normalized(
    peptide: str,
    antigen_sequence: str,
    all_scores: list[float] | None = None,
) -> float:
    """Normalized processing score in [0, 1]."""
```

### Substrate Specificity (Position-Specific Matrices)

Literature-grounded (Choe et al. 2006, Biniossek et al. 2011, MEROPS C01.034):

**P2 (S2 pocket) — Dominant specificity determinant, weighted 3x:**
- Strong preference: Leu (score 2.5), Val (2.0), Phe (1.5), Ile (1.3)
- Strongly excluded: Asp (−3.0), Glu (−2.5), Arg (−2.5), Lys (−2.5)

**P1 (substrate residue):**
- Moderate selectivity, broadly tolerant except Pro (−2.0)
- Preferred: Gly (0.5), Lys (0.4), Thr (0.3)

**P1' (residue after cut):**
- Small/polar preferred: Ala (0.5), Ser (0.4), Glu (0.3)
- Excluded: Pro (−2.0), Trp (−0.5)

### Algorithm

1. Locate peptide in antigen sequence
2. Score N-terminal cleavage site (residues immediately before peptide start):
   - P2 position gets 3x weight (dominant signal)
   - P1 position gets 1x weight
   - P1' position gets 1x weight
3. Score C-terminal cleavage site (residues immediately after peptide end)
4. Sum both sites
5. Normalize to [0, 1] either by percentile across all peptides or fixed sigmoid

### Integration in Scoring Pipeline

In `score_all_peptides()` (line 669-707):
- Computes raw processing scores for all peptides
- Pre-computes all scores for percentile normalization
- Normalizes via `score_processing_normalized()` to [0, 1]
- Stores as `processing` column in output
- **NOT currently included in composite score** (experimental feature)

---

## 5. POPULATION COVERAGE COMPUTATION (Bui Formula)

### Location
- **File:** `src/assembly/population_coverage.py` (full module, lines 1-189)
- **Integration:** Called in `construct_builder.py` line 607-617

### What It Does

Computes what fraction of a given population can present at least one epitope from a multi-epitope construct, based on published HLA allele frequencies and peptide-allele binding data.

### Key Function Signatures

```python
def compute_coverage(
    epitope_alleles: dict[str, list[str]],
    population: str,  # "European", "East_Asian", "African", "South_Asian"
) -> float:
    """Fraction of population that can present at least one epitope."""

def compute_coverage_table(
    epitope_alleles: dict[str, list[str]],
) -> dict[str, float]:
    """Coverage across all populations."""

def get_epitope_alleles(
    peptides: list[str],
    predictions_df: pd.DataFrame,
    threshold: float = 10.0,
) -> dict[str, list[str]]:
    """Extract alleles each peptide binds."""

def format_coverage_report(
    epitope_alleles: dict[str, list[str]],
    construct_name: str = "",
) -> str:
    """Generate formatted coverage report."""
```

### Formula (Bui et al. 2006)

```
Coverage = 1 - ∏(1 - p_phenotype_i)
where product runs over all UNIQUE alleles that bind at least one epitope
p_phenotype_i = 1 - (1 - p_allele)^2  [Hardy-Weinberg]
```

### HLA Allele Frequencies Database

**DRB1 alleles (9 alleles, 7 populations):**
- HLA-DRB1*01:01, *03:01, *04:01, *04:05, *07:01, *09:01, *11:01, *13:01, *15:01
- European, East_Asian, African, South_Asian

**DQ alleles (2 alleles, proxy alpha/beta linkage):**
- HLA-DQA1*01:01/DQB1*05:01
- HLA-DQA1*01:02/DQB1*06:02

**DP allele (1 allele):**
- HLA-DPA1*01:03/DPB1*04:01

### Example Output

```
Population Coverage — CONSTRUCT
Epitopes: 10, covering 12 unique alleles

  European    : 95.2%
  East_Asian  : 87.3%
  African     : 91.8%
  South_Asian : 93.1%
```

---

## 6. DISEASE PROFILE SYSTEM (Generalization Architecture)

### Location
- **Code:** `src/data/disease_profile.py` (full module, lines 1-107)
- **ITP Profile:** `data/diseases/itp.json`
- **MS Profile:** `data/diseases/ms.json`

### What It Does

Enables disease-agnostic pipeline via JSON configuration files. Adding a new disease requires ONLY creating a new JSON profile and gold-standard file — no code changes.

### Function Signatures

```python
def load_disease_profile(
    disease_id: str,
    profiles_dir: str | os.PathLike = _PROFILES_DIR,
) -> dict[str, Any]:
    """Load and validate disease profile JSON."""

def get_antigens(profile: dict) -> list[dict]:
    """Return list of antigen dicts."""

def get_primary_antigen(profile: dict) -> str:
    """Return primary antigen UniProt ID."""

def get_iedb_filter(profile: dict) -> str:
    """Return IEDB disease filter string."""

def get_calibration_spec(profile: dict) -> dict:
    """Return calibration specification."""
```

### Profile JSON Schema

**Required Fields:**
- `disease_id`: Short ID (e.g., "itp", "ms")
- `disease_name`: Full name (e.g., "Immune Thrombocytopenic Purpura")
- `iedb_disease_filter`: Disease filter string for IEDB API
- `primary_antigen`: UniProt ID of main target
- `antigens`: List of all target antigens with metadata
- `gold_standard_path`: Path to gold-standard JSON
- `calibration`: Per-peptide calibration targets and thresholds

**ITP Profile Example:**
```json
{
  "disease_id": "itp",
  "disease_name": "Immune Thrombocytopenic Purpura",
  "primary_antigen": "P05106",
  "signal_peptide_length": 26,
  "antigens": [
    {"uniprot_id": "P08514", "gene": "ITGA2B", "name": "..."},
    {"uniprot_id": "P05106", "gene": "ITGB3", "name": "..."},
    ...
  ],
  "calibration": {
    "peptides": {
      "Peptide_2": {
        "check": "criteria",
        "criteria": {"itp_proximity": 0.99, "jmx": 0.99},
        "note": "..."
      }
    }
  }
}
```

**MS Profile Example:**
```json
{
  "disease_id": "ms",
  "disease_name": "Multiple Sclerosis",
  "primary_antigen": "P02686",
  "antigens": [
    {"uniprot_id": "P02686", "gene": "MBP", ...},
    {"uniprot_id": "P60201", "gene": "PLP1", ...},
    {"uniprot_id": "Q16653", "gene": "MOG", ...}
  ],
  "calibration": {
    "peptides": {
      "MBP83-99": {
        "check": "criteria",
        "criteria": {"itp_proximity": 0.40, "jmx": 0.50},
        "note": "..."
      }
    }
  }
}
```

### Calibration Specification

Defines per-peptide validation targets:

**"criteria" check:** Validates specific score thresholds
**"rank_and_mhc" check:** Validates rank percentile + MHC zone
**"rank" check:** Default — just validates rank percentile

---

## 7. UPDATED SCORING WEIGHTS & CRITERIA

### Location
- **File:** `src/scoring/scorer.py` (lines 43-69)
- **Documentation:** `docs/tolerogenic_criteria.md`

### Current Weights (v2+)

```python
DEFAULT_WEIGHTS: dict[str, float] = {
    "mhc_zone":         0.20,
    "hla_promiscuity":  0.20,
    "itp_proximity":    0.30,
    "treg_tcr":         0.15,  # NEW: Replaces IL-10
    "ifng":             0.07,
    "solubility":       0.08,
    "jmx":              0.00,  # Zero (no variance for self-proteins)
}
```

### Previous Weights (v1)

```python
DEFAULT_WEIGHTS: dict[str, float] = {
    "mhc_zone":         0.20,
    "hla_promiscuity":  0.20,
    "itp_proximity":    0.25,  # Was 0.25
    "il10":             0.15,  # WAS 0.15 (REMOVED)
    "ifng":             0.07,
    "solubility":       0.08,
    "jmx":              0.05,  # WAS 0.05
}
```

### Changes Explained

1. **IL-10 → Treg TCR (0.15 → 0.15):** Direct replacement — same weight but new, evidence-based criterion
2. **ITP proximity (0.25 → 0.30):** Increased — strongest ITP-specific evidence
3. **JMX (0.05 → 0.00):** Set to zero — binary 9-mer matching gives 1.00 for all human proteins (zero variance)

### Criterion Documentation

**Criterion 1 — MHC Zone:** Moderate binders (2–10% rank) optimal for Treg induction
**Criterion 2 — HLA Promiscuity:** Fraction of alleles where rank ≤ 10%
**Criterion 3 — Disease Proximity:** Substring overlap with validated epitopes
**Criterion 4 — Treg TCR Contact:** [NEW] Percentile frequency of TCR-facing motif in human proteome
**Criterion 5 — IFN-gamma Penalty:** 1.0 − P(IFN-γ inducer) via IFNepitope2
**Criterion 6 — Solubility:** GRAVY score normalization
**Criterion 7 — JMX Self-Similarity:** 9-mer overlap with human proteome (currently 0 weight)

---

## 8. NEW DISEASE PROFILES & VALIDATION

### Location
- **Files:** 
  - `data/diseases/itp.json`
  - `data/diseases/ms.json`
- **Gold Standards:**
  - `data/processed/itp_gold_standard.json`
  - `data/processed/ms_gold_standard.json`

### ITP Profile

**Disease:** Immune Thrombocytopenic Purpura
**Primary Antigen:** P05106 (ITGB3 — Integrin beta-3 / GPIIIa)
**Secondary Antigens:** P08514, P07359, P13224, P14770, P40197 (GPIIb complex, GPIb-IX-V complex)
**IEDB Filter:** "thrombocytopenic purpura"
**Signal Peptide Offset:** 26 amino acids (critical for position verification)

**Calibration Targets:**
- Peptide_2 (aa6–20): GDCNCTKDDSVMCIG — must achieve itp_proximity ≥ 0.99 + jmx ≥ 0.99
- Peptide_82 (aa711–725): NPIYKSAVTTVVNP — must rank in top 20% AND mhc_zone > 0.5

### MS Profile

**Disease:** Multiple Sclerosis
**Primary Antigen:** P02686 (MBP — Myelin Basic Protein)
**Secondary Antigens:** P60201 (PLP1), Q16653 (MOG)
**IEDB Filter:** "multiple sclerosis"

**Calibration Targets:**
- MBP83-99 (17-mer): Overlapping 15-mers must achieve itp_proximity ≥ 0.40 + jmx ≥ 0.50
- MBP30-44: Same criteria

**Note:** 17-mers not in 15-mer scan — calibration checks proximity scores of overlapping 15-mers

---

## 9. SCORING COMMAND-LINE INTERFACE

### Location
- **Entry Point:** `src/scoring/scorer.py` (lines 877-1031, `if __name__ == "__main__"`)
- **Usage:**
  ```bash
  python -m src.scoring.scorer --disease itp
  python -m src.scoring.scorer --disease ms
  ```

### Features

1. **Disease Profile Loading:** Loads all settings, antigens, calibration specs
2. **Antigen Sequence Fetching:** Retrieves from UniProt with length verification
3. **Predictions Loading:** Reads Phase 2 MHC-II predictions (CSV)
4. **Core Peptide Enrichment:** Merges core_peptide data from cached TSVs if missing
5. **Gold Standard Validation:** Position-verified against antigen sequence
6. **Batch Scoring:** All peptides scored on all criteria
7. **Calibration Check:** Validates against disease-specific targets
8. **Top 20 Report:** Formatted output of highest-scoring peptides
9. **Construct Assembly:** Triggers Phase 4 mRNA generation

### Output Files

- `data/processed/{disease_id}_tolerogenic_scores.csv` — Full scored peptide list
- `data/processed/{disease_id}_mrna_constructs.csv` — Multi-epitope constructs
- `data/processed/{disease_id}_mrna_constructs_peptide_detail.csv` — Per-peptide construct data

---

## 10. CONSTRUCT ASSEMBLY IMPROVEMENTS

### Location
- **File:** `src/assembly/construct_builder.py` (lines 1-699)

### Key New Features

**JMX Proxy Scoring (Criterion 7):**
- Function: `score_jmx_proxy(peptide: str) -> float` (lines 112-135)
- Generates all 9-mer windows from peptide
- Checks fraction found in human proteome index
- Returns [0, 1] score

**B-Cell Epitope Safety Filter:**
- Function: `score_bcell_risk(peptide, antigen_sequence) -> tuple[bool, float]` (lines 142-206)
- Uses Parker hydrophilicity scale
- NEW: Protein-relative thresholding (top 80% for source protein)
- Returns (is_risky, penalty: -0.15 if flagged)

**Junction Epitope Detection:**
- Function: `detect_junction_epitopes(construct, linker, known_binders)` (lines 231-270)
- Scans 15-mers spanning linker junctions
- Flags if any match known strong binders (rank < 2%)

**Construct-Level Scoring:**
- Function: `score_full_construct(peptides, scores_df, ...)` (lines 277-335)
- Average composite score of components
- +0.2 bonus: ≥3 epitopes from same antigen
- +0.1 bonus: Spatial clustering (gold-standard regions close together)
- −0.15 penalty: Per junction with new strong binders

**Codon Optimization:**
- Function: `optimize_codons(aa_sequence, seed=42, gc_max_window=0.62)` (lines 342-401)
- Human codon usage frequencies (Kazusa database)
- Sliding-window GC constraint (max 62% in 54nt windows)
- Deterministic seeded randomness for reproducibility

**mRNA Generation:**
- Function: `build_mrna(aa_sequence: str) -> dict` (lines 404-444)
- Structure: 5'UTR (Kozak + ATG) | CDS | Stop | 3'UTR (β-globin) | polyA(120)
- Manufacturing notes: m1Ψ, CleanCap, dsRNA-depleted LNP

**Experimental Tier Assignment:**
- Function: `assign_experimental_tier(peptide, scores, gold_standard) -> int` (lines 451-472)
- Tier 1: composite ≥ 0.65 + exact gold-standard match
- Tier 2: composite ≥ 0.55 + itp_proximity > 0
- Tier 3: Everything else

**Peptide Diversity Selection:**
- Function: `select_diverse_peptides(scores_df, antigen_seq, top_n=10)` (lines 475-515)
- Selects top peptides with positional diversity
- Excludes overlapping 15-mers within min_distance (default 10 aa)

---

## 11. NEW ANALYSIS DOCUMENTS

### Location
- `docs/esm2_upgrade_log.md` — Complete upgrade investigation log (16 KB)
- `docs/itp_prototype_v1.2_report.md` — ITP results with new scoring (11 KB)
- `docs/itp_prototype_v1.3_report.md` — Latest ITP results (12 KB)
- `docs/ms_prototype_v1.0_report.md` — MS baseline results (12 KB)
- `docs/ms_prototype_v1.1_report.md` — MS with Treg TCR scoring (11 KB)
- `docs/preprint_v2.md` — Updated preprint (26 KB)
- `docs/preprint_v3.md` — Latest preprint version (24 KB)
- `docs/weight_sensitivity_analysis.md` — Weight tuning analysis (5.7 KB)

### Key Findings in Reports

**ITP Prototype v1.3:**
- Peptide_82: Rank 2/257 (top 0.8%) ✓
- Peptide_2: Rank 53/257 (top 20.6%) ✓
- Both pass calibration with new weights
- 156/257 peptides shifted >5 rank positions (meaningful reordering)

**MS Prototype v1.1:**
- MBP131-145: Rank 1/113 ✓
- MBP30-44: Rank 2/113 ✓
- Both validated tolerogenic candidates from ATX-MS-1467
- 41/113 peptides shifted >5 positions

---

## 12. WEIGHT SENSITIVITY ANALYSIS

### Location
- **File:** `docs/weight_sensitivity_analysis.md`

### Analysis Scope

Tests calibration robustness to weight variations:
- Treg TCR: ±50% variation (±0.075)
- JMX: 0.00 to 0.05 (testing if it should be zero)
- Solubility: ±30% (0.056 to 0.104)

### Key Finding

**Treg TCR weight robust:** Calibration passes across the tested range. The criterion provides consistent ranking signal independent of ±50% weight adjustments.

---

## SUMMARY TABLE

| Feature | File | Type | Status |
|---------|------|------|--------|
| Treg TCR Contact (C4) | scorer.py | New criterion | Production |
| IL-10 data leakage | esm2_upgrade_log.md | Investigation | Documented |
| ESM-2 embeddings | esm_embeddings.py | Infrastructure | Available |
| Cathepsin S processing | processing.py | New scoring | Integrated |
| Population coverage (Bui) | population_coverage.py | New module | Production |
| Disease profiles | disease_profile.py | Architecture | Production |
| ITP profile | data/diseases/itp.json | Configuration | Active |
| MS profile | data/diseases/ms.json | Configuration | Active |
| Updated weights | scorer.py (61-69) | Parameters | Production |
| Construct assembly | construct_builder.py | Enhanced | Production |
| CLI interface | scorer.py (877+) | Interface | Production |
| Analysis documents | docs/ | Reports | Reference |

---

## VERIFICATION CHECKLIST

### ITP Pipeline (v1.3)
- ✅ Calibration: Peptide_82 rank 2, Peptide_2 rank 53
- ✅ Treg TCR variance: 0.277 − 0.999 (good discrimination)
- ✅ 156/257 peptides reranked (meaningful changes)
- ✅ Top 10 stable (9/10 overlap maintained)

### MS Pipeline (v1.1)
- ✅ Calibration: MBP131-145 rank 1, MBP30-44 rank 2
- ✅ Both validated tolerogenic candidates
- ✅ 41/113 peptides reranked
- ✅ Top 10 stable (8/10 overlap maintained)

### Integration Tests
- ✅ Disease profiles load and validate
- ✅ Position verification works (signal peptide offset)
- ✅ Construct assembly generates mRNA
- ✅ Population coverage computes correctly
- ✅ Codon optimization produces GC-balanced sequences

---

## ARCHITECTURE DIAGRAM

```
┌─────────────────────────────────────────────────────────┐
│           DISEASE-AGNOSTIC PIPELINE v2                  │
├─────────────────────────────────────────────────────────┤
│                                                           │
│  Disease Profile (JSON)                                  │
│  ├─ disease_id, primary_antigen, calibration specs      │
│  └─ gold_standard_path, IEDB filter                     │
│                                                           │
│  ↓                                                        │
│                                                           │
│  Phase 1-2: Data & Prediction (unchanged)               │
│  └─ UniProt sequences, IEDB epitopes, NetMHCIIpan      │
│                                                           │
│  ↓                                                        │
│                                                           │
│  Phase 3: ENHANCED SCORING                              │
│  ├─ Criterion 1: MHC Zone (0.20 weight)                 │
│  ├─ Criterion 2: HLA Promiscuity (0.20)                 │
│  ├─ Criterion 3: Disease Proximity (0.30)               │
│  ├─ Criterion 4: Treg TCR Contact (0.15) [NEW]          │
│  ├─ Criterion 5: IFN-gamma Penalty (0.07)               │
│  ├─ Criterion 6: Solubility/GRAVY (0.08)                │
│  ├─ Criterion 7: JMX Self-Sim (0.00)                    │
│  ├─ Bonus: Cathepsin S processing likelihood            │
│  ├─ Bonus: B-cell risk flagging (protein-relative)      │
│  └─ Output: Ranked peptide scores CSV                   │
│                                                           │
│  ↓                                                        │
│                                                           │
│  Calibration Check                                       │
│  └─ Validate against disease-specific targets (JSON)    │
│                                                           │
│  ↓                                                        │
│                                                           │
│  Phase 4: CONSTRUCT ASSEMBLY                            │
│  ├─ Positional diversity selection (top N peptides)     │
│  ├─ Multi-epitope assembly (GPGPG or AAY linker)        │
│  ├─ Construct-level scoring (bonuses + penalties)       │
│  ├─ Population coverage (Bui formula, 4 populations)    │
│  ├─ Junction epitope checking                           │
│  ├─ Codon optimization (GC-balanced, seed=42)           │
│  ├─ mRNA generation (5'UTR + CDS + 3'UTR + polyA)       │
│  ├─ Experimental tier assignment (T1/T2/T3)             │
│  └─ Output: mRNA constructs + peptide details CSV       │
│                                                           │
└─────────────────────────────────────────────────────────┘
```

---

## MAJOR CHANGES FROM v1 TO v2+

### Conceptual
1. **From ITP-specific to disease-agnostic:** Profiles enable any autoimmune disease
2. **From IL-10 model to TCR biology:** Evidence-based replacement after data leakage investigation
3. **From binary JMX to frequency-based scoring:** TCR-contact motif approach has better variance

### Technical
1. **New modules:** esm_embeddings.py, processing.py, disease_profile.py, population_coverage.py
2. **Enhanced existing:** scorer.py (new criteria, weights), construct_builder.py (JMX proxy, B-cell filter)
3. **Data files:** Two disease profiles, two gold standards
4. **Documentation:** Multiple analysis reports + upgrade log

### Performance
1. **ITP:** Calibration passes (Peptide_82 rank 2, Peptide_2 rank 53)
2. **MS:** Calibration passes (MBP131-145 rank 1, MBP30-44 rank 2)
3. **Stability:** Top 10 rankings largely preserved (8-9/10 overlap)
4. **Discrimination:** Treg TCR scores range 0.27-1.00 (good variance)

---

END OF COMPREHENSIVE FEATURE SUMMARY
