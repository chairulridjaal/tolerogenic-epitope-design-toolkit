# Tolerogenic Scoring Criteria

This document defines the scoring criteria used in Phase 3 of the 
Tolerogenic Epitope Design Toolkit. Each criterion is derived from 
published experimental literature. Confidence levels reflect the 
strength and directness of the evidence as it applies to ITP 
specifically.

Every criterion is implemented as a normalized score in [0, 1].
The final tolerogenic score is a weighted sum across all criteria.
Weights are documented in `src/scoring/scorer.py` and are explicitly
marked as tunable — they reflect our current best interpretation of 
the literature, not experimentally validated parameters.

The pipeline is calibrated against seven experimentally validated 
GPIIIa immunodominant peptides from ITP patients (see 
`data/processed/itp_gold_standard.json`), two of which 
(Peptide 2 and Peptide 82) have been confirmed as tolerogenic 
candidates in humanized HLA-DR15 mice. A correctly weighted scorer 
must rank these peptides in the top 20% of all ITGB3 candidates.

---

## Criterion 1 — MHC Binding in the Tolerogenic Zone

**Confidence: High**

**Biological rationale:**
T cell activation outcome is determined not just by whether a peptide
binds MHC, but by the *strength* of the resulting TCR signal. A 
quantitative "Goldilocks" model is supported across multiple studies:
moderate-affinity peptide-MHC complexes favor regulatory T cell 
induction, while very high affinity drives either aggressive effector 
activation, activation-induced cell death (AICD), or — in existing 
Tregs — FoxP3 destabilization and conversion to a pro-inflammatory 
phenotype. Very low affinity peptides fail to be presented reliably.

Specific IC50 thresholds from the Tregitope literature define the 
zones: tolerogenic candidates exhibit IC50 values in the 
100–10,000 nM range (moderate-to-high affinity), while 
supraphysiological binders (IC50 < 100 nM) correlate with 
pro-inflammatory outcomes. In percentile rank terms (as returned by NetMHCIIpan): the 
tolerogenic zone corresponds to ranks of 2–10%. Ranks below 2% 
represent extreme binders (penalized). Ranks above 10% represent 
weak or non-binders (also penalized, as reliable presentation is 
required).

**Implementation:**
For each peptide-allele prediction, assign a score based on 
percentile rank:
- Rank 2–10%: score 1.0 (optimal)
- Rank 10–20%: score 0.5 (acceptable weak binder)
- Rank <2%: score 0.2 (penalized — extreme binder)
- Rank >20%: score 0.0 (insufficient binding)

Aggregate across all alleles by taking the mean score.

**Primary citations:**
- Evavold BD & Allen PM (1991). Induction of T-cell anergy by altered T-cell-receptor ligand on live antigen-presenting cells. *Nature* **356**:604–607. [doi:10.1038/356604a0](https://doi.org/10.1038/356604a0)
- De Groot AS *et al.* (2008). Activation of natural regulatory T cells by IgG Fc–derived peptide “Tregitopes.” *Blood* **112**(7):3303–3311. [doi:10.1182/blood-2008-02-138073](https://doi.org/10.1182/blood-2008-02-138073)

---

## Criterion 2 — HLA Promiscuity / Population Coverage

**Confidence: High**

**Biological rationale:**
Tolerogenic epitopes that bind promiscuously across multiple HLA 
alleles have two advantages. First, they provide therapeutic coverage 
across a genetically diverse patient population. Second, promiscuous 
binding is itself a property of known Tregitopes — naturally 
occurring regulatory T cell epitopes found in conserved human 
proteins. Parasitic Tregitopes exploiting this same mechanism bind 
HLA-DR across >80% of the human population, hijacking the host 
regulatory repertoire.

ITP-specific HLA associations (DRB1*04:05 for anti-GPIIb/IIIa, 
DRB1*08:03 for anti-GPIb-IX) confirm that disease susceptibility 
is HLA-linked, reinforcing that a broadly promiscuous tolerogenic 
epitope is preferable to one restricted to a single allele.

**Implementation:**
Count the number of alleles in HLA_PANEL for which the peptide 
achieves percentile rank ≤ 10%. Normalize by panel size (12 alleles).

Score = (alleles with rank ≤ 10%) / len(HLA_PANEL)

A peptide binding 8 of 12 alleles scores 0.67. A peptide binding 
all 12 scores 1.0.

**Primary citations:**
- Serra P & Santamaria P (2019). Antigen-specific therapeutic approaches for autoimmunity. *Nature Biotechnology* **37**:238–251. [doi:10.1038/s41587-019-0015-4](https://doi.org/10.1038/s41587-019-0015-4)

---

## Criterion 3 — ITP-Validated Peptide Proximity

**Confidence: High**

**Biological rationale:**
Seven 15-mer peptides from the GPIIIa (ITGB3) subunit have been 
experimentally mapped as immunodominant T cell epitopes in ITP 
patients, eliciting recall responses in over 80% of patients tested. 
Two of these — Peptide 2 (aa6–20: GDCNCTKDDSVMCIG) and Peptide 82 
(aa711–725: NPIYKSAVTTVVNP) — have been validated as tolerogenic 
candidates in humanized HLA-DR15 mice, successfully inducing 
antigen-specific Tregs and suppressing anti-GPIIb/IIIa autoantibody 
production.

Note on bystander suppression: tolerogenic peptides are not required 
to share exact specificity with disease-driving effector T cells. 
Tregs induced by any peptide from the target tissue will migrate to 
the site of inflammation and suppress locally via IL-10 and TGF-β 
secretion (bystander suppression). However, peptides with direct 
experimental evidence of tolerogenicity in ITP receive a strong 
priority bonus, as they represent the highest-confidence candidates.

**Implementation:**
Compute substring overlap between the candidate peptide and each of 
the 7 gold-standard sequences. Scoring tiers:

- Exact match or full containment with Peptide 2 or 82: score 1.0
- Exact match or full containment with any of the 7 peptides: 
  score 0.8
- Partial overlap (≥ 9 residues shared) with any of the 7: 
  score 0.4
- No overlap: score 0.0

**Gold standard sequences** 
(source: `data/processed/itp_gold_standard.json`):

| ID         | Position   | Sequence (core 15-mer)     | Notes                                      |
|------------|------------|----------------------------|--------------------------------------------|
| Peptide 2  | aa6–20     | GDCNCTKDDSVMCIG ★          | Confirmed tolerogenic in DR15 mice         |
| Peptide 44 | aa331–345  | WNIQNPWSIKRKRK             | Used with RK wrapper for solubility        |
| Peptide 48 | aa361–375  | SYWNFGNNVTLHKKS            | —                                          |
| Peptide 56 | aa421–435  | TEEKIQLIVQANIDQ            | —                                          |
| Peptide 70 | aa591–605  | GCVYIEDEVHRCYGR            | —                                          |
| Peptide 77 | aa661–675  | CRKSLVSSFAECRTF            | —                                          |
| Peptide 82 | aa711–725  | NPIYKSAVTTVVNP ★           | Confirmed tolerogenic (often RK-flanked)   |

★ = validated tolerogenic candidates in Hall et al. (2019).

**Primary citations:**
- Sukati H *et al.* (2007). Mapping helper T-cell epitopes on platelet membrane glycoprotein IIIa in chronic autoimmune thrombocytopenic purpura. *Blood* **109**(10):4528–4538. [doi:10.1182/blood-2006-09-044388](https://doi.org/10.1182/blood-2006-09-044388)
- Hall LS *et al.* (2019). Combination peptide immunotherapy suppresses antibody and helper T-cell responses to the major human platelet autoantigen glycoprotein IIb/IIIa in HLA-transgenic mice. *Haematologica* **104**(5):1079–1087. [doi:10.3324/haematol.2017.179424](https://doi.org/10.3324/haematol.2017.179424)

---

## Criterion 4 — IL-10 Induction Potential

**Confidence: Medium-High**

**Biological rationale:**
IL-10 is the primary immunosuppressive cytokine required for 
peripheral tolerance. It is produced by Tr1 regulatory cells and 
is the key mediator through which peptide-induced tolerance operates 
in vivo. Peptides that favor IL-10 production over IFN-gamma 
production are more compatible with a tolerogenic outcome. The 
cytokine balance is partially determined by peptide sequence — 
specific dipeptide compositions correlate with IL-10 vs. IFN-gamma 
induction probability.

**Implementation:**
A local Random Forest model is trained on the exact 73 features from
Table S1 of Nagpal et al. (2017): 16 amino acid composition + 57
dipeptide composition features. Training data: Table S4 (394 positive,
848 negative peptides). The model outputs P(IL-10 inducer) in [0, 1].

The local RF model reproduces the original 73-feature set from
Nagpal et al. (2017) with AUC 0.91 on the independent validation set
(Table S5); accuracy is 0.67 (slightly lower than the paper's internal
CV of ~0.81, as expected for an independent holdout).

Train once: ``python -m src.scoring.train_il10_model``

**Primary citations:**
- Nagpal G *et al.* (2017). Computer-aided designing of immunosuppressive peptides based on IL-10 inducing potential. *Scientific Reports* **7**:42851. [doi:10.1038/srep42851](https://doi.org/10.1038/srep42851)

---

## Criterion 5 — IFN-gamma Penalty

**Confidence: Medium-High**

**Biological rationale:**
IFN-gamma is the defining cytokine of Th1 effector responses — 
the pro-inflammatory arm of CD4+ T cell immunity. Peptides that 
strongly induce IFN-gamma production are driving exactly the 
opposite of tolerance. Computational tools can predict IFN-gamma 
induction probability from sequence features.

**Implementation:**
Uses IFNepitope2 (Dhall et al., 2024), the official successor to the
2013 tool. Key improvements: host-specific models (human/mouse),
hybrid ML + BLAST method, AUROC 0.90 (human). Runs locally via the
``ifnepitope2`` Python package (ExtraTrees model, DPC features).

Score = 1.0 − P(IFN-γ inducer)

High IFN-gamma prediction → low tolerogenic score.

**Primary citations:**
- Dhanda SK, Vir P & Raghava GPS (2013). Designing of interferon-gamma inducing MHC class-II binders. *Biology Direct* **8**:30. [doi:10.1186/1745-6150-8-30](https://doi.org/10.1186/1745-6150-8-30)
- Dhall A *et al.* (2024). IFNepitope2: improved prediction of interferon-gamma inducing peptides. *Scientific Reports*.

---

## Criterion 6 — Solubility (GRAVY Score)

**Confidence: Medium**

**Biological rationale:**
Peptide solubility is a prerequisite for tolerogenic efficacy. 
Insoluble or highly hydrophobic peptides aggregate at the injection 
site and fail to reach tolerogenic dendritic cells in lymphoid 
organs. Studies of tolerogenic peptide design specifically show that 
adding charged lysine/arginine wrappers to improve solubility is 
necessary for successful IL-10 induction in vivo — and conversely, 
hydrophobic LF motifs that reduce solubility abolish the 
tolerogenic effect.

The Grand Average of Hydropathy (GRAVY) score is computable 
directly from sequence with no API call. Negative values indicate 
hydrophilic, soluble peptides. Positive values indicate hydrophobic, 
aggregation-prone peptides.

**Implementation:**
Compute GRAVY score using the Kyte-Doolittle scale:

GRAVY = sum(hydropathy values for each amino acid) / peptide length

Normalize: score = 1.0 if GRAVY < −0.5, linear decay to 0.0 
at GRAVY = +1.0.

Computable with no external dependency using a lookup table.

**Primary citations:**
- Kyte J & Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* **157**:105–132. [doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0)

---

## Criterion 7 — Human Proteome Cross-Conservation

**Confidence: Medium**

**Biological rationale:**
Peptides with high sequence similarity to other human self-proteins 
are more likely to engage natural regulatory T cells (nTregs), 
because those Tregs were selected in the thymus to recognize 
self-peptides. The JanusMatrix (JMX) tool quantifies this 
cross-conservation. High JMX scores predict Treg activation; low 
JMX scores predict effector T cell activation. This is mechanistically 
related to why Tregitopes (which are highly conserved across human 
IgG molecules) are so potent at inducing tolerance.

**Implementation:**
Submit to JanusMatrix web server 
(https://janusmatrix.essentialfacts.com) and retrieve JMX score. 
Normalize to [0, 1] based on the score distribution across all 
candidate peptides.

Note: JanusMatrix has no public API. This criterion is flagged as 
manually computable for exploratory use, with a fallback of 
sequence identity to the human UniProt reviewed proteome using 
BLAST or local alignment.

**Primary citations:**
- Moise L *et al.* (2013). The two-faced T cell epitope: examining the host-microbe interface with JanusMatrix. *Human Vaccines & Immunotherapeutics* **9**(7):1577–1586. [doi:10.4161/hv.24615](https://doi.org/10.4161/hv.24615)

---

## Composite Scoring

The final tolerogenic score is a weighted sum:
```
score = (
    w1 * mhc_zone_score +
    w2 * hla_promiscuity_score +
    w3 * itp_proximity_score +
    w4 * il10_score +
    w5 * ifng_penalty_score +
    w6 * solubility_score +
    w7 * jmx_score
)
```

**Default weights** (tunable, calibrated against gold standard):

| Criterion | Weight | Rationale |
|---|---|---|
| MHC zone | 0.20 | Necessary condition for presentation |
| HLA promiscuity | 0.20 | Population coverage — practical priority |
| ITP proximity | 0.25 | Highest direct experimental evidence |
| IL-10 induction | 0.15 | Key mechanistic output |
| IFN-gamma penalty | 0.10 | Safety filter |
| Solubility | 0.05 | Practical requirement |
| JMX cross-conservation | 0.05 | Mechanistic support, lower confidence |

Total: 1.00

Weights are stored as constants in `src/scoring/scorer.py` and 
can be adjusted by the user. Any weight adjustment should be 
documented with biological justification.

---

## Calibration Requirement

Before any results are reported, the scoring pipeline must pass 
this calibration check:

Run the full scorer on all 774 ITGB3 15-mer candidates. 
Peptide 2 (GDCNCTKDDSVMCIG) and Peptide 82 (NPIYKSAVTTVVNP) 
must both rank in the top 20% by composite score. If they do not, 
weights must be adjusted and the adjustment documented here with 
justification.

---

## Limitations and Honest Uncertainties

1. All scoring criteria are derived from preclinical models, 
   primarily mouse EAE (MS model) and NOD mouse (T1D model). 
   Direct experimental validation in ITP models is limited to 
   Criterion 3.

2. The tolerogenic zone for Criterion 1 (percentile rank 2–10%) 
   is an approximation. The exact optimal range likely varies by 
   HLA allele and patient population.

3. The local IL-10 RF model and IFNepitope2 were trained on general
   T cell epitope datasets, not ITP-specific data. Their predictions
   for ITP antigens are extrapolations. The IL-10 model achieves
   AUC 0.91 on the independent validation set (S5) but accuracy of
   0.67 (vs. the paper's internal CV ~0.81) — probability rankings
   are reliable, binary thresholds less so.

4. Criterion 7 (JMX) has no public API and limited validation 
   outside the EpiVax research group's own publications.

5. This pipeline produces computational hypotheses only. No 
   criterion has been validated experimentally in the context 
   of this specific pipeline or for ITP specifically. All outputs 
   require expert immunological review before any downstream use.