# A disease-agnostic computational pipeline for tolerogenic epitope design: case studies in immune thrombocytopenia and multiple sclerosis

## Abstract

Antigen-specific tolerogenic therapy offers the prospect of selectively restoring immune self-tolerance in autoimmune disease without broad immunosuppression, but no integrated computational pipeline exists for rational design of tolerogenic epitope candidates. We present an open-source pipeline that scores peptide candidates from autoimmune disease target antigens against seven literature-grounded criteria — MHC-II binding zone, HLA allele promiscuity, disease-validated peptide proximity, Treg TCR-contact self-similarity, IFN-gamma penalty, solubility, and human proteome cross-conservation — and assembles top candidates into codon-optimized multi-epitope mRNA vaccine constructs.

Applied to immune thrombocytopenia (ITP) targeting GPIIIa/ITGB3 (257 candidates from 774 15-mers) and multiple sclerosis (MS) targeting myelin basic protein (113 candidates from 290 15-mers), the pipeline independently recovers known immunodominant protein regions without gold standard labels. In MS, all five of the top independently-ranked 15-mers are overlapping windows from the MBP83-99 immunodominant region targeted by the ATX-MS-1467 clinical product (Kappos et al. 2000, Streeter et al. 2015). In ITP, the C-terminal cytoplasmic tail of ITGB3 (aa749-765), adjacent to the validated tolerogenic Peptide 82 region, occupies independent ranks 1-2 and 4, with the highest Treg TCR-contact self-similarity (0.99) and HLA promiscuity (0.67-0.75) in the dataset.

During development, we discovered that the IL-10 induction model originally used as Criterion 4 (Nagpal et al. 2017, trained on 394 peptides) achieves AUC 0.50 on truly non-overlapping validation data — random chance — due to 36.6% sequence overlap between the published training and validation sets. Neither the original composition features nor ESM-2 protein language model embeddings (320-dimensional) rescue generalization. We replaced this criterion with Treg TCR-contact self-similarity, a non-ML score measuring the frequency of TCR-facing residues (MHC-II binding core positions P2, P3, P5, P7, P8) in the human proteome, grounded in the Tregitope framework (De Groot et al. 2008, Moise et al. 2013). This criterion identified and corrected a specific false positive: ENPVVHFFKNIVTPR (MBP aa217-231), previously the top-ranked novel MS candidate due to inflated IL-10 scores, was correctly demoted to rank 55 based on its foreign-like TCR contacts (score 0.56 vs median 0.77).

The pipeline adds a new disease in two JSON configuration files with zero code changes and reproduces from empty data directories in approximately nine minutes. Code and data are available at https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit.

---

## 1. Introduction

Autoimmune diseases collectively affect 5-8% of the global population and are characterised by breakdown of immune tolerance to self-antigens. Current treatments suppress immunity broadly: corticosteroids and splenectomy for ITP, interferon-beta and natalizumab for MS. Antigen-specific tolerogenic therapy — administering self-antigen peptides in a context that promotes regulatory T cell induction rather than effector activation — offers disease-specific immune re-education without compromising systemic immunity (Serra and Santamaria 2019, Larche and Wraith 2005).

Clinical proof-of-concept exists for both diseases studied here. In ITP, Hall et al. (2019) demonstrated that GPIIIa peptides aa32-46 and aa737-751 (UniProt P05106 precursor numbering; mature protein aa6-20 and aa711-725) suppressed anti-GPIIb/IIIa antibody responses and induced FoxP3+ Tregs in HLA-DR15 transgenic mice. The tolerogenic mechanism was IL-10 independent — suppression was mediated by classical CD4+CD25+FoxP3+ regulatory T cells (Hall et al. 2019, Figure 5). In MS, ATX-MS-1467, a mixture of four MBP peptides including the immunodominant MBP83-99 region, completed Phase I and Phase II clinical trials demonstrating reduction in MBP-specific T cell responses and Th2 cytokine skewing (Streeter et al. 2015, Kappos et al. 2000).

Despite this clinical progress, no integrated computational pipeline exists for systematic design of tolerogenic epitope candidates. Existing tools address components — NetMHCIIpan for MHC-II binding prediction (Reynisson et al. 2020), IFNepitope2 for cytokine prediction (Dhall et al. 2024), JanusMatrix for cross-conservation analysis (Moise et al. 2013) — but are not assembled into a disease-agnostic workflow incorporating tolerogenic scoring criteria grounded in the relevant experimental literature.

We present the Tolerogenic Epitope Design Toolkit, an open-source Python pipeline that integrates MHC binding prediction, seven tolerogenic scoring criteria, and multi-epitope mRNA construct assembly. We apply it to ITP and MS as case studies, characterise the relationship between computational predictions and experimentally validated tolerogenic peptides, report the discovery and correction of a data leakage problem in a widely-used cytokine prediction model, and identify novel candidate regions for experimental follow-up.

---

## 2. Methods

### 2.1 Disease profile system

Each disease is specified by a JSON configuration file containing target antigens (UniProt accessions), IEDB disease filter, gold standard peptide path, and calibration specification. Adding a new disease requires two files with no code modifications.

### 2.2 Antigen sequences and epitope data

Protein sequences were retrieved from UniProt via the REST API. For ITP: ITGA2B (P08514), ITGB3 (P05106), GP1BA (P07359), GP1BB (P13224), GP9 (P14770), GP5 (P40197). For MS: MBP (P02686), PLP1 (P60201), MOG (Q16653). Known epitopes were retrieved from the IEDB query API filtered by disease annotation and human host organism. Gold standard peptide positions were verified against UniProt sequences at pipeline startup; a position mismatch raises an immediate error with guidance on signal peptide offset correction.

### 2.3 MHC-II binding prediction

Overlapping 15-mer peptides were generated from the primary antigen (ITGB3 for ITP, 774 peptides from 788 residues; MBP for MS, 290 peptides from 304 residues). Binding predictions were obtained via the IEDB tools API using NetMHCIIpan 4.0 (eluted ligand predictor) across 12 HLA class II alleles: HLA-DRB1*01:01, *03:01, *04:01, *04:05, *07:01, *09:01, *11:01, *13:01, *15:01, HLA-DQA1*01:01/DQB1*05:01, HLA-DQA1*01:02/DQB1*06:02, and HLA-DPA1*01:03/DPB1*04:01. Peptides with percentile rank ≤10% for at least one allele (590 for ITP, 260 for MS) were retained for scoring; gold standard peptides were additionally included regardless of binding rank.

### 2.4 Tolerogenic scoring criteria

Seven criteria, each normalised to [0, 1], are combined in a weighted composite score (Table 1).

**Table 1.** Tolerogenic scoring criteria and weights.

| # | Criterion | Weight | Implementation |
|---|-----------|--------|----------------|
| 1 | MHC binding zone | 0.20 | Percentile rank 2-10% scores 1.0 (tolerogenic zone); <2% scores 0.2 (penalised); 10-20% scores 0.5; >20% scores 0.1 |
| 2 | HLA promiscuity | 0.20 | Fraction of 12-allele panel with rank ≤10% |
| 3 | Disease-validated proximity | 0.30 | Substring overlap with experimentally validated peptides |
| 4 | Treg TCR-contact self-similarity | 0.08 | Frequency of TCR-facing residues in human proteome (see 2.5) |
| 5 | IFN-gamma penalty | 0.07 | Inverted IFNepitope2 probability (Dhall et al. 2024) |
| 6 | Solubility | 0.08 | GRAVY hydropathy score (Kyte and Doolittle 1982) |
| 7 | Human proteome cross-conservation | 0.07 | Fraction of 9-mer windows in Swiss-Prot human proteome |

The >20% MHC binding zone receives 0.1 rather than 0.0 to accommodate ITP autoepitopes, which are characteristically low-affinity MHC binders. Sukati et al. (2007, Table 6) demonstrated that most immunodominant GPIIIa peptides have no predicted high-affinity HLA-DR binding — a finding we confirm with NetMHCIIpan (Peptide 2: percentile rank 98% for HLA-DRB1*15:01). Their immunodominance arises from thymic escape: poor MHC presentation prevents thymic deletion of autoreactive T cells (Liu et al. 1995).

### 2.5 Treg TCR-contact self-similarity (Criterion 4)

In MHC-II presentation, the 9-mer binding core occupies the groove with anchor residues (P1, P4, P6, P9) facing the MHC and TCR-contact residues (P2, P3, P5, P7, P8) facing the T cell receptor. The TCR contacts determine whether the responding T cell is regulatory or effector. Tregitopes — naturally occurring regulatory T cell epitopes — are characterised by TCR-facing residues that are highly conserved in the human proteome (De Groot et al. 2008, Moise et al. 2013).

We computed a TCR-contact frequency index from the Swiss-Prot reviewed human proteome (20,431 proteins, 10,406,322 unique 9-mers): for each 9-mer, the 5 TCR-facing residues were extracted and their frequency counted (2,356,207 unique motifs, covering 73.6% of the 3,200,000 possible 5-mer space). For each candidate peptide, the 9-mer binding core(s) from NetMHCIIpan were retrieved, TCR motifs extracted, and their frequency converted to a percentile score. Higher frequency indicates a more self-like TCR surface associated with natural Treg engagement.

This criterion replaced an IL-10 induction model (see Section 3.5). It requires no ML model, no training data, and carries no generalization risk.

### 2.6 IFN-gamma penalty (Criterion 5)

IFN-gamma induction probability was predicted using IFNepitope2 (Dhall et al. 2024), run locally via the ifnepitope2 Python package with a scikit-learn compatibility shim. Dipeptide composition features were computed and scored by the bundled ExtraTrees classifier. The tolerogenic score is 1.0 minus the predicted probability.

An overlap audit confirmed that 33 of 370 candidate peptides (8.9%) appear in the IEDB T-cell dataset from which IFNepitope2 was trained. This overlap is expected (ITGB3 is a studied autoantigen) and does not represent train-test leakage: IFNepitope2 score distributions are indistinguishable between overlapping (mean 0.553 ± 0.130) and non-overlapping (mean 0.568 ± 0.121) peptides.

### 2.7 B-cell epitope safety filter

Linear B-cell epitope risk was assessed using the Parker hydrophilicity scale (Parker et al. 1986) with a sliding 7-residue window and protein-relative thresholding: a peptide is flagged only if its maximum window score is in the top 20th percentile for its source protein. This avoids false positives for intrinsically disordered proteins like MBP where baseline hydrophilicity is uniformly high.

### 2.8 Construct assembly

Top-ranked peptides were selected using a positional diversity filter (minimum 10-residue separation) to prevent overlapping 15-mers from dominating the construct. Peptides were joined with GPGPG flexible linkers. Codon optimisation used weighted random sampling from human Kazusa codon usage frequencies with an 18-codon sliding GC window constraint (maximum 62%), targeting 50-60% GC. Constructs include Kozak 5'UTR, human beta-globin 3'UTR, 120-nucleotide poly-A tail, and annotations for N1-methylpseudouridine modification and dsRNA-depleted LNP formulation.

### 2.9 Calibration

Disease-specific calibration checks are defined in the disease profile. For ITP: Peptide 82 (ALLIWKLLITIHDRK, genuine HLA-DR binder per Sukati 2007 Table 6) must rank in the top 20% with mhc_zone > 0.5; Peptide 2 (TTRGVSSCQQCLAVS, thymic escape autoepitope) must achieve itp_proximity ≥ 1.0 and jmx ≥ 1.0 regardless of composite rank. For MS: MBP83-99 and MBP30-44 must achieve proximity ≥ 0.4 and jmx ≥ 0.5.

### 2.10 Independent score analysis

To separate findings dependent on prior gold standard knowledge from genuinely predictive results, we computed an independent score excluding Criterion 3 (disease-validated proximity, weight 0.30). The remaining six criteria were renormalised to sum to 1.0. All results described as "independent" use this proximity-excluded score.

---

## 3. Results

### 3.1 Calibration

All calibration checks passed for both diseases (Table 2). In ITP, Peptide 77 (DDCVVRFQYYEDSSG) ranked first overall — the only gold standard peptide with both high MHC affinity (zone 1.00) and strong Treg TCR-contact self-similarity (0.79). Peptide 82 ranked second (composite 0.672). Peptide 2 ranked 53rd; its low MHC binding (rank 98% for DRB1*15:01) is the expected biology for a thymic escape autoepitope. In MS, MBP131-145 and MBP30-44 (both ATX-MS-1467 components) ranked first and second, with Treg TCR-contact scores of 0.95 and 0.93 respectively.

**Table 2.** Gold standard peptide rankings.

| Disease | Peptide | Rank | Composite | Treg TCR | MHC zone | Validated |
|---------|---------|------|-----------|----------|----------|-----------|
| ITP | Peptide_77 | 1/257 | 0.733 | 0.79 | 1.00 | No |
| ITP | Peptide_82 | 2/257 | 0.672 | 0.60 | 0.73 | Hall 2019 |
| ITP | Peptide_44 | 3/257 | 0.632 | 0.70 | 1.00 | No |
| ITP | Peptide_2 | 53/257 | 0.507 | 0.50 | 0.10 | Hall 2019 |
| MS | MBP131-145 | 1/113 | 0.788 | 0.95 | 1.00 | ATX-MS-1467 Phase I |
| MS | MBP30-44 | 2/113 | 0.781 | 0.93 | 1.00 | ATX-MS-1467 Phase I |
| MS | MBP140-154 | 20/113 | 0.552 | 0.50 | 0.13 | ATX-MS-1467 Phase I |
| MS | MBP83-99 | 21/113 | 0.543 | 0.50 | 0.10 | Phase I/II |

The ITP gold standard splits into two mechanistic classes: three conventional MHC binders (ranks 1-3) and four thymic escape peptides (ranks 53-193). This split is consistent with Sukati et al. (2007) Table 6 and validates the pipeline's sensitivity to mechanistically distinct autoepitope classes.

### 3.2 Independent regional recovery — MS

When disease proximity is excluded, all five top-ranked MBP peptides are overlapping 15-mers from the MBP83-99 immunodominant window (Table 3). This region is the primary target of ATX-MS-1467, which has shown clinical activity (Kappos et al. 2000, Streeter et al. 2015). The MBP131-145 region (another ATX-MS-1467 component) appears at independent ranks 9, 12, and 14. The pipeline identified both regions without any gold standard information, driven by MHC binding zone scores in the 2-10% tolerogenic range, high Treg TCR-contact self-similarity (0.83-0.96 for MBP83-99 overlapping 15-mers), and uniform proteome cross-conservation (JMX 1.00 for all MBP peptides).

No exact gold standard peptide appears in the top 20% by independent score. The pipeline identifies the correct protein regions, not the exact clinical peptides — consistent with the biology of MHC-II antigen processing, where the precise peptide within a region is shaped by protease activity and endosomal chemistry not modelled here.

**Table 3.** Top 5 MBP peptides by independent score (proximity excluded).

| Rank | Peptide | Position | MHC | HLA | Treg TCR | Ind. score | Region |
|------|---------|----------|-----|-----|----------|------------|--------|
| 1 | RPHLIRLFSRDAPGR | aa88-102 | 1.00 | 0.42 | 0.83 | 0.788 | MBP83-99 |
| 2 | LIRLFSRDAPGREDN | aa91-105 | 1.00 | 0.33 | 0.94 | 0.777 | MBP83-99 |
| 3 | IRLFSRDAPGREDNT | aa92-106 | 1.00 | 0.25 | 0.96 | 0.762 | MBP83-99 |
| 4 | SRPHLIRLFSRDAPG | aa87-101 | 1.00 | 0.33 | 0.83 | 0.757 | MBP83-99 |
| 5 | HLIRLFSRDAPGRED | aa90-104 | 1.00 | 0.25 | 0.96 | 0.755 | MBP83-99 |

### 3.3 Independent regional recovery — ITP

For ITP, the independent analysis reveals a qualitatively different pattern consistent with the known biology of ITP autoepitopes. Six of seven gold standard peptides rank in the bottom 20% of independent scores — expected, because ITP immunodominant peptides are low-affinity MHC binders by design (Sukati et al. 2007).

However, the protein regions adjacent to validated peptides emerge independently. The C-terminal cytoplasmic tail of ITGB3 (aa749-765), 14 residues downstream from the validated Peptide 82 region (aa737-751), occupies independent ranks 1, 2, and 4, with the highest Treg TCR-contact scores in the dataset (0.99) and the strongest HLA promiscuity (0.67-0.75, binding 8-9 of 12 alleles). The region adjacent to Peptide 77 (aa661-676) occupies ranks 3 and 5.

**Table 4.** Top 5 ITGB3 peptides by independent score (proximity excluded).

| Rank | Peptide | Position | MHC | HLA | Treg TCR | Ind. score | Region |
|------|---------|----------|-----|-----|----------|------------|--------|
| 1 | KEFAKFEEERARAKW | aa751-765 | 0.90 | 0.67 | 0.99 | 0.820 | Near Pep82 |
| 2 | DRKEFAKFEEERARA | aa749-763 | 0.73 | 0.75 | 0.99 | 0.810 | Near Pep82 |
| 3 | RDEIESVKELKDTGK | aa662-676 | 0.89 | 0.58 | 0.99 | 0.802 | Near Pep77 |
| 4 | RKEFAKFEEERARAK | aa750-764 | 0.73 | 0.75 | 0.99 | 0.798 | Near Pep82 |
| 5 | CRDEIESVKELKDTG | aa661-675 | 0.87 | 0.50 | 0.99 | 0.782 | Near Pep77 |

### 3.4 Construct outputs

Both disease constructs incorporate 10 positionally diverse epitopes with GC content within the 50-60% therapeutic target (Table 5).

**Table 5.** Multi-epitope mRNA construct summary.

| Construct | Epitopes | AA length | mRNA length | GC% | Score |
|-----------|----------|-----------|-------------|-----|-------|
| ITP-GPGPG-10ep | 10 | 195 | 768 nt | 53.7% | 0.887 |
| MS-GPGPG-10ep | 10 | 195 | 768 nt | 57.3% | 0.885 |

Both constructs include manufacturing annotations for N1-methylpseudouridine modification, CleanCap AG, and dsRNA-depleted LNP formulation for spleen/liver targeting.

### 3.5 IL-10 induction model audit

The pipeline originally used Criterion 4 as an IL-10 induction score based on a Random Forest classifier trained on Nagpal et al. (2017) data (S4: 394 positive, 848 negative peptides; 73 dipeptide composition features). Initial evaluation on the published validation set (S5: 1,082 peptides) yielded AUC 0.91.

However, we discovered that 396 sequences (36.6% of S5) are shared with S4, and all 396 are positive examples. When these memorised sequences are excluded, the AUC on the 686 truly non-overlapping S5 sequences (65 positive, 621 negative) is 0.50 — random chance. The model memorises training sequences (99.5% correct, mean probability 0.81) but has zero generalization to unseen peptides.

We tested ESM-2 protein language model embeddings (esm2_t6_8M_UR50D, 320-dimensional mean-pooled representations) as an alternative feature representation. Three classifiers trained on ESM-2 embeddings achieved non-overlapping S5 AUC of 0.51 (logistic regression), 0.46 (random forest), and 0.51 (MLP). The failure to generalize is therefore not a feature representation problem but a training data limitation: 394 positive examples may be insufficient for any model to learn sequence-general patterns of IL-10 induction, or IL-10 induction may not be purely sequence-determined.

This criterion was replaced with Treg TCR-contact self-similarity (Section 2.5), which requires no training data and carries no generalization risk. The replacement identified and corrected a specific false positive: ENPVVHFFKNIVTPR (MBP aa217-231), which had the highest IL-10 score in the MS dataset (0.88), was correctly demoted from an apparent top candidate to rank 55 based on its foreign-like TCR contacts (score 0.56 vs median 0.77).

### 3.6 Weight sensitivity analysis

We evaluated eight alternative weight configurations. MS results were robust: MBP131-145 and MBP30-44 ranked 1-2 under all configurations except those removing proximity entirely. ITP results were weight-sensitive: Peptide 82 ranged from rank 2 (current weights) to rank 233 (no proximity), consistent with the thymic escape biology — ITP autoepitopes score poorly on independent criteria by design.

The regional recovery findings are robust across configurations: the MBP83-99 cluster and ITGB3 C-terminal tail dominate independent rankings regardless of how the six non-proximity criteria are weighted.

---

## 4. Discussion

The pipeline demonstrates two distinct modes of validated target recovery across mechanistically different autoimmune diseases. In MS, the six independent scoring criteria (MHC binding, HLA promiscuity, Treg TCR-contact, IFN-gamma, solubility, proteome conservation) directly identify the ATX-MS-1467 target regions because those peptides are genuine moderate-affinity MHC binders with self-like TCR contacts. In ITP, the same criteria identify neighboring rather than overlapping regions because the validated tolerogenic peptides are low-affinity binders whose immunodominance arises from thymic escape — a fundamentally different mechanism that MHC-based scoring cannot capture by design.

This disease-specific contrast validates the pipeline's sensitivity to the underlying biology rather than a one-size-fits-all ranking. The ITP result in particular illustrates a known limitation of MHC-based epitope prediction for autoimmune disease: as Sukati et al. (2007) established, immunodominant ITP epitopes are selected by failed thymic deletion rather than efficient MHC presentation. Computational pipelines that discard low-affinity peptides will systematically miss this class of autoepitopes.

The IL-10 model audit represents a cautionary finding for the field. The Nagpal et al. (2017) model achieves AUC 0.91 on its published validation set, a number that has been cited as evidence of strong predictive performance. We show that this metric is entirely driven by 396 memorised training sequences (36.6% overlap) and that true independent performance is random chance. This data leakage pattern — where supplementary "validation" sets share sequences with training sets — may affect other published peptide property predictors and warrants systematic audit.

The replacement Treg TCR-contact criterion addresses the actual biological question: whether a peptide's TCR-facing residues are common in the human proteome, predicting engagement of natural regulatory T cells. This criterion identified a specific false positive (ENPVVHFFKNIVTPR) that the broken IL-10 model was promoting, demonstrating that criterion quality affects construct composition, not just numerical scores.

### Limitations

All outputs are computational hypotheses requiring experimental validation. The Treg TCR-contact criterion uses proteome-level frequency as a proxy for natural Treg cross-reactivity; the full JanusMatrix algorithm additionally considers HLA-restricted analysis. B-cell epitope risk uses the Parker hydrophilicity scale rather than structure-based prediction. The IFN-gamma model (IFNepitope2) was audited and found clean but is trained on general epitope data, not ITP- or MS-specific datasets. The pipeline currently scores only the primary antigen per disease; comprehensive vaccine design would include multiple antigens (e.g., PLP1 and MOG for MS). Peptide 2's tolerogenic mechanism in ITP operates through thymic escape and FoxP3+ Treg induction (Hall et al. 2019) without requiring high MHC affinity — a pathway that no purely computational pipeline can rank by design.

---

## 5. Data and code availability

The Tolerogenic Epitope Design Toolkit is available at https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit under MIT license. All input data derives from public sources (UniProt, IEDB, Nagpal et al. 2017 supplementary tables). Disease profiles for ITP and MS are included. The complete pipeline reproduces from empty data directories in approximately 9 minutes (ITP) or 3 minutes (MS), dominated by IEDB API prediction calls. Adding a new autoimmune disease requires two JSON configuration files.

---

## References

Clemente-Casares X et al. (2016). Expanding antigen-specific regulatory networks to treat autoimmunity. *Nature* 530:434-440.

De Groot AS et al. (2008). Activation of natural regulatory T cells by IgG Fc-derived peptide "Tregitopes." *Blood* 112(7):3303-3311.

Dhall A et al. (2024). IFNepitope2: improved prediction of interferon-gamma inducing peptides. *Scientific Reports*.

Evavold BD, Allen PM (1991). Induction of T-cell anergy by altered T-cell-receptor ligand on live antigen-presenting cells. *Nature* 356:604-607.

Hall LS et al. (2019). Combination peptide immunotherapy suppresses antibody and helper T-cell responses to GPIIb/IIIa in HLA-transgenic mice. *Haematologica* 104(5):1079-1087.

Kappos L et al. (2000). Induction of a non-encephalitogenic type 2 T helper autoimmune response in MS after administration of an altered peptide ligand. *Nature Medicine* 6:1176-1182.

Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157:105-132.

Larche M, Wraith DC (2005). Peptide-based therapeutic vaccines for allergic and autoimmune diseases. *Nature Medicine* 11:S69-76.

Liu GY et al. (1995). Low avidity recognition of self-antigen by T cells permits escape from central tolerance. *Immunity* 3:407-415.

Livingston BD et al. (2002). A rational strategy to design multiepitope immunogens based on multiple Th lymphocyte epitopes. *Journal of Immunology* 168:5499-5506.

Moise L et al. (2013). The two-faced T cell epitope: examining the host-microbe interface with JanusMatrix. *Human Vaccines and Immunotherapeutics* 9(7):1577-1586.

Nagpal G et al. (2017). Computer-aided designing of immunosuppressive peptides based on IL-10 inducing potential. *Scientific Reports* 7:42851.

Parker JMR, Guo D, Hodges RS (1986). New hydrophilicity scale derived from HPLC peptide retention data. *Biochemistry* 25:5425-5432.

Reynisson B et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation. *Nucleic Acids Research* 48:W449-454.

Serra P, Santamaria P (2019). Antigen-specific therapeutic approaches for autoimmunity. *Nature Biotechnology* 37:238-251.

Streeter HB et al. (2015). Preclinical development and first-in-human study of ATX-MS-1467 for immunotherapy of MS. *Journal of Autoimmunity* 65:104-111.

Sukati H et al. (2007). Mapping helper T-cell epitopes on platelet membrane glycoprotein IIIa in chronic autoimmune thrombocytopenic purpura. *Blood* 109(10):4528-4538.
