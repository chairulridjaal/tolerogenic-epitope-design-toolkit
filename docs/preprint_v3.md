# A disease-agnostic computational pipeline for tolerogenic epitope selection and mRNA construct design: validation across immune thrombocytopenia and multiple sclerosis

## Abstract

Antigen-specific tolerogenic therapy offers selective restoration of immune self-tolerance without broad immunosuppression, but no integrated computational pipeline exists for rational design of tolerogenic epitope candidates from arbitrary autoimmune disease target antigens. We present an open-source pipeline that scores peptide candidates against six weighted criteria — MHC-II binding zone, HLA allele promiscuity, disease-validated peptide proximity, Treg TCR-contact self-similarity, IFN-gamma penalty, and solubility — supplemented by cathepsin S processing likelihood and population coverage statistics, then assembles top candidates into codon-optimised multi-epitope mRNA vaccine constructs.

Applied to immune thrombocytopenia (ITP, targeting ITGB3/GPIIIa, 257 scored candidates) and multiple sclerosis (MS, targeting MBP, 113 candidates), the pipeline independently recovers known immunodominant protein regions using five criteria that require no gold standard knowledge. In MS, all five top independently-ranked peptides are overlapping 15-mers from the MBP83-99 immunodominant region targeted by the ATX-MS-1467 clinical product. In ITP, the ITGB3 C-terminal cytoplasmic tail (aa749-765) adjacent to the validated Peptide 82 occupies independent ranks 1-2 and 4, with the highest Treg TCR-contact self-similarity (0.99) and HLA promiscuity (0.67-0.75) in the dataset. Cathepsin S processing scores provide additional peptide-level resolution: among overlapping 15-mers from the ITGB3 Peptide 77 region, the experimentally validated peptide has the highest processing score (0.55 vs 0.02-0.19 for neighbours), and among MBP83-99 overlapping 15-mers, the top independently-ranked peptide also has the highest processing score (0.63).

During development, we discovered and reported a data leakage problem in IL-10pred (Nagpal et al. 2017): 36.6% sequence overlap between the published training and validation sets inflates the reported AUC from 0.50 (true independent performance) to 0.91. ESM-2 protein language model embeddings do not rescue generalisation, confirming the limitation is in the training data. We replaced this criterion with Treg TCR-contact self-similarity, a non-ML score measuring the frequency of TCR-facing residues (MHC-II binding core positions P2, P3, P5, P7, P8) in the human proteome. This replacement identified and corrected a false positive: ENPVVHFFKNIVTPR (MBP aa217-231), previously the top-ranked novel MS candidate, was correctly demoted from rank 2 to rank 67 based on its foreign-like TCR contacts.

Constructs include population coverage estimates (ITP: European 90%, East Asian 74%, African 71%, South Asian 82%; MS: European 80%, East Asian 62%, African 60%, South Asian 69%). The pipeline adds a new disease in two JSON configuration files with zero code changes and reproduces from empty data directories in approximately nine minutes. Code and data are available at https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit. All outputs are computational hypotheses requiring experimental validation.

---

## 1. Introduction

Autoimmune diseases collectively affect 5-8% of the global population and are characterised by breakdown of immune tolerance to self-antigens. Current treatments suppress immunity broadly: corticosteroids and splenectomy for ITP, interferon-beta and natalizumab for MS. Antigen-specific tolerogenic therapy — administering self-antigen peptides in a context promoting regulatory T cell induction rather than effector activation — offers disease-specific immune re-education without compromising systemic immunity (Serra and Santamaria 2019, Larche and Wraith 2005).

Clinical proof-of-concept exists for both diseases studied here. In ITP, Hall et al. (2019) demonstrated that GPIIIa peptides aa32-46 and aa737-751 (UniProt P05106 precursor numbering) suppressed anti-GPIIb/IIIa antibody responses and induced FoxP3+ Tregs in HLA-DR15 transgenic mice. The tolerogenic mechanism was IL-10 independent — suppression was mediated by classical CD4+CD25+FoxP3+ regulatory T cells (Hall et al. 2019, Figure 5). In MS, ATX-MS-1467, a mixture of four MBP peptides including the immunodominant MBP83-99 region, completed Phase I and Phase II clinical trials demonstrating reduction in MBP-specific T cell responses and Th2 cytokine skewing (Streeter et al. 2015, Kappos et al. 2000).

Despite this progress, no integrated computational pipeline exists for systematic tolerogenic epitope design. Existing tools address components — NetMHCIIpan for MHC-II binding (Reynisson et al. 2020), IFNepitope2 for cytokine prediction (Dhall et al. 2024), JanusMatrix for cross-conservation (Moise et al. 2013) — but are not assembled into a disease-agnostic workflow incorporating tolerogenic scoring criteria grounded in experimental literature. We present the Tolerogenic Epitope Design Toolkit, apply it to ITP and MS, characterise the relationship between computational predictions and validated tolerogenic peptides, report the discovery of a data leakage problem in a widely-used cytokine prediction model, and identify novel candidate regions for experimental follow-up.

---

## 2. Methods

### 2.1 Disease profile system and data sources

Each disease is specified by a JSON configuration file containing target antigens, IEDB disease filter, gold standard peptide path, and calibration specification. Adding a new disease requires two files with no code changes. Protein sequences were retrieved from UniProt; known epitopes from the IEDB query API filtered by disease and human host. Gold standard peptide positions are verified against UniProt sequences at startup, with signal peptide offset detection.

### 2.2 MHC-II binding prediction

Overlapping 15-mers were scanned from the primary antigen (774 for ITGB3, 290 for MBP). Binding predictions were obtained via the IEDB tools API using NetMHCIIpan 4.0 (eluted ligand predictor) across 12 HLA class II alleles covering global diversity. Peptides with percentile rank ≤10% for at least one allele were retained; gold standard peptides were additionally included regardless of binding rank.

### 2.3 Tolerogenic scoring criteria

Six criteria, each normalised to [0, 1], are combined in a weighted composite score (Table 1). Two additional metrics — cathepsin S processing likelihood and population coverage — are computed as informational columns.

**Table 1.** Scoring criteria and weights.

| # | Criterion | Weight | Implementation | Citation |
|---|-----------|--------|----------------|----------|
| 1 | MHC binding zone | 0.20 | Percentile rank 2-10% optimal; <2% penalised; >20% scores 0.1 | De Groot 2008 |
| 2 | HLA promiscuity | 0.20 | Fraction of 12-allele panel with rank ≤10% | Serra 2019 |
| 3 | Disease proximity | 0.30 | Substring overlap with validated peptides | Sukati 2007, Hall 2019 |
| 4 | Treg TCR-contact | 0.15 | Frequency of TCR-facing 5-mer motif in human proteome | De Groot 2008, Moise 2013 |
| 5 | IFN-gamma penalty | 0.07 | Inverted IFNepitope2 probability | Dhall 2024 |
| 6 | Solubility | 0.08 | GRAVY hydropathy score | Kyte 1982 |
| — | Processing | info | Cathepsin S cleavage site likelihood (PSSM) | Choe 2006, Biniossek 2011 |
| — | Population coverage | info | Bui formula with AFND allele frequencies | Bui 2006 |

The >20% MHC binding zone receives 0.1 rather than 0.0 to accommodate ITP autoepitopes, which are characteristically low-affinity MHC binders whose immunodominance arises from thymic escape (Sukati et al. 2007, Liu et al. 1995).

### 2.4 Treg TCR-contact self-similarity (Criterion 4)

In MHC-II presentation, the 9-mer binding core has anchor residues (P1, P4, P6, P9) facing the MHC and TCR-contact residues (P2, P3, P5, P7, P8) facing the T cell receptor. Tregitopes are characterised by TCR contacts highly conserved in the human proteome (De Groot et al. 2008, Moise et al. 2013). We computed a frequency index from the Swiss-Prot human proteome (20,431 proteins, 10,406,322 unique 9-mers): for each 9-mer, the 5 TCR-facing residues were extracted and counted (2,356,207 unique motifs, 73.6% of the 3.2M possible 5-mer space). For each candidate, the binding core(s) from NetMHCIIpan were retrieved, TCR motifs extracted, and frequency converted to a percentile score.

This criterion replaced an IL-10 induction model found to have no independent predictive power (see Section 3.5). It requires no ML model and carries no generalisation risk.

We evaluated whether the full 9-mer JanusMatrix proxy could provide additional signal, but all three approaches tested — binary presence, protein-count frequency, and sequence similarity at ≥7/9 identity — produced zero or near-zero variance for human self-protein peptides. The 9-mer sequence space (20^9 ≈ 5×10^11) is too sparse relative to the human proteome's 10.4M unique 9-mers for any metric to discriminate. The 5-mer TCR motif space (20^5 = 3.2M, 73.6% covered) is dense enough for frequency-based discrimination.

### 2.5 Cathepsin S processing likelihood

To address the gap between regional recovery and exact peptide prediction, we scored each 15-mer by the cathepsin S cleavage likelihood at its N-terminal and C-terminal boundaries. Cathepsin S is the dominant protease in the MHC-II pathway in dendritic cells. A position-specific scoring matrix was derived from published substrate profiling data (Choe et al. 2006, Biniossek et al. 2011, MEROPS C01.034): the P2 position dominates specificity (Leu 35%, Val 26%, Asp 0%). Scores are normalised as percentiles within each antigen. This metric is reported as an informational column, not a weighted criterion.

### 2.6 Population coverage

Population coverage was computed using the Bui et al. (2006) formula with allele frequencies from the Allele Frequency Net Database for European, East Asian, African, and South Asian populations. Coverage = 1 − ∏(1 − p_phenotype_i), where p_phenotype = 1 − (1 − p_allele)² and the product runs over all alleles binding at least one construct epitope.

### 2.7 B-cell safety filter and construct assembly

B-cell epitope risk was assessed using the Parker hydrophilicity scale with protein-relative thresholding (top 20th percentile for the source protein). Top-ranked peptides were selected with a positional diversity filter (minimum 10-residue separation). Constructs use GPGPG flexible linkers, human Kazusa codon-optimised mRNA with 18-codon sliding GC window constraint (50-60% target), Kozak 5'UTR, beta-globin 3'UTR, 120nt poly-A tail, with annotations for N1-methylpseudouridine and dsRNA-depleted LNP formulation.

### 2.8 Calibration

Disease-specific calibration checks are defined in the disease profile. ITP: Peptide 82 (genuine HLA-DR binder) must rank in top 20% with mhc_zone > 0.5; Peptide 2 (thymic escape autoepitope) must achieve proximity ≥ 1.0 and jmx ≥ 1.0 regardless of composite rank. MS: MBP83-99 and MBP30-44 must achieve proximity ≥ 0.4 and jmx ≥ 0.5.

### 2.9 Independent score analysis

To separate findings dependent on prior knowledge from genuinely predictive results, we computed an independent score excluding Criterion 3 (disease proximity, weight 0.30). The remaining five criteria were renormalised to sum to 1.0.

---

## 3. Results

### 3.1 Calibration

All calibration checks passed (Table 2). In ITP, Peptide 77 ranked first — the only gold standard peptide with both high MHC affinity (zone 1.00) and strong Treg TCR-contact self-similarity (0.79). The gold standard splits into two mechanistic classes: three conventional MHC binders (ranks 1, 2, 4) and four thymic escape peptides (ranks 87-211), consistent with Sukati et al. (2007) Table 6. In MS, MBP131-145 and MBP30-44 (ATX-MS-1467 components) ranked first and second, with Treg TCR-contact scores of 0.95 and 0.93.

**Table 2.** Gold standard peptide rankings.

| Disease | Peptide | Rank | Composite | Treg TCR | MHC | Processing | Validated |
|---------|---------|------|-----------|----------|-----|------------|-----------|
| ITP | Peptide_77 | 1/257 | 0.718 | 0.79 | 1.00 | 0.55 | No |
| ITP | Peptide_82 | 2/257 | 0.644 | 0.60 | 0.73 | 0.75 | Hall 2019 |
| ITP | Peptide_44 | 4/257 | 0.611 | 0.70 | 1.00 | 0.79 | No |
| ITP | Peptide_2 | 87/257 | 0.472 | 0.50 | 0.10 | 0.93 | Hall 2019 |
| MS | MBP131-145 | 1/113 | 0.785 | 0.95 | 1.00 | 0.17 | ATX-MS-1467 |
| MS | MBP30-44 | 2/113 | 0.776 | 0.93 | 1.00 | 0.42 | ATX-MS-1467 |
| MS | MBP140-154 | 22/113 | 0.517 | 0.50 | 0.13 | 0.23 | ATX-MS-1467 |
| MS | MBP83-99 | 25/113 | 0.508 | 0.50 | 0.10 | 0.01 | Phase I/II |

### 3.2 Independent regional recovery — MS

When disease proximity is excluded, all five top-ranked MBP peptides are overlapping 15-mers from the MBP83-99 immunodominant window (Table 3), driven by MHC binding zone scores in the tolerogenic range and high Treg TCR-contact self-similarity (0.83-0.96). No exact gold standard peptide appears in the top 20% by independent score.

**Table 3.** Top 5 MBP peptides by independent score.

| Rank | Peptide | Position | MHC | HLA | Treg TCR | Ind. score |
|------|---------|----------|-----|-----|----------|------------|
| 1 | RPHLIRLFSRDAPGR | aa88-102 | 1.00 | 0.42 | 0.83 | 0.771 |
| 2 | LIRLFSRDAPGREDN | aa91-105 | 1.00 | 0.33 | 0.94 | 0.771 |
| 3 | IRLFSRDAPGREDNT | aa92-106 | 1.00 | 0.25 | 0.96 | 0.758 |
| 4 | HLIRLFSRDAPGRED | aa90-104 | 1.00 | 0.25 | 0.96 | 0.751 |
| 5 | SRPHLIRLFSRDAPG | aa87-101 | 1.00 | 0.33 | 0.83 | 0.740 |

### 3.3 Independent regional recovery — ITP

Six of seven ITP gold standard peptides rank in the bottom 20% of independent scores — expected for thymic escape autoepitopes (Sukati et al. 2007). The ITGB3 C-terminal tail (aa749-765) occupies independent ranks 1, 2, and 4, with Treg TCR-contact scores of 0.99 and HLA promiscuity of 0.67-0.75. The Peptide 77 region (aa661-676) occupies ranks 3 and 5.

**Table 4.** Top 5 ITGB3 peptides by independent score.

| Rank | Peptide | Position | MHC | HLA | Treg TCR | Ind. score |
|------|---------|----------|-----|-----|----------|------------|
| 1 | KEFAKFEEERARAKW | aa751-765 | 0.90 | 0.67 | 0.99 | 0.819 |
| 2 | DRKEFAKFEEERARA | aa749-763 | 0.73 | 0.75 | 0.99 | 0.809 |
| 3 | RDEIESVKELKDTGK | aa662-676 | 0.89 | 0.58 | 0.99 | 0.801 |
| 4 | RKEFAKFEEERARAK | aa750-764 | 0.73 | 0.75 | 0.99 | 0.797 |
| 5 | CRDEIESVKELKDTG | aa661-675 | 0.87 | 0.50 | 0.99 | 0.781 |

### 3.4 Cathepsin S processing resolves peptide clusters

The processing score provides peptide-level resolution within high-scoring regions. Among the four overlapping 15-mers from the ITP Peptide 77 cluster (aa687-702), the experimentally validated peptide (DDCVVRFQYYEDSSG) has the highest processing score: 0.55 versus 0.19, 0.02, and 0.16 for its neighbours. Among the five overlapping MBP83-99 15-mers, RPHLIRLFSRDAPGR (independently ranked first) also has the highest processing score (0.63 versus 0.15, 0.50, 0.06, 0.51).

ITP Peptide 2 (a thymic escape autoepitope with MHC zone 0.10) has a processing score of 0.93 — the highest among all gold standard peptides. This is biologically coherent: poor MHC binding allows autoreactive T cells to escape thymic deletion, but the peptide IS efficiently generated by cathepsin S processing, ensuring it reaches the cell surface at low but immunologically relevant concentrations.

### 3.5 IL-10 model data leakage

The pipeline originally used an IL-10 induction model (Nagpal et al. 2017) as Criterion 4. We discovered that 396 sequences (36.6% of the validation set S5) are shared with the training set S4, and all 396 are positive examples. On the 686 truly non-overlapping S5 sequences (65 positive, 621 negative), AUC is 0.50 — random chance. ESM-2 protein language model embeddings (esm2_t6_8M_UR50D, 320 dimensions, tested with logistic regression, random forest, and MLP classifiers) yield AUC 0.51, confirming the limitation is in the training data rather than the feature representation. This data leakage pattern — where supplementary "validation" sets share sequences with training sets — may affect other published peptide property predictors.

The replacement Treg TCR-contact criterion identified a specific false positive: ENPVVHFFKNIVTPR (MBP aa217-231), which had the highest IL-10 score (0.88), was demoted from rank 2 to rank 67 based on its foreign-like TCR contacts (score 0.56 versus median 0.77).

### 3.6 Construct outputs and population coverage

**Table 5.** Multi-epitope mRNA construct summary.

| Construct | Epitopes | AA | mRNA | GC% | European | East Asian | African | South Asian |
|-----------|----------|-----|------|-----|----------|------------|---------|-------------|
| ITP-GPGPG-10ep | 10 | 195 | 768 nt | 53.7% | 90% | 74% | 71% | 82% |
| MS-GPGPG-10ep | 10 | 195 | 768 nt | 57.9% | 80% | 62% | 60% | 69% |

African coverage is lowest for both constructs, reflecting the higher HLA-DRB1 diversity in African populations relative to the 12-allele panel. Panel expansion with African-enriched alleles (e.g., DRB1*13:02, DRB1*11:04) would improve equity.

### 3.7 Weight sensitivity

MS results were robust: MBP131-145 and MBP30-44 ranked 1-2 under all eight weight configurations tested except those removing proximity entirely. ITP results were weight-sensitive (Peptide 82: rank 2 under current weights, rank 233 with no proximity), consistent with the thymic escape biology. Regional recovery was robust across all configurations.

---

## 4. Discussion

The pipeline demonstrates two modes of validated target recovery reflecting the distinct biology of each disease. In MS, five independent scoring criteria directly identify the ATX-MS-1467 target region because those peptides are genuine moderate-affinity MHC binders with self-like TCR contacts. In ITP, the same criteria identify neighbouring rather than overlapping regions because validated tolerogenic peptides are low-affinity binders whose immunodominance arises from thymic escape.

Cathepsin S processing scores provide a new dimension of peptide-level discrimination. The finding that experimentally validated peptides tend to have the highest processing scores within their clusters — DDCVVRFQYYEDSSG scoring 0.55 versus ≤0.19 for its neighbours, and RPHLIRLFSRDAPGR scoring 0.63 versus ≤0.51 for its neighbours — suggests that antigen processing efficiency contributes to immunodominance hierarchy within a region, even when MHC binding potential is similar across overlapping windows. This metric, derived from published cathepsin S substrate profiling data without machine learning, may have broader applicability in epitope prediction.

The IL-10 model audit represents a cautionary finding. The reported AUC of 0.91 (Nagpal et al. 2017) is entirely driven by 396 memorised sequences (36.6% overlap). This pattern may be widespread in peptide property prediction literature and warrants systematic audit of published training/validation set overlaps.

The Treg TCR-contact criterion addresses the question that IL-10 was meant to proxy: whether a peptide's TCR surface engages regulatory or effector T cells. By measuring proteome-level frequency of TCR-facing residue motifs at the 5-mer level (where sequence space density permits discrimination), rather than binary 9-mer presence (where it does not), the criterion achieves meaningful variance (0.28-1.00) for human self-protein peptides and identified a specific false positive that the broken IL-10 model was masking.

### Limitations

All outputs are computational hypotheses requiring experimental validation. The Treg TCR-contact criterion uses proteome frequency as a proxy for Treg cross-reactivity; the full JanusMatrix algorithm additionally considers HLA-restricted analysis. Cathepsin S processing scores use published substrate preferences, not a trained model, and do not account for protein 3D structure or competing substrates. B-cell risk uses Parker hydrophilicity rather than structural prediction. The IFN-gamma model (IFNepitope2) was audited and found clean but trained on general epitope data. The pipeline scores only the primary antigen per disease; comprehensive vaccine design would include multiple antigens. Individual patient TCR repertoire and in vivo tolerogenic efficacy are beyond computational prediction.

---

## 5. Data and code availability

The Tolerogenic Epitope Design Toolkit is available at https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit under MIT licence. All input data derives from public sources. Disease profiles for ITP and MS are included. The pipeline reproduces from empty data directories in approximately 9 minutes (ITP) or 3 minutes (MS).

---

## References

Biniossek ML et al. (2011). Identification of protease specificity by combining proteome-derived peptide libraries and quantitative proteomics. *Molecular and Cellular Proteomics* 10:M111.009613.

Bui HH et al. (2006). Predicting population coverage of T-cell epitope-based diagnostics and vaccines. *BMC Bioinformatics* 7:153.

Choe Y et al. (2006). Substrate profiling of cysteine proteases using a combinatorial peptide library identifies functionally unique specificities. *Journal of Biological Chemistry* 281:12824-12832.

Clemente-Casares X et al. (2016). Expanding antigen-specific regulatory networks to treat autoimmunity. *Nature* 530:434-440.

De Groot AS et al. (2008). Activation of natural regulatory T cells by IgG Fc-derived peptide "Tregitopes." *Blood* 112:3303-3311.

Dhall A et al. (2024). IFNepitope2: improved prediction of interferon-gamma inducing peptides. *Scientific Reports*.

Evavold BD, Allen PM (1991). Induction of T-cell anergy by altered T-cell-receptor ligand on live antigen-presenting cells. *Nature* 356:604-607.

Hall LS et al. (2019). Combination peptide immunotherapy suppresses antibody and helper T-cell responses to GPIIb/IIIa in HLA-transgenic mice. *Haematologica* 104:1079-1087.

Kappos L et al. (2000). Induction of a non-encephalitogenic type 2 T helper autoimmune response in multiple sclerosis. *Nature Medicine* 6:1176-1182.

Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157:105-132.

Larche M, Wraith DC (2005). Peptide-based therapeutic vaccines for allergic and autoimmune diseases. *Nature Medicine* 11:S69-76.

Liu GY et al. (1995). Low avidity recognition of self-antigen by T cells permits escape from central tolerance. *Immunity* 3:407-415.

Livingston BD et al. (2002). A rational strategy to design multiepitope immunogens. *Journal of Immunology* 168:5499-5506.

Moise L et al. (2013). The two-faced T cell epitope: JanusMatrix. *Human Vaccines and Immunotherapeutics* 9:1577-1586.

Nagpal G et al. (2017). Computer-aided designing of immunosuppressive peptides based on IL-10 inducing potential. *Scientific Reports* 7:42851.

Parker JMR et al. (1986). New hydrophilicity scale derived from HPLC peptide retention data. *Biochemistry* 25:5425-5432.

Reynisson B et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0. *Nucleic Acids Research* 48:W449-454.

Serra P, Santamaria P (2019). Antigen-specific therapeutic approaches for autoimmunity. *Nature Biotechnology* 37:238-251.

Streeter HB et al. (2015). Preclinical development and first-in-human study of ATX-MS-1467 for immunotherapy of MS. *Journal of Autoimmunity* 65:104-111.

Sukati H et al. (2007). Mapping helper T-cell epitopes on platelet membrane glycoprotein IIIa in chronic autoimmune thrombocytopenic purpura. *Blood* 109:4528-4538.
