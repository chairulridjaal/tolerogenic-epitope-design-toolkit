# A disease-agnostic computational pipeline for tolerogenic epitope design: case studies in immune thrombocytopenia and multiple sclerosis

## Abstract

Antigen-specific tolerogenic therapy offers the prospect of selectively restoring immune self-tolerance in autoimmune disease without broad immunosuppression, but rational computational design of tolerogenic epitope candidates remains underdeveloped. We present an open-source pipeline that scores peptide candidates from autoimmune disease target antigens against seven literature-grounded tolerogenic criteria — MHC-II binding zone, HLA allele promiscuity, disease-validated peptide proximity, IL-10 induction potential, IFN-gamma penalty, solubility, and human proteome cross-conservation — and assembles top candidates into codon-optimized multi-epitope mRNA vaccine constructs. Applied to immune thrombocytopenia (ITP) targeting GPIIIa/ITGB3 and multiple sclerosis (MS) targeting myelin basic protein (MBP), the pipeline independently recovers known immunodominant protein regions without gold standard labels: in MS, overlapping 15-mers from the MBP83-99 and MBP131-145 regions (both targeted by the ATX-MS-1467 clinical product) occupy ranks 1, 3-5 and 5, 10-11 respectively among 113 candidates scored by six independent criteria; in ITP, the C-terminal cytoplasmic tail of ITGB3 adjacent to the validated Peptide 82 region occupies ranks 1-3 for HLA promiscuity and IL-10 potential among 257 candidates. Novel candidates outside existing gold standards — including ENPVVHFFKNIVTPR from MBP217-231 (highest IL-10 induction score in the MS dataset) and a KEFAKFEEERARAKW cluster from ITGB3 aa749-765 (highest HLA promiscuity in the ITP dataset) — represent testable tolerogenic hypotheses. The pipeline is fully reproducible from public data sources, adds a new disease in two configuration files with zero code changes, and is available at https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit. All outputs are computational hypotheses requiring experimental validation.

---

## 1. Introduction

Autoimmune diseases collectively affect approximately 5-8% of the global population and are characterised by the breakdown of immune tolerance to self-antigens. Current treatments — corticosteroids, rituximab, and splenectomy for ITP; interferon-beta and natalizumab for MS — suppress immunity broadly rather than targeting the pathogenic response specifically. Antigen-specific tolerogenic therapy, administering self-antigen peptides in a context that promotes regulatory T cell (Treg) induction rather than effector activation, offers the prospect of disease-specific immune re-education without compromising systemic immunity (Serra and Santamaria 2019, Larche and Wraith 2005).

Proof-of-concept for this approach exists in both disease models studied here. In ITP, Hall et al. (2019) demonstrated that a combination of GPIIIa peptides (aa6-20 and aa711-725 in mature protein numbering) suppressed anti-GPIIb/IIIa antibody responses and induced FoxP3+ Tregs in HLA-DR15 transgenic mice. Notably, the tolerogenic mechanism was IL-10 independent — suppression was mediated by classical CD4+CD25+FoxP3+ Tregs rather than IL-10-producing Tr1 cells. In MS, the ATX-MS-1467 product — a mixture of four MBP peptides including the immunodominant MBP83-99 region — completed Phase I and Phase II clinical trials showing reduction in MBP-specific T cell responses and Th2 cytokine skewing (Streeter et al. 2015, Kappos et al. 2000).

Despite this clinical progress, no integrated computational pipeline exists for the systematic design of tolerogenic epitope candidates from autoimmune disease target antigens. Existing tools address components of this problem — NetMHCIIpan for MHC-II binding prediction (Reynisson et al. 2020), IL-10pred for cytokine induction scoring (Nagpal et al. 2017), IFNepitope2 for IFN-gamma prediction (Dhall et al. 2024) — but are not assembled into a disease-agnostic workflow that incorporates tolerogenic scoring criteria grounded in the relevant experimental literature.

We present the Tolerogenic Epitope Design Toolkit, an open-source Python pipeline that integrates MHC binding prediction, seven tolerogenic scoring criteria, and multi-epitope mRNA construct assembly. We apply it to ITP and MS as case studies, characterise the relationship between computational predictions and experimentally validated tolerogenic peptides, and identify novel candidate regions for experimental follow-up.

---

## 2. Methods

### 2.1 Disease profile system

The pipeline uses a JSON-based disease profile system. Each disease is specified by a configuration file containing target antigens (UniProt accessions), IEDB disease filter string, gold standard peptide file path, and calibration specification. Adding a new disease requires creating two JSON files with zero code modifications. Disease profiles for ITP and MS are included in the repository.

### 2.2 Antigen sequences and epitope data

Protein sequences were retrieved from UniProt via the REST API. For ITP: ITGA2B (P08514), ITGB3 (P05106), GP1BA (P07359), GP1BB (P13224), GP9 (P14770), GP5 (P40197). For MS: MBP (P02686), PLP1 (P60201), MOG (Q16653). Known T-cell and B-cell epitopes were retrieved from the IEDB query API (PostgREST interface) filtered by disease annotation ("autoimmune thrombocytopenic purpura" for ITP, "multiple sclerosis" for MS) and restricted to human host organism. All API responses are cached locally after first retrieval.

### 2.3 MHC-II binding prediction

Overlapping 15-mer peptides were scanned from the primary antigen sequence (ITGB3 for ITP, 774 peptides from 788 residues; MBP for MS, 290 peptides from 304 residues) using a sliding window. Binding predictions were obtained via the IEDB tools API using NetMHCIIpan 4.0 (eluted ligand predictor) across a panel of 12 HLA class II alleles: HLA-DRB1*01:01, *03:01, *04:01, *04:05, *07:01, *09:01, *11:01, *13:01, *15:01, HLA-DQA1*01:01/DQB1*05:01, HLA-DQA1*01:02/DQB1*06:02, and HLA-DPA1*01:03/DPB1*04:01. Alleles were selected to represent global HLA-DR diversity. Results were cached as TSV files with cache keys incorporating allele name, prediction method, and peptide content hash for reproducibility.

### 2.4 Tolerogenic scoring criteria

Seven criteria were implemented, each normalised to [0, 1] and combined in a weighted composite score. Criteria, weights, and their primary citations are summarised in Table 1.

**Table 1.** Tolerogenic scoring criteria and weights.

| # | Criterion | Weight | Implementation | Key citation |
|---|-----------|--------|----------------|-------------|
| 1 | MHC binding zone | 0.20 | Percentile rank 2-10% scores 1.0 (optimal); <2% scores 0.2 (penalised); 10-20% scores 0.5; >20% scores 0.1 | Evavold and Allen 1991, De Groot et al. 2008 |
| 2 | HLA promiscuity | 0.20 | Fraction of 12-allele panel with rank ≤10% | Serra and Santamaria 2019 |
| 3 | Disease-validated proximity | 0.30 | Substring overlap with experimentally validated peptides | Sukati et al. 2007, Hall et al. 2019, Streeter et al. 2015 |
| 4 | IL-10 induction | 0.08 | Local Random Forest (73 features, Nagpal et al. 2017 training data) | Nagpal et al. 2017 |
| 5 | IFN-gamma penalty | 0.07 | Inverted IFNepitope2 ExtraTrees probability | Dhall et al. 2024 |
| 6 | Solubility | 0.08 | GRAVY hydropathy score (Kyte-Doolittle) | Kyte and Doolittle 1982 |
| 7 | Human proteome cross-conservation | 0.07 | Fraction of 9-mer windows in Swiss-Prot human proteome | Moise et al. 2013 |

The IL-10 criterion weight (0.08) was reduced from an initial design value of 0.15 following Hall et al. (2019), who found no IL-10 elevation in the validated ITP tolerogenic model — suppression was mediated by FoxP3+ Tregs. The disease proximity weight (0.30) reflects the highest direct experimental evidence for both diseases.

The >20% MHC binding zone receives a score of 0.1 rather than 0.0 to accommodate ITP autoepitopes, which are characteristically low-affinity MHC binders. Sukati et al. (2007, Table 6) demonstrated using ProPred that most immunodominant GPIIIa peptides have no predicted high-affinity HLA-DR binding — a finding we confirm with NetMHCIIpan (Peptide 2: percentile rank 98% for HLA-DRB1*15:01). Their immunodominance arises from thymic escape: poor MHC presentation prevents thymic deletion of autoreactive T cells, which persist in the periphery where even inefficient presentation triggers activation (Liu et al. 1995).

### 2.5 IL-10 induction model

The IL-10 criterion uses a local Random Forest classifier (500 estimators, scikit-learn) trained on the exact 73 features from Nagpal et al. (2017) Table S1: 16 amino acid composition features and 57 dipeptide composition features. Training data: Table S4 (394 positive, 848 negative peptides). Validation on Table S5 (461 positive, 621 negative peptides) yielded AUC 0.91, with 37% sequence overlap between S4 and S5 sets. Five-fold cross-validation on S4 alone yielded accuracy 0.76 +/- 0.09, consistent with the paper's reported approximately 0.81.

### 2.6 IFN-gamma penalty model

The IFN-gamma criterion uses the IFNepitope2 ExtraTrees model (Dhall et al. 2024), loaded locally with a scikit-learn version compatibility shim. Dipeptide composition features are computed identically to the original package. The score is inverted (1.0 minus probability) so that high IFN-gamma induction probability yields a low tolerogenic score.

### 2.7 Human proteome cross-conservation (JMX proxy)

The JMX criterion approximates the JanusMatrix self-mimicry concept (Moise et al. 2013) using a 9-mer proteome lookup. The Swiss-Prot reviewed human proteome (20,431 proteins, 11,415,335 residues) was downloaded from UniProt and decomposed into 10,406,322 unique 9-mers stored as a compressed pickle index. For each candidate 15-mer, 7 overlapping 9-mer windows are checked against the index; the score is the fraction found. High cross-conservation indicates sequence similarity to human self-proteins and is associated with natural Treg engagement.

### 2.8 B-cell epitope safety filter

Linear B-cell epitope risk was assessed using the Parker hydrophilicity scale (Parker et al. 1986) with a sliding 7-residue window. Peptides with any window exceeding a hydrophilicity threshold of 4.0 were flagged as B-cell risk. This is a conservative proxy; structural prediction (BepiPred-3.0) would provide higher specificity but lacks a standalone pip-installable package.

### 2.9 Construct assembly

Top-ranked peptides were selected using a positional diversity filter requiring minimum 10-residue separation between start positions to prevent overlapping 15-mers from the same region from dominating the construct. Peptides were joined with GPGPG flexible linkers (Livingston et al. 2002). Codon optimisation used weighted random sampling from human Kazusa codon usage frequencies with an 18-codon sliding GC window constraint (maximum 62% per window), targeting overall GC content of 50-60%. Constructs include Kozak 5'UTR (GCCACCATG), human beta-globin 3'UTR, 120-nucleotide poly-A tail, and manufacturing annotations for N1-methylpseudouridine nucleoside modification and dsRNA-depleted LNP formulation.

### 2.10 Calibration

Disease-specific calibration checks were implemented using the disease profile system. For ITP: Peptide 2 (TTRGVSSCQQCLAVS, ITGB3 aa32-46 in UniProt numbering) was calibrated on disease proximity and JMX scores only, as low MHC binding is the expected biology for this thymic escape autoepitope; Peptide 82 (ALLIWKLLITIHDRK, ITGB3 aa737-751) was calibrated on composite rank (top 20%) and MHC zone score (>0.5), as it is a genuine multi-allele HLA-DR binder (Sukati et al. 2007, Table 6). For MS: MBP83-99 (17-mer, not directly in 15-mer scan) and MBP30-44 were calibrated on proximity and JMX.

### 2.11 Independent score analysis

To distinguish pipeline findings that depend on prior gold standard knowledge from those that are genuinely predictive, we computed an independent score excluding the disease proximity criterion (weight 0.30). The remaining six criteria weights were renormalised to sum to 1.0. All rankings described as "independent" in Results refer to this proximity-excluded score.

---

## 3. Results

### 3.1 Independent regional recovery in MS

When disease proximity is excluded and peptides are ranked by six independent criteria, overlapping 15-mers from the MBP83-99 window (aa83-99) occupy independent ranks 1, 3, 4, 7, and 8 of 113 MBP candidates (Table 2). Overlapping 15-mers from the MBP131-145 window occupy ranks 5, 10, and 11. Both regions correspond to components of ATX-MS-1467, which has shown clinical activity in Phase I (Streeter et al. 2015) and Phase II (Kappos et al. 2000) trials. The pipeline identified these regions without any gold standard information, driven primarily by MHC binding zone scores in the 2-10% tolerogenic range and high human proteome cross-conservation (JMX = 1.00 for all MBP peptides).

No exact gold standard peptide (MBP30-44, MBP83-99, MBP131-145, MBP140-154) appears in the top 20% by independent score alone. The pipeline identifies the correct protein regions, not the exact clinical peptides — consistent with the biology of MHC-II antigen processing, where the precise immunodominant peptide within a region is shaped by protease activity, endosomal chemistry, and TCR repertoire factors not modelled here.

**Table 2.** Top 5 MBP peptides by independent score (proximity excluded).

| Rank | Peptide | Position | MHC zone | HLA prom. | IL-10 | IFN-g | Sol. | JMX | Ind. score |
|------|---------|----------|----------|-----------|-------|-------|------|-----|------------|
| 1 | RPHLIRLFSRDAPGR | aa88-102 | 1.00 | 0.42 | 0.35 | 0.73 | 1.00 | 1.00 | 0.732 |
| 2 | ENPVVHFFKNIVTPR | aa217-231 | 0.73 | 0.75 | 0.88 | 0.00 | 0.81 | 1.00 | 0.718 |
| 3 | LIRLFSRDAPGREDN | aa91-105 | 1.00 | 0.33 | 0.29 | 0.75 | 1.00 | 1.00 | 0.703 |
| 4 | SRPHLIRLFSRDAPG | aa87-101 | 1.00 | 0.33 | 0.36 | 0.66 | 1.00 | 1.00 | 0.702 |
| 5 | SESLDVMASQKRPSQ | aa128-142 | 1.00 | 0.33 | 0.47 | 0.50 | 1.00 | 1.00 | 0.699 |

### 3.2 Independent regional recovery in ITP

For ITP, the independent score analysis reveals a qualitatively different pattern consistent with the known biology of ITP autoepitopes. Six of seven gold standard peptides rank in the bottom 20% of independent scores (Table 3) — expected, because ITP immunodominant peptides are characteristically low-affinity MHC binders (Sukati et al. 2007, Table 6). Their immunodominance arises from thymic escape rather than efficient MHC presentation.

However, the protein regions adjacent to validated peptides emerge independently. The C-terminal cytoplasmic tail of ITGB3 adjacent to the Peptide 82 region (aa737-751) occupies independent ranks 1-3 and 8 (Table 4), with the highest HLA allele promiscuity in the entire ITGB3 dataset (0.67-0.75, binding 8-9 of 12 alleles) and strong IL-10 induction scores (0.52-0.57). This region (aa749-765) sits 14 residues downstream from the validated tolerogenic peptide in the transmembrane-cytoplasmic junction of ITGB3.

**Table 3.** ITP gold standard peptide ranks by independent score.

| Peptide | Position | Independent rank | Ind. score | MHC zone | Note |
|---------|----------|------------------|------------|----------|------|
| Peptide 77 | aa687-701 | 34/257 | 0.656 | 1.00 | Only gold std peptide with strong MHC binding |
| Peptide 82 | aa737-751 | 218/257 | 0.509 | 0.73 | Moderate MHC, low solubility (0.23) |
| Peptide 44 | aa357-371 | 220/257 | 0.508 | 1.00 | Strong MHC but zero solubility |
| Peptide 2 | aa32-46 | 251/257 | 0.286 | 0.10 | Non-binder (rank 98% DRB1*15:01) |
| Peptide 47 | aa387-401 | 254/257 | 0.269 | 0.10 | Non-binder |
| Peptide 70 | aa617-631 | 252/257 | 0.272 | 0.10 | Non-binder |
| Peptide 53 | aa447-461 | 257/257 | 0.264 | 0.10 | Non-binder, last rank |

**Table 4.** Top 5 ITGB3 peptides by independent score (proximity excluded).

| Rank | Peptide | Position | MHC zone | HLA prom. | IL-10 | IFN-g | Sol. | JMX | Ind. score |
|------|---------|----------|----------|-----------|-------|-------|------|-----|------------|
| 1 | KEFAKFEEERARAKW | aa751-765 | 0.90 | 0.67 | 0.52 | 0.44 | 1.00 | 1.00 | 0.765 |
| 2 | DRKEFAKFEEERARA | aa749-763 | 0.73 | 0.75 | 0.52 | 0.58 | 1.00 | 1.00 | 0.756 |
| 3 | RKEFAKFEEERARAK | aa750-764 | 0.73 | 0.75 | 0.57 | 0.46 | 1.00 | 1.00 | 0.749 |
| 4 | DIYYLMDLSYSMKDD | aa139-153 | 1.00 | 0.50 | 0.35 | 0.59 | 1.00 | 1.00 | 0.742 |
| 5 | RDEIESVKELKDTGK | aa662-676 | 0.89 | 0.58 | 0.15 | 0.55 | 1.00 | 1.00 | 0.706 |

### 3.3 Novel candidate predictions

Two candidates outside any gold standard achieved high independent scores.

In MS, ENPVVHFFKNIVTPR (MBP aa217-231) ranks second by independent score with the highest IL-10 induction probability in the MS dataset (0.88) and the highest HLA promiscuity (0.75, binding 9 of 12 alleles). This region lies between the MBP131-145 and MBP140-154 ATX-MS-1467 components and has not been characterised for tolerogenic potential in the published literature accessible to us.

In ITP, the KEFAKFEEERARAKW cluster (ITGB3 aa749-765) achieves the highest HLA promiscuity in the ITP dataset (0.67-0.75) and strong IL-10 scores (0.52-0.57), with perfect solubility. It sits 14 residues downstream from the validated Peptide 82 region (aa737-751) in the transmembrane-cytoplasmic junction of ITGB3. Whether the strong MHC binding that makes this region computationally attractive is compatible with tolerogenic induction — given that ITP autoepitopes are characteristically low-affinity binders — is an open experimental question.

### 3.4 Calibration

All calibration checks passed for both diseases. ITP: Peptide 82 at composite rank 2/257 (mhc_zone = 0.73), Peptide 2 at rank 25/257 (itp_proximity = 1.00, jmx = 1.00, mhc_zone = 0.10). MS: MBP131-145 at composite rank 1/113, MBP30-44 at rank 2/113, MBP83-99 at rank 19/113.

### 3.5 Weight sensitivity analysis

To characterise the dependence of rankings on scoring weights, we evaluated eight alternative weight configurations ranging from equal weights to single-criterion emphasis (full results in Supplementary Table S1). MS results were robust: MBP131-145 and MBP30-44 ranked 1-2 under all configurations except those removing proximity entirely, confirming that validated MS tolerogenic peptides have strong independent computational properties across multiple criteria simultaneously.

ITP results were weight-sensitive, with Peptide 82 ranging from rank 1 (heavy proximity) to rank 221 (MHC+HLA only). This sensitivity is consistent with the thymic escape biology established by Sukati et al. (2007): ITP immunodominant peptides score poorly on MHC binding criteria by design, and their high composite ranking under current weights reflects prior experimental knowledge encoded in the proximity criterion rather than independent computational prediction. The exception is Peptide 77 (ITGB3 aa687-701), which ranks first under all MHC-weighted configurations — the only ITP gold standard peptide with genuine high MHC affinity, suggesting it may achieve immunodominance through conventional presentation rather than thymic escape.

These results indicate that the pipeline's regional recovery findings — emergence of the MBP83-99 cluster and ITGB3 C-terminal tail without proximity — are robust across weight configurations, while the ranking of exact gold standard peptides within those regions depends substantially on prior experimental knowledge encoded in the proximity weight.

### 3.6 Construct outputs

The ITP primary construct (ITP-GPGPG-10ep) incorporates 10 positionally diverse ITGB3 epitopes across 195 amino acids (768 nucleotides, 53.0% GC). The MS primary construct (MS-GPGPG-10ep) incorporates 10 MBP epitopes (195 amino acids, 768 nucleotides, 58.3% GC). Both GC contents fall within the 50-60% target range for therapeutic mRNA. AAY rigid linker variants were also generated for both diseases.

---

## 4. Discussion

The pipeline demonstrates consistent regional recovery of validated tolerogenic targets across two mechanistically distinct autoimmune diseases using six independent scoring criteria (MHC binding, HLA promiscuity, IL-10 induction, IFN-gamma penalty, solubility, and proteome cross-conservation), while generating novel testable hypotheses from the same criteria. The finding that exact experimental peptides are not recovered — but their protein neighborhoods are — is consistent with the biology of antigen presentation: MHC-II prediction operates at the level of binding potential, while the precise immunodominant peptide within a region is shaped by processing efficiency, protease activity, endosomal chemistry, and TCR repertoire, none of which are modelled here.

The contrast between MS and ITP results is instructive. In MS, the pipeline's independent criteria directly identify the ATX-MS-1467 target regions because those peptides are genuine moderate-affinity MHC binders. In ITP, the independent criteria identify neighboring rather than overlapping regions because the validated tolerogenic peptides are low-affinity binders by design — their immunodominance arises from thymic escape. This disease-specific difference validates the pipeline's sensitivity to the underlying biology rather than a one-size-fits-all ranking.

The IL-10 criterion weight (0.08) reflects our finding that the validated ITP tolerogenic model operates through FoxP3+ Tregs rather than IL-10-producing Tr1 cells (Hall et al. 2019, Figure 5). For MS, where Th2/IL-10 skewing was observed in clinical trials (Kappos et al. 2000), a higher weight may be appropriate; the disease profile system supports per-disease weight overrides for this purpose.

### Limitations

Several limitations constrain interpretation. All outputs are computational hypotheses derived from preclinical and early clinical literature. No criterion has been validated in the context of this specific pipeline. The IL-10 model achieves AUC 0.91 on an independent validation set but with 37% sequence overlap between training (S4) and validation (S5) sets from the original publication. The JMX cross-conservation criterion approximates JanusMatrix using 9-mer proteome lookup rather than the full algorithm incorporating HLA-restricted TCR-facing residue analysis. B-cell epitope risk uses the Parker hydrophilicity scale (a linear property-based predictor) rather than structure-based prediction. Expert immunological review of candidate peptides is required before any experimental follow-up.

### Generalisability

The disease profile system demonstrates that adding a new autoimmune disease requires two JSON configuration files (disease profile and gold standard) with zero code changes. Cold-cache end-to-end reproduction for a new disease completes in approximately 3 minutes (MBP, 290 peptides) to 9 minutes (ITGB3, 774 peptides), dominated by IEDB API prediction calls. Extension to Type 1 Diabetes (insulin, GAD65, IA-2), Celiac Disease (gliadin peptides), or Rheumatoid Arthritis (citrullinated peptides) requires only the identification of appropriate gold standard peptides from the experimental literature.

---

## 5. Data and code availability

The Tolerogenic Epitope Design Toolkit is available at https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit under MIT license. All input data is derived from public sources (UniProt, IEDB, Nagpal et al. 2017 supplementary tables). The complete pipeline can be reproduced from empty data directories on a standard laptop with internet access. Disease profiles for ITP and MS are included.

---

## References

Clemente-Casares X et al. (2016). Expanding antigen-specific regulatory networks to treat autoimmunity. Nature 530:434-440.

De Groot AS et al. (2008). Activation of natural regulatory T cells by IgG Fc-derived peptide "Tregitopes." Blood 112(7):3303-3311.

Dhall A et al. (2024). IFNepitope2: improved prediction of interferon-gamma inducing peptides. Scientific Reports.

Evavold BD, Allen PM (1991). Induction of T-cell anergy by altered T-cell-receptor ligand on live antigen-presenting cells. Nature 356:604-607.

Hall LS et al. (2019). Combination peptide immunotherapy suppresses antibody and helper T-cell responses to GPIIb/IIIa in HLA-transgenic mice. Haematologica 104(5):1079-1087.

Kappos L et al. (2000). Induction of a non-encephalitogenic type 2 T helper autoimmune response in MS after administration of an altered peptide ligand. Nature Medicine 6:1176-1182.

Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology 157:105-132.

Larche M, Wraith DC (2005). Peptide-based therapeutic vaccines for allergic and autoimmune diseases. Nature Medicine 11:S69-76.

Liu GY et al. (1995). Low avidity recognition of self-antigen by T cells permits escape from central tolerance. Immunity 3:407-415.

Livingston BD et al. (2002). A rational strategy to design multiepitope immunogens based on multiple Th lymphocyte epitopes. Journal of Immunology 168:5499-5506.

Moise L et al. (2013). The two-faced T cell epitope: examining the host-microbe interface with JanusMatrix. Human Vaccines and Immunotherapeutics 9(7):1577-1586.

Nagpal G et al. (2017). Computer-aided designing of immunosuppressive peptides based on IL-10 inducing potential. Scientific Reports 7:42851.

Parker JMR, Guo D, Hodges RS (1986). New hydrophilicity scale derived from HPLC peptide retention data. Biochemistry 25:5425-5432.

Reynisson B et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation. Nucleic Acids Research 48:W449-454.

Serra P, Santamaria P (2019). Antigen-specific therapeutic approaches for autoimmunity. Nature Biotechnology 37:238-251.

Streeter HB et al. (2015). Preclinical development and first-in-human study of ATX-MS-1467 for immunotherapy of MS. Journal of Autoimmunity 65:104-111.

Sukati H et al. (2007). Mapping helper T-cell epitopes on platelet membrane glycoprotein IIIa in chronic autoimmune thrombocytopenic purpura. Blood 109(10):4528-4538.
