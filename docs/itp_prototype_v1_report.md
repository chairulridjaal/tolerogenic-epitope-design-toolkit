# ITP Tolerogenic mRNA Vaccine Prototype v1.0

## ITGB3 (GPIIIa) Multi-Epitope Constructs — Review Report

### Pipeline Summary

| Metric | Value |
|---|---|
| Source antigen | ITGB3 / GPIIIa (P05106), 788 aa |
| Peptides scored | 260 (253 MHC binders + 7 gold standard) |
| Scoring criteria | 7 (MHC zone, HLA promiscuity, ITP proximity, IL-10, IFN-γ, solubility, JMX) |
| Calibration | **PASS** — Peptide 82 rank 1/260, Peptide 2 rank 2/260 |
| Constructs generated | 2 variants (GPGPG and AAY linker) |

---

### Construct Summary Table

| Rank | Construct ID | Peptides | Adjusted Score | JMX Range | B-cell Risk | Tier | Key Notes |
|------|-------------|----------|----------------|-----------|-------------|------|-----------|
| 1 | ITP-ITGB3-GPGPG-5ep | Pep82, Pep2, Pep44, Pep70, Pep56 | 0.957 | 0.70–1.00 | 2/5 flagged | T1x2, T2x3 | Primary candidate. Flexible linker, 93 aa, 462 nt. +0.3 bonuses. |
| 2 | ITP-ITGB3-AAY-5ep | Pep82, Pep2, Pep44, Pep70, Pep56 | 0.957 | 0.70–1.00 | 2/5 flagged | T1x2, T2x3 | Rigid linker variant, 85 aa, 438 nt. Better cathepsin processing. |

---

### Detailed Breakdown

---

#### Construct 1: ITP-ITGB3-GPGPG-5ep (Primary Candidate)

**Full amino acid sequence (93 aa):**

```
NPIYKSAVTTVVNP-GPGPG-GDCNCTKDDSVMCIG-GPGPG-WNIQNPWSIKRKRK-GPGPG-GCVYIEDEVHRCYGR-GPGPG-TEEKIQLIVQANIDQ
```

```
NPIYKSAVTTVVNPGPGPGGDCNCTKDDSVMCIGGPGPGWNIQNPWSIKRKRKGPGPGGCVYIEDEVHRCYGRGPGPGTEEKIQLIVQANIDQ
```

**Codon-optimized mRNA sequence (462 nt):**

```
GCCACCATGAACCCCATCTACAAGAGCGCCGTGACCACCGTGGTGAACCCCGGCCCCGGC
CCCGGCGGCGACTGCAACTGCACCAAGGACGACAGCGTGATGTGCATCGGCGGCCCCGGC
CCCGGCTGGAACATCCAGAACCCCTGGAGCATCAAGCGGAAGCGGAAGGGCCCCGGCCCCG
GCGGCTGCGTGTACATCGAGGACGAGGTGCACCGGTGCTACGGCCGGGGCCCCGGCCCCG
GCACCGAGGAGAAGATCCAGCTGATCGTGCAGGCCAACATCGACCAGTGAGCTCGCTTTC
TTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAG + polyA(120)
```

**Manufacturing notes:**

- All U residues to 1-methylpseudouridine (m1 psi) during synthesis
- 5' cap: CleanCap AG or m7GpppAm
- dsRNA-depleted LNP formulation for spleen/liver targeting
- GC content: 69.9%

**Component scores:**

| Epitope | MHC Zone | HLA Prom. | ITP Prox. | IL-10 | IFN-g | Solubility | JMX | Composite | B-cell | Tier |
|---------|----------|-----------|-----------|-------|-------|------------|-----|-----------|--------|------|
| NPIYKSAVTTVVNP | 0.75 | 0.58 | 1.00 | 0.23 | 0.43 | 0.60 | **1.00** | **0.680** | Clean | **T1** |
| GDCNCTKDDSVMCIG | 0.75 | 0.58 | 1.00 | 0.12 | 0.51 | 0.76 | 0.85 | **0.673** | Risk | **T1** |
| WNIQNPWSIKRKRK | 0.65 | 0.50 | 0.80 | 0.58 | 0.41 | 1.00 | 0.70 | **0.661** | Clean | T2 |
| GCVYIEDEVHRCYGR | 0.65 | 0.50 | 0.80 | 0.42 | 0.56 | 1.00 | 0.70 | **0.647** | Risk | T2 |
| TEEKIQLIVQANIDQ | 0.65 | 0.50 | 0.80 | 0.29 | 0.59 | 0.92 | 0.70 | **0.624** | Clean | T2 |

**Construct-level bonuses:** +0.2 (multi-epitope, 3 or more from same antigen) + +0.1 (spatial clustering across multiple gold-standard regions). **Junction flags: 0** (no unintended strong MHC binders at linker junctions).

**Why this construct ranked #1:** Contains both tolerogenic-validated peptides (Peptide 82 and Peptide 2) alongside three additional immunodominant epitopes spanning the GPIIIa extracellular domain. GPGPG flexible linkers minimize junctional neo-epitope formation. The lead epitope (Peptide 82) has perfect human self-mimicry (JMX 1.0) and no B-cell risk.

**Recommendation:** Strong Tier 1 candidate — anchored by two experimentally validated tolerogenic peptides (Hall et al. 2019), with high self-mimicry, no junction epitope flags, and broad HLA coverage. Prioritize for preclinical testing in HLA-DR15 humanized mice.

---

#### Construct 2: ITP-ITGB3-AAY-5ep (Rigid Linker Variant)

**Full amino acid sequence (85 aa):**

```
NPIYKSAVTTVVNPAAYGDCNCTKDDSVMCIGAAYWNIQNPWSIKRKRKAAYGCVYIEDEVHRCYGRAAYTEEKIQLIVQANIDQ
```

**mRNA:** 438 nt, GC content 63.9%. Same component peptides, manufacturing notes, and scoring as Construct 1.

**Difference from Construct 1:** AAY rigid linker (vs GPGPG flexible). The AAY motif contains a cathepsin cleavage site that promotes intracellular processing and may enhance MHC-II loading efficiency in dendritic cells. Shorter total length (85 vs 93 aa) may improve translational efficiency.

**Recommendation:** Consider as an alternative if in vivo antigen processing of the GPGPG construct is suboptimal. The rigid AAY linker may improve epitope liberation in the endosomal compartment.

---

### Component Epitope Notes

**Epitope 1 — NPIYKSAVTTVVNP (Peptide 82, Tier 1)**

Confirmed tolerogenic in HLA-DR15 humanized mice (Hall et al. 2019). Suppressed anti-GPIIb/IIIa autoantibody production. Perfect JMX self-mimicry (1.0). Moderate IL-10 induction (0.23), low IFN-gamma risk (0.43 = tolerogenic-favorable). No B-cell epitope risk. The strongest overall candidate.

**Epitope 2 — GDCNCTKDDSVMCIG (Peptide 2, Tier 1)**

Confirmed Treg induction in HLA-DR15 mice. Strong solubility (0.76). B-cell risk flagged — the cysteine-rich sequence (3x Cys) may form surface-exposed disulfide knots that attract antibodies. Consider RK-flanking (arginine-lysine wrapper) to improve delivery, as used in the original Hall et al. study.

**Epitope 3 — WNIQNPWSIKRKRK (Peptide 44, Tier 2)**

Immunodominant in ITP patient T-cell recall assays (Sukati et al. 2007). Highest IL-10 induction score of all components (0.58). The RK-rich C-terminus provides excellent solubility (1.0). Strong tolerogenic profile despite lacking direct validation.

**Epitope 4 — GCVYIEDEVHRCYGR (Peptide 70, Tier 2)**

Immunodominant epitope from the GPIIIa cysteine-rich domain (aa 591-605). B-cell risk flagged due to hydrophilic Glu-Asp patch. Good IFN-gamma safety (0.56).

**Epitope 5 — TEEKIQLIVQANIDQ (Peptide 56, Tier 2)**

Immunodominant at aa 421-435. Balanced scoring across all criteria. Strong solubility (0.92), moderate IFN-gamma safety (0.59). No B-cell risk. Solid supporting epitope for bystander suppression.

---

### Calibration Confirmation

```
Peptide_82: Rank 1/260 (top 0.4%) — PASS
Peptide_2:  Rank 2/260 (top 0.8%) — PASS
Overall: PASS
```

Pipeline runs cleanly end-to-end. All 7 scoring criteria produce real variation. No neutral 0.5 placeholders remain for any criterion with available data.

---

### References

1. Sukati H et al. (2007). Mapping helper T-cell epitopes on platelet membrane glycoprotein IIIa in chronic autoimmune thrombocytopenic purpura. *Blood* 109(10):4528-4538.
2. Hall LS et al. (2019). Combination peptide immunotherapy suppresses antibody and helper T-cell responses to the major human platelet autoantigen glycoprotein IIb/IIIa in HLA-transgenic mice. *Haematologica* 104(5):1079-1087.
3. Nagpal G et al. (2017). Computer-aided designing of immunosuppressive peptides based on IL-10 inducing potential. *Scientific Reports* 7:42851.
4. Dhall A et al. (2024). IFNepitope2: improved prediction of interferon-gamma inducing peptides. *Scientific Reports*.
5. Moise L et al. (2013). The two-faced T cell epitope: examining the host-microbe interface with JanusMatrix. *Human Vaccines & Immunotherapeutics* 9(7):1577-1586.
6. Reynisson B et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation. *Nucleic Acids Research*.
