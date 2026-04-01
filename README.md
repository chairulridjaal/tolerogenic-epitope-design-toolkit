# Tolerogenic Epitope Design Toolkit

A computational pipeline for designing antigen-specific tolerogenic vaccines, starting with Immune Thrombocytopenic Purpura (ITP).

Current treatments for ITP — steroids, IVIg, splenectomy, TPO receptor agonists — suppress the immune system broadly or compensate for platelet destruction without addressing the root cause: autoantibodies targeting platelet surface glycoproteins. This toolkit takes a different approach. It identifies the specific peptide fragments (epitopes) of platelet antigens most likely to induce immune tolerance rather than immune activation, then assembles them into candidate mRNA constructs for tolerogenic vaccine design.

The architecture is disease-agnostic. ITP is the first target. The same pipeline generalizes to Multiple Sclerosis, Type 1 Diabetes, Celiac Disease, Lupus, Rheumatoid Arthritis, and other autoimmune conditions by swapping the input antigen configuration.

> **Status:** Phases 0–3 complete. The full tolerogenic scoring engine is live, 100% offline, and calibrated against experimentally validated ITP epitopes (Peptide 2 + Peptide 82 from Hall et al. 2019). Real IL-10 and IFN-γ induction predictions are now integrated. Phase 4 (construct assembly) is next.

---

## The Problem

In ITP, B cells produce autoantibodies against platelet surface glycoproteins — primarily **GPIIb/IIIa** (integrin alphaIIb-beta3) and the **GPIb-IX-V complex**. These autoantibodies flag platelets for destruction by splenic macrophages. The bone marrow keeps producing platelets; the immune system keeps destroying them. Platelet counts fall from a normal range of 150,000-400,000/uL to below 20,000, sometimes below 10,000, where spontaneous bleeding becomes life-threatening.

The underlying failure is a breakdown in **immune tolerance** — the learned non-response to self-antigens. Both central tolerance (thymic deletion of self-reactive cells) and peripheral tolerance (Treg-mediated suppression) have failed for these specific antigens.

## The Approach

A tolerogenic vaccine presents self-antigens in an immunological context that promotes regulatory T cell (Treg) activity rather than effector T cell activation. The critical design question is: **which peptide fragments of the target antigen, presented on which MHC Class II molecules, are most likely to induce tolerance?**

This toolkit automates that question through four completed stages and one in progress:

1. **Data retrieval** — Fetch protein sequences from UniProt and known epitope data from IEDB
2. **Epitope prediction** — Scan antigen sequences for peptides that bind MHC Class II molecules across diverse HLA alleles (using NetMHCIIpan for T-cell epitopes and BepiPred for B-cell epitopes)
3. **Tolerogenic scoring** — Literature-grounded ranking using seven criteria, including **local retrained Random Forest model for IL-10 induction** (73 features from Nagpal et al. 2017) and **local IFNepitope2 model for IFN-γ penalty** (Dhall et al. 2024). Calibrated against the seven experimentally validated GPIIIa epitopes from Sukati et al. (2007) and Hall et al. (2019). 100% offline. (complete)
4. **Construct assembly** — Combine top-ranked epitopes with appropriate linkers into candidate multi-epitope mRNA constructs

### Why MHC Variation Matters

Different people carry different HLA alleles, which means different peptides from the same antigen get displayed on their cell surfaces. The autoimmune response in one ITP patient may be driven by entirely different epitopes than in another. The pipeline accounts for this by running predictions across a panel of common HLA-DRB1, HLA-DQB1, and HLA-DPB1 alleles and computing population coverage statistics. Candidate constructs are selected to maximize coverage across real-world HLA frequency distributions.

---

## Project Structure

```
tolerogenic-epitope-design-toolkit/
├── src/
│   ├── data/           # API connectors for UniProt and IEDB
│   ├── prediction/     # Epitope prediction wrappers (NetMHCIIpan, BepiPred)
│   ├── scoring/        # Tolerogenic scoring module
│   └── assembly/       # Multi-epitope mRNA construct builder
├── data/
│   ├── raw/            # Downloaded sequences and epitope databases
│   └── processed/      # Cleaned tables, parsed binding data
├── notebooks/          # Analysis notebooks for each pipeline stage
├── docs/               # Background reading, design decisions, references
├── tests/              # Unit and integration tests
├── requirements.txt
└── README.md
```

## Installation

```bash
git clone https://github.com/chairulridjaal/tolerogenic-epitope-design-toolkit.git
cd tolerogenic-epitope-design-toolkit
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

### Dependencies

Core dependencies (specified in `requirements.txt`):

- `requests` — API calls to UniProt and IEDB
- `biopython` — Sequence parsing and manipulation
- `pandas` — Tabular data handling
- `numpy` — Numerical operations
- `matplotlib` / `seaborn` — Visualization
- `pytest` — Testing
- `scikit-learn`, `joblib` — local IL-10 & IFN-γ models

External tools called by the prediction module:

- [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) — MHC Class II binding prediction (requires academic license)
- [BepiPred 3.0](https://services.healthtech.dtu.dk/services/BepiPred-3.0/) — B-cell epitope prediction

---

## Roadmap

### Phase 0 — Foundation
- Repository structure and development environment
- Background reading on tolerogenic vaccine design
- Conceptual documentation of the pipeline

### Phase 1 — Data Layer
- UniProt API connector: fetch GPIIb/IIIa (ITGA2B + ITGB3) and GPIb-IX-V sequences
- IEDB API connector: retrieve known ITP-associated epitopes as benchmark data
- Data parsing, cleaning, and local storage
- Notebook: antigen sequence landscape

### Phase 2 — Epitope Prediction Engine
- NetMHCIIpan integration for CD4+ T-cell epitope prediction
- BepiPred integration for B-cell epitope prediction
- Peptide scanning across target antigens for a panel of common HLA alleles
- Benchmark against known ITP epitopes from Phase 1
- Notebook: epitope landscape of GPIIb/IIIa

### Phase 3 — Tolerogenic Scoring
- Local IL-10 Random Forest model (73 features from original paper)
- Local IFNepitope2 integration for IFN-γ penalty
- Gold-standard calibration against real ITP tolerogenic peptides
- Full composite scoring with population coverage

### Phase 4 — Construct Assembly (CURRENT)
- Multi-epitope linker design from vaccine literature
- mRNA sequence generation with nucleoside modification annotations
- End-to-end pipeline: antigen sequence in, ranked mRNA constructs out
- Notebook: full pipeline demonstration

### Phase 5 — Generalization
- Add a second autoimmune disease (e.g., Multiple Sclerosis with MBP, PLP, MOG antigens) to validate disease-agnostic architecture
- Streamlit web interface for non-programmers
- Outreach to immunology researchers for collaboration
- Preprint or methods paper

---

## Target Antigens (ITP)

| Protein | UniProt ID | Role | Relevance |
|---------|-----------|------|-----------|
| Integrin alpha-IIb (GPIIb) | P08514 | Platelet adhesion and aggregation | Primary autoantibody target in most ITP cases |
| Integrin beta-3 (GPIIIa) | P05106 | Forms heterodimer with GPIIb | Primary autoantibody target in most ITP cases |
| Glycoprotein Ib alpha (GPIb-alpha) | P07359 | Initial platelet adhesion to vessel walls | Secondary autoantibody target |
| Glycoprotein IX | P14770 | Part of GPIb-IX-V receptor complex | Secondary autoantibody target |
| Glycoprotein V | P40197 | Part of GPIb-IX-V receptor complex | Secondary autoantibody target |

## Generalizable to Other Autoimmune Diseases

The pipeline accepts any self-antigen as input. Planned extensions:

| Disease | Target Antigens |
|---------|----------------|
| Multiple Sclerosis | MBP, PLP, MOG |
| Type 1 Diabetes | Insulin, GAD65, IA-2 |
| Celiac Disease | Gluten-derived gliadin peptides |
| Systemic Lupus Erythematosus | dsDNA, Smith antigen, Ro/La |
| Rheumatoid Arthritis | Citrullinated peptides (CCP) |

---

## Background

For readers unfamiliar with the immunology:

**Immune tolerance** is the learned ability of the immune system to not attack the body's own tissues. It is maintained by central mechanisms (thymic deletion of self-reactive lymphocytes) and peripheral mechanisms (suppression by regulatory T cells). Autoimmune diseases occur when tolerance breaks down for specific self-antigens.

**Epitopes** are short peptide fragments (typically 9-25 amino acids) of a protein that can be recognized by the immune system. Not every fragment of a protein is immunologically relevant — only those that physically fit into the binding groove of an MHC molecule get presented to T cells.

**MHC Class II molecules** (called HLA in humans) present extracellular antigen fragments to CD4+ T cells. The shape of the binding groove varies between HLA alleles, so different people present different peptides from the same protein. This is why population coverage analysis is essential for vaccine design.

**Tolerogenic vaccines** present self-antigens in a context that promotes regulatory T cell activity rather than effector T cell activation. The immunological context — absence of danger signals, specific adjuvant formulations, route of administration — determines whether the immune response is tolerogenic or immunogenic.

---

## Contributing

This project is in early development. If you are an immunologist, computational biologist, or developer interested in tolerogenic vaccine design, please open an issue or reach out.

Areas where collaboration would be especially valuable:
- Validation of tolerogenic scoring criteria against experimental data
- Access to patient-derived epitope binding data
- Clinical perspective on tolerogenic vaccine feasibility for ITP

## License

TBD

## References

Key literature informing this project (to be expanded in `docs/`):

1. Provan, D. et al. "International consensus report on the investigation and management of primary immune thrombocytopenia." *Blood* (2010).
2. Cines, D.B. & Bussel, J.B. "How I treat idiopathic thrombocytopenic purpura (ITP)." *Blood* (2005).
3. Serra, P. & Santamaria, P. "Antigen-specific therapeutic approaches for autoimmunity." *Nature Biotechnology* (2019).
4. Clemente-Casares, X. et al. "Expanding antigen-specific regulatory networks to treat autoimmunity." *Nature* (2016).
5. Reynisson, B. et al. "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation." *Nucleic Acids Research* (2020).
6. Nagpal et al. (2017) *Scientific Reports* — IL-10pred (basis for local RF model)
7. Dhall et al. (2024) *Scientific Reports* — IFNepitope2