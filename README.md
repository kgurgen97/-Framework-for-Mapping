# STRICT-ADAPTIVE

**Stringent, Threshold-Regulated Identification of Consistent Taxa with ADAPTIVE variance scaling** — a multi-criterion framework for mapping placebo responsiveness across gut microbial taxa, with leave-one-subject-out (LOSO) validation.

This repository accompanies the manuscript:

> Karapetyan G., Pepoyan A. *A Multi-criterion Framework for Mapping Placebo Responsiveness across 18,724 Gut Microbial Taxa: a Reproducible, LOSO-validated Stability Core.*

## Overview

The workflow identifies gut microbial taxa that remain stable ("placebo-nonresponsive") under placebo exposure, using paired pre/post PhyloChip™ G3 microarray profiles from individuals with familial Mediterranean fever (FMF) and healthy controls. Rather than reporting a single candidate list, it validates reproducibility by LOSO resampling and reports a **stable core** — taxa that survive removal of any single participant.

The pipeline performs:

1. **Preprocessing** — paired difference Δ = post − pre per OTU, per sex, per cohort.
2. **STRICT-ADAPTIVE criterion** — a joint, sex-stratified filter combining paired Wilcoxon testing (BH-FDR), Cliff's δ effect size (|δ| < 0.147), and adaptive, healthy-referenced variance/IQR/prevalence constraints. A taxon must be stable simultaneously in FMF and healthy cohorts.
3. **Duplication diagnosis** — confirms that cohort-specific selections differ, so the joint criterion is a deliberate, stricter choice.
4. **CLR robustness** — checks that the stability ranking is preserved under a compositional (centred-log-ratio) transform.
5. **Threshold justification** — parameter sweep over Cliff's δ and the variance scaling constant.
6. **LOSO validation + stable core** — leave-one-subject-out resampling, consensus curve, and extraction of the reproducible core (≥ 0.7 selection frequency).
7. **Genus profile trees** — quantitative branch-length metrics with a tip-count-controlled permutation test (replaces qualitative "compactness" claims).
8. **Figures** — Δ overview, consensus curve, sex-overlap Venn, and taxonomic composition of the core.
9. **Benchmark** — Python equivalents of ALDEx2 and ANCOM-BC as compositionally-aware stress tests.

## Key results (full data)

| Quantity | Value |
|---|---|
| Total OTUs analysed | 18,724 |
| Candidate placebo-nonresponsive (joint) | female 1,919 (10.2%); male 861 (4.6%) |
| **Reproducible stable core (LOSO ≥ 0.7)** | **female 146 (8%); male 514 (60%)**; shared 3 |
| LOSO median Jaccard | female 0.33; male 0.47 |
| CLR robustness (Spearman ρ, Cliff's δ raw vs CLR) | female 0.71; male 0.79 |
| ALDEx2-style recall of the male core | 0.95 |
| Genus tree shape (male vs female) | no significant difference (all permutation p > 0.05) |

The larger female candidate list reflects **lower statistical power** (n = 4–5 pairs), not greater stability; the reproducible core is larger and denser in males.

## Repository structure

```
strict-adaptive/
├── README.md
├── requirements.txt
├── .gitignore
├── notebook/
│   └── STRICT_ADAPTIVE.ipynb        # full, re-runnable analysis
├── data/
│   ├── pairs_only_P_H.xlsx          # paired pre/post intensity matrix (main input)
│   ├── Dataset_1_All_taxa_studied.csv
│   ├── Dataset_2..5_*_placebo_resistant.csv
│   └── resistant_taxa_taxonomic_tree.newick
├── figures/                         # publication figures (300 dpi)
└── results/                         # created when the notebook runs (sa_outputs/)
```

## Data

Raw and normalised microarray data are publicly available in GEO under accession **GSE111835**.
The main analysis input (`data/pairs_only_P_H.xlsx`) contains, for each of 18,724 OTUs, an annotation block (`OTU_ID, Genus, Species, Phylum, Class, Order, Family, Lineage, rep_gene`) followed by paired intensity columns named `P##a/P##b` (FMF cohort) and `H##a/H##b` (healthy cohort), where `a` = before and `b` = after one month of placebo.

> Note: `pairs_only_P_H.xlsx` is ~17 MB. If GitHub upload is inconvenient, it can be regenerated from GSE111835; the CSV datasets are smaller.



## Running

```bash
jupyter notebook notebook/STRICT_ADAPTIVE.ipynb
```

Run cells top to bottom. Section 0 installs pinned packages; Section 1 auto-detects the input file (`.xlsx` or `.csv`). The core criterion (Section 2) and LOSO (Section 6) are the slowest steps (a few minutes each on the full 18,724-taxon matrix). All outputs — CSV tables and figures — are written to `sa_outputs/`.

## Notes on methodology

- The criterion is **joint across cohorts**: healthy samples serve as a reference within a single per-sex analysis, so the procedure yields one candidate set per sex. Per-cohort (healthy/FMF) versions of a given sex therefore contain the same taxa with cohort-specific metric columns; they are provided for transparency and are not independent selections.
- The genus neighbour-joining trees are **profile trees** built from standardised stability features, not molecular phylogenies; their branching is not bootstrap-supported and is used only descriptively.
- ANCOM-BC2 and ALDEx2 are implemented here as Python equivalents for benchmarking; because PhyloChip yields intensities rather than sequencing counts, these serve as stress tests, not gold standards.

## Citation

If you use this workflow, please cite the associated manuscript and the underlying software (NumPy, SciPy, pandas, statsmodels, scikit-bio, Biopython).

