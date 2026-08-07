# STRICT-ADAPTIVE

**Stringent, Threshold-Regulated Identification of Consistent Taxa with ADAPTIVE variance scaling** - a multi-criterion framework for identifying placebo-nonresponsive gut microbial taxa, with leave-one-subject-out (LOSO) reproducibility analysis.

This repository accompanies the manuscript:

> Karapetyan G., Pepoyan A. *A Reproducible Framework for Identifying Placebo-Nonresponsive Gut Microbial Taxa.*

## Overview

The workflow identifies gut microbial taxa that remain stable ("placebo-nonresponsive") under placebo exposure, using paired pre/post PhyloChip™ G3 microarray profiles from individuals with familial Mediterranean fever (FMF) and healthy controls. Rather than reporting a single candidate list, it assesses internal reproducibility by LOSO resampling and reports a **stable core** - taxa retained across repeated removal of individual participants. LOSO agreement is treated as an internal reproducibility/stability check, not as independent, external validation (see manuscript Methods 2.5 and Limitations).

The pipeline performs:

1. **Preprocessing** - paired difference Δ = post − pre per OTU, per participant, on the RMA-normalised intensity scale.
2. **STRICT-ADAPTIVE criterion** - a joint, sex-stratified filter combining paired Wilcoxon testing (BH-FDR), Cliff's δ effect size (|δ| < 0.147), and adaptive, healthy-referenced variance/IQR/prevalence constraints. A taxon must be stable simultaneously in FMF and healthy cohorts.
3. **Duplication diagnosis** - confirms that cohort-specific selections differ, so the joint criterion is a deliberate, stricter choice.
4. **CLR robustness** - checks that the between-cohort effect-size ordering (Cliff's δ) is preserved under a compositional (centred-log-ratio) transform.
5. **Threshold justification** - parameter sweep over Cliff's δ and the variance scaling constant; consensus-curve check that the LOSO threshold (0.7) sits on a stable plateau for the male core.
6. **LOSO reproducibility analysis + stable core** - leave-one-subject-out resampling, consensus curve, and extraction of the reproducible core (≥ 0.7 selection frequency). The female core is retained for transparency but flagged as threshold-sensitive and not on the same evidential footing as the male core (see manuscript Results 3.3).
7. **Genus profile trees** - quantitative branch-length metrics with a tip-count-controlled permutation test (replaces qualitative "compactness" claims).
8. **Figures** - Δ overview, consensus curve, sex-overlap Venn, and taxonomic composition of the male core.
9. **Benchmark** - (a) lightweight CLR-based paired Wilcoxon/t-test sensitivity analyses, distinct from the ALDEx2/ANCOM-BC software; (b) an independent, platform-appropriate benchmark against **limma** (the canonical microarray differential-expression method - see `r_benchmark/`), comparing raw limma p-values for core vs. non-core taxa via a rank-based (Mann-Whitney) test.

## Key results (full data)

| Quantity | Value |
|---|---|
| Total OTUs analysed | 18,724 |
| Candidate placebo-nonresponsive (joint) | female 1,919 (10.2%); male 861 (4.6%) |
| **Reproducible stable core (LOSO ≥ 0.7)** | **female 146 (8%, threshold-sensitive); male 514 (60%, stable plateau)**; shared 3 |
| LOSO median Jaccard | female 0.33; male 0.47 |
| CLR robustness (Spearman ρ, between-cohort effect-size ordering raw vs CLR) | female 0.71; male 0.79 |
| limma rank-based benchmark (core vs non-core raw p-values, Mann-Whitney) | p = 4.3×10⁻⁶ to 1.5×10⁻³⁷ across all four sex×cohort strata |
| Genus tree shape (male vs female) | no significant difference (all permutation p > 0.05) |

The larger female candidate list reflects **lower statistical power** (n = 4–5 pairs), not greater stability. The male core is the primary result of this study; the female core is retained for transparency but is not treated as having the same evidential status (see manuscript Results 3.3 for the full discussion).

## Repository structure

```
strict-adaptive/
├── README.md
├── requirements.txt
├── .gitignore
├── notebook/
│   └── STRICT_ADAPTIVE.ipynb        # full, re-runnable Python analysis
├── r_benchmark/
│   ├── run_limma.R                  # step 1: limma paired moderated t-tests
│   └── compare_core_pvals.R         # step 2: rank-based benchmark (headline numbers, Results 3.5)
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

## Requirements

**Python 3.10+** for the main pipeline:
```bash
pip install -r requirements.txt
```

**R** (with the `limma` Bioconductor package) for the independent benchmark in `r_benchmark/`. `limma` is a lightweight, commonly pre-installed Bioconductor package; if needed:
```r
BiocManager::install("limma")
```
The R scripts read `pairs_only_P_H.csv` (not `.xlsx`) to avoid an `readxl`/`rlang` binary-compatibility issue we encountered on some R installations. Export the CSV once from Python:
```python
import pandas as pd
pd.read_excel("data/pairs_only_P_H.xlsx").to_csv("pairs_only_P_H.csv", index=False)
```

## Running

**Main pipeline:**
```bash
jupyter notebook notebook/STRICT_ADAPTIVE.ipynb
```
Run cells top to bottom. Section 0 installs pinned packages; Section 1 auto-detects the input file (`.xlsx` or `.csv`). The core criterion and LOSO reproducibility analysis are the slowest steps (a few minutes each on the full 18,724-taxon matrix). All outputs - CSV tables and figures, including `sa_outputs/stable_core_female.csv` and `stable_core_male.csv` - are written to `sa_outputs/`.

**limma benchmark** (after the notebook has produced `sa_outputs/`):
```bash
cd r_benchmark
Rscript run_limma.R            # step 1: fits limma models, writes limma_toptable_*.csv
Rscript compare_core_pvals.R   # step 2: rank-based comparison - reproduces the Results 3.5 numbers
```

## Notes on methodology

- The criterion is **joint across cohorts**: healthy samples serve as a reference within a single per-sex analysis, so the procedure yields one candidate set per sex. Per-cohort (healthy/FMF) versions of a given sex therefore contain the same taxa with cohort-specific metric columns; they are provided for transparency and are not independent selections.
- LOSO agreement is an **internal reproducibility/stability check**: at these sample sizes, any two LOSO folds necessarily share most of their data (75–80% of participants in the female strata), so fold-to-fold agreement shows the candidate set is not driven by any single participant, but does not by itself rule out pipeline-level bias reproduced across folds.
- The genus neighbour-joining trees are **profile trees** built from standardised stability features, not molecular phylogenies; their branching is not bootstrap-supported and is used only descriptively.
- ALDEx2 and ANCOM-BC are not run directly (they require sequencing read counts / library sizes, which PhyloChip intensities do not have); `limma`, the standard microarray differential-expression method, is used instead as the platform-appropriate contemporary benchmark.

## Citation

If you use this workflow, please cite the associated manuscript and the underlying software (NumPy, SciPy, pandas, statsmodels, scikit-bio, Biopython, limma).

## License

Not yet specified.
