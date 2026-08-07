#!/usr/bin/env Rscript
# ============================================================
# STRICT-ADAPTIVE / limma benchmark - step 2 (headline result)
#
# Rank-based comparison: does the STRICT-ADAPTIVE stable core
# have systematically higher (less-significant) limma raw
# p-values than non-core taxa? This is the benchmark reported
# in the manuscript (Results, Section 3.5) -- it avoids the
# degenerate FDR-threshold comparison (which is uninformative at
# n = 4-10 across 18,724 simultaneous tests; see run_limma.R)
# and instead asks whether the underlying evidence-of-change
# ranking is shifted for core taxa, using a one-sided
# Mann-Whitney test (H1: core taxa have larger raw p-values,
# i.e. less evidence of change).
#
# PREREQUISITES: run run_limma.R first (produces the
# limma_toptable_*.csv files read below). Also requires the
# STRICT-ADAPTIVE stable-core taxon lists exported from the
# Python notebook (sa_outputs/stable_core_female.csv and
# stable_core_male.csv, each with a "taxon" column of OTU IDs).
# ============================================================

cat("Loading limma per-taxon results and STRICT-ADAPTIVE cores...\n")

tt_pf <- read.csv("limma_toptable_FMF_female.csv")
tt_pm <- read.csv("limma_toptable_FMF_male.csv")
tt_hf <- read.csv("limma_toptable_Healthy_female.csv")
tt_hm <- read.csv("limma_toptable_Healthy_male.csv")

core_female <- read.csv("sa_outputs/stable_core_female.csv")
core_male   <- read.csv("sa_outputs/stable_core_male.csv")
core_female_ids <- as.character(core_female$taxon)
core_male_ids   <- as.character(core_male$taxon)

cat(sprintf("Female core: %d taxa | Male core: %d taxa\n", length(core_female_ids), length(core_male_ids)))

# recover OTU_ID in the same row order used by run_limma.R
df <- read.csv("pairs_only_P_H.csv", check.names = FALSE)
otu_id <- as.character(df$OTU_ID)

attach_ids <- function(tt) { tt$OTU_ID <- otu_id; tt }
tt_pf <- attach_ids(tt_pf)
tt_pm <- attach_ids(tt_pm)
tt_hf <- attach_ids(tt_hf)
tt_hm <- attach_ids(tt_hm)

rank_test <- function(tt, core_ids, label) {
  is_core <- tt$OTU_ID %in% core_ids
  p_core <- tt$P.Value[is_core]
  p_rest <- tt$P.Value[!is_core]
  wt <- wilcox.test(p_core, p_rest, alternative = "greater")
  cat(sprintf("[%s] n_core=%d n_rest=%d | median_p_core=%.4f median_p_rest=%.4f | Mann-Whitney p=%.3e\n",
              label, length(p_core), length(p_rest), median(p_core), median(p_rest), wt$p.value))
  invisible(wt)
}

cat("\n--- Rank-based test: core vs non-core limma raw p-values ---\n")
cat("(H1: core taxa have systematically LARGER p-values, i.e. less evidence of change)\n\n")
cat("These are the numbers reported in Results 3.5 of the manuscript.\n\n")

cat("FEMALE core vs FMF-female stratum:\n")
rank_test(tt_pf, core_female_ids, "female core vs FMF-female p")
cat("FEMALE core vs Healthy-female stratum:\n")
rank_test(tt_hf, core_female_ids, "female core vs Healthy-female p")

cat("\nMALE core vs FMF-male stratum:\n")
rank_test(tt_pm, core_male_ids, "male core vs FMF-male p")
cat("MALE core vs Healthy-male stratum:\n")
rank_test(tt_hm, core_male_ids, "male core vs Healthy-male p")

cat("\nDone. These four Mann-Whitney p-values are the ones reported in the manuscript.\n")
