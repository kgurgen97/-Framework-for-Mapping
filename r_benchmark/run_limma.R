#!/usr/bin/env Rscript
# ============================================================
# STRICT-ADAPTIVE / limma benchmark - step 1
#
# limma is the canonical differential-expression method for
# microarray intensity data and is used here as the
# platform-appropriate contemporary benchmark for PhyloChip
# intensities (Manuscript, Methods 2.6; Results 3.5), in place
# of count-based sequencing tools (ALDEx2, ANCOM-BC, MaAsLin2),
# whose bias-correction assumes a library size that PhyloChip
# intensities do not have.
#
# For each sex x cohort stratum: intensities are log2-transformed,
# paired post-minus-baseline contrasts are computed per
# participant, and a moderated one-sample t-test (limma eBayes)
# is fitted per taxon.
#
# INPUT: pairs_only_P_H.csv in the working directory (export the
# xlsx to CSV first, e.g. in Python: df.to_csv('pairs_only_P_H.csv', index=False)
# -- the readxl package can be incompatible with some R/Bioconductor
# environments; reading CSV via base R avoids that dependency entirely).
#
# OUTPUT: limma_toptable_<stratum>.csv (per-taxon results, 4 files)
#         limma_stable_female.txt, limma_stable_male.txt (joint
#         "not significantly changed" sets, for reference only --
#         the headline benchmark result uses the RANK-based
#         comparison in compare_core_pvals.R, not these threshold
#         lists; see note at the end of this script).
# ============================================================

suppressMessages(library(limma))

cat("Loading data...\n")
df <- read.csv("pairs_only_P_H.csv", check.names = FALSE)
cat("Rows:", nrow(df), "\n")

otu_id <- as.character(df$OTU_ID)

# participant sex groups (same IDs as the Python STRICT-ADAPTIVE pipeline)
P_FEMALE <- c(10, 36, 43, 44, 8)
P_MALE   <- c(11, 16, 21, 31, 39, 45, 47, 48, 49, 51)
H_FEMALE <- c(12, 22, 23, 27)
H_MALE   <- c(3, 7, 9, 14, 15, 16, 20, 21, 25, 29)

get_pair_cols <- function(prefix, ids) {
  out <- list()
  for (id in ids) {
    a <- paste0(prefix, id, "a")
    b <- paste0(prefix, id, "b")
    if (a %in% colnames(df) && b %in% colnames(df)) {
      out[[as.character(id)]] <- c(a, b)
    }
  }
  out
}

run_limma_paired <- function(pairs, label) {
  ids <- names(pairs)
  n <- length(ids)
  if (n < 3) {
    cat(sprintf("[%s] too few pairs (n=%d), skipping\n", label, n))
    return(NULL)
  }
  a_cols <- sapply(pairs, `[`, 1)
  b_cols <- sapply(pairs, `[`, 2)

  # log2 transform (RMA intensities are linear-scale here; floor at 1 before log)
  A <- log2(pmax(as.matrix(df[, a_cols]), 1))
  B <- log2(pmax(as.matrix(df[, b_cols]), 1))

  # paired difference per participant (rows = taxa, cols = participants)
  Delta <- B - A
  colnames(Delta) <- ids

  # moderated one-sample t-test: intercept-only design tests mean(Delta) != 0
  design <- matrix(1, nrow = n, ncol = 1)
  colnames(design) <- "Intercept"
  fit <- lmFit(Delta, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "Intercept", number = Inf, sort.by = "none")

  stable <- otu_id[tt$adj.P.Val >= 0.05]
  cat(sprintf("[%s] n_pairs=%d, taxa_tested=%d, limma-stable(FDR>=0.05)=%d (%.1f%%), min adj.P=%.4f, min raw P=%.4f\n",
              label, n, nrow(tt), length(stable), 100 * length(stable) / nrow(tt),
              min(tt$adj.P.Val, na.rm = TRUE), min(tt$P.Value, na.rm = TRUE)))

  list(stable = stable, table = tt)
}

pf <- get_pair_cols("P", P_FEMALE)
pm <- get_pair_cols("P", P_MALE)
hf <- get_pair_cols("H", H_FEMALE)
hm <- get_pair_cols("H", H_MALE)

cat("\n--- Running limma paired moderated t-tests ---\n")
res_pf <- run_limma_paired(pf, "FMF-female")
res_pm <- run_limma_paired(pm, "FMF-male")
res_hf <- run_limma_paired(hf, "Healthy-female")
res_hm <- run_limma_paired(hm, "Healthy-male")

female_limma_stable <- intersect(res_pf$stable, res_hf$stable)
male_limma_stable   <- intersect(res_pm$stable, res_hm$stable)
cat(sprintf("\n[reference only, not the headline result] Joint limma-stable at FDR>=0.05 (both cohorts): female=%d, male=%d\n",
            length(female_limma_stable), length(male_limma_stable)))
cat("Note: at these sample sizes, essentially all taxa pass the FDR>=0.05 gate in every stratum\n")
cat("(genome-wide correction across 18,724 tests is too conservative to be discriminating here).\n")
cat("The benchmark reported in the manuscript (Results 3.5) instead compares the RAW p-value\n")
cat("DISTRIBUTION of STRICT-ADAPTIVE core taxa against non-core taxa -- run compare_core_pvals.R\n")
cat("next, after this script, to reproduce those numbers.\n")

writeLines(female_limma_stable, "limma_stable_female.txt")
writeLines(male_limma_stable, "limma_stable_male.txt")
write.csv(res_pf$table, "limma_toptable_FMF_female.csv", row.names = TRUE)
write.csv(res_pm$table, "limma_toptable_FMF_male.csv", row.names = TRUE)
write.csv(res_hf$table, "limma_toptable_Healthy_female.csv", row.names = TRUE)
write.csv(res_hm$table, "limma_toptable_Healthy_male.csv", row.names = TRUE)

cat("\nDone. Files written to the current working directory.\n")
cat("Next: run compare_core_pvals.R to reproduce the manuscript's Results 3.5 benchmark numbers.\n")
