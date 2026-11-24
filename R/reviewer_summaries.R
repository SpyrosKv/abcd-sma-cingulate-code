paper1_reviewer_summaries <- function(cfg) {
  baseline <- paper1_baseline(cfg)
  out <- file.path(cfg$work_dir, "Reviewer_Revision_Summaries")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  nct <- file.path(cfg$work_dir, "R2_NCT_10000")
  completion <- file.path(nct, "RESTRICTED", "COMPLETE.rds")
  if (!file.exists(completion)) stop("Complete the six validated r2_nct comparisons first.")
  complete <- readRDS(completion)
  if (!identical(unname(tools::md5sum(names(complete$checksums))), unname(complete$checksums)))
    stop("R2 result files changed after completion.")
  omni <- read.csv(file.path(nct, "21_NCT10000_omnibus.csv"))
  edge <- read.csv(file.path(nct, "21_NCT10000_edges.csv"))
  read <- function(n) read.csv(paper1_find(baseline$parent, n), check.names = FALSE)
  family <- read("17_family_selection_NCT_omnibus.csv")
  similarity <- read("17_family_selection_network_comparisons.csv")
  if (!all(c("repetition", "modality", "visit_A", "visit_B", "p_structure", "p_global", "permutations") %in% names(family)) ||
      !all(c("network", "repetition", "edge_correlation", "max_abs_difference") %in% names(similarity)))
    stop("Unexpected family-sensitivity table columns.")
  if (nrow(family) != 30L || !identical(as.integer(sort(unique(family$repetition))), 1:5) ||
      anyNA(family$permutations) || any(family$permutations != 1000L) ||
      anyDuplicated(family[c("modality", "visit_A", "visit_B", "repetition")]))
    stop("Expected six comparisons across five additional family selections.")
  if (nrow(similarity) != 30L || length(unique(similarity$network)) != 6L ||
      anyDuplicated(similarity[c("network", "repetition")]) ||
      any(!is.finite(similarity$edge_correlation)) || any(!is.finite(similarity$max_abs_difference)))
    stop("Family-sensitivity graph-comparison rows are incomplete.")
  summary <- list()
  for (m in c("FA", "T2w")) for (pair in combn(c("ses-02A", "ses-04A", "ses-06A"), 2, simplify = FALSE)) {
    z <- family[family$modality == m & family$visit_A == pair[1] & family$visit_B == pair[2], , drop = FALSE]
    if (nrow(z) != 5L) stop("Family sensitivity rows are incomplete.")
    for (test in c("p_structure", "p_global")) {
      p <- z[[test]]
      if (any(!is.finite(p))) stop("Missing family-sensitivity p-value.")
      summary[[length(summary) + 1L]] <- data.frame(modality = m, visit_A = pair[1], visit_B = pair[2],
        test = test, repetitions = 5L, raw_p_below_05 = sum(p < .05),
        minimum_raw_p = min(p), maximum_raw_p = max(p),
        threshold_pattern = if (all(p < .05)) "below .05 in all five selections" else
          if (all(p >= .05)) "not below .05 in any selection" else "crosses .05 across selections",
        permutations_per_selection = 1000L,
        interpretation = "descriptive sensitivity; overlapping samples, not independent replications")
    }
  }
  write.csv(do.call(rbind, summary), file.path(out, "22_family_selection_test_consistency.csv"), row.names = FALSE)
  ranges <- lapply(split(similarity, similarity$network), function(z) data.frame(network = z$network[1],
    edge_correlation_min = min(z$edge_correlation), edge_correlation_max = max(z$edge_correlation),
    maximum_absolute_edge_difference = max(z$max_abs_difference)))
  write.csv(do.call(rbind, ranges), file.path(out, "22_family_selection_network_similarity.csv"), row.names = FALSE)
  write.csv(omni, file.path(out, "22_primary_NCT10000_for_manuscript.csv"), row.names = FALSE)
  changed <- lapply(split(edge, edge$comparison), function(z) data.frame(comparison = z$comparison[1],
    edge_tests = nrow(z), BH_significant = sum(z$p_BH_within153 < .05),
    Holm_significant = sum(z$p_Holm_within153 < .05)))
  write.csv(do.call(rbind, changed), file.path(out, "22_edge_correction_summary.csv"), row.names = FALSE)
  writeLines(c(
    "Reviewer 1: six primary paired NCT comparisons were repeated with 10000 permutations; samples, adjustment, estimator and BH/Holm families were retained.",
    "Use the new omnibus and edge p-values together; do not select between old and new p-values by significance.",
    "Domain tests and family/other sensitivity tests retain their R1 1000-permutation settings; their resolution limitations still apply to any edge-level non-rejection.",
    "Reviewer 3: omnibus tests and domain-specific strength tests answer distinct questions and were treated as separate reporting families. Unadjusted omnibus p-values remain exploratory and do not control multiplicity across all omnibus claims.",
    "Reviewer 3: with one participant selected per family, the pooled family random intercept effectively represents participant-specific stable variation. Its prediction was not subtracted from the network inputs.",
    "Describe network similarity and p-value sensitivity separately. A threshold crossing is not proof that the underlying effect changed, and correlations near one do not guarantee identical inference.",
    "These files are local author-review outputs. Check DUC disclosure requirements before sharing any tables; no automatic public export is performed."
  ), file.path(out, "22_interpretation_notes.txt"))
  message("Reviewer summaries created: ", out)
  invisible(out)
}
