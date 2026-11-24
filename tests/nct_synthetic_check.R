# Optional API check using generated Gaussian numbers only; never ABCD records.
source("run_analysis.R")
paper1_missing_packages(c("NetworkComparisonTest", "bootnet", "qgraph"))
set.seed(3241)
A <- matrix(rnorm(160 * 18), 160, 18)
B <- .5 * A + matrix(rnorm(160 * 18), 160, 18)
colnames(A) <- colnames(B) <- paste0("node", 1:18)
K <- 19L
fit <- NetworkComparisonTest::NCT(A, B, gamma = .5, it = K,
  paired = TRUE, binary.data = FALSE, weighted = TRUE,
  test.edges = TRUE, edges = "all", progressbar = FALSE, p.adjust.methods = "none")
e <- paper1_extract_edge_p(fit, colnames(A))
stopifnot(nrow(e) == 153L, min(e$p_raw) >= 1 / (K + 1) - 1e-12,
  abs(fit$nwinv.pval - (sum(fit$nwinv.perm >= fit$nwinv.real) + 1) / (K + 1)) < 1e-12,
  abs(fit$glstrinv.pval - (sum(fit$glstrinv.perm >= fit$glstrinv.real) + 1) / (K + 1)) < 1e-12)
message("Synthetic NCT API check passed; no cohort data used.")
