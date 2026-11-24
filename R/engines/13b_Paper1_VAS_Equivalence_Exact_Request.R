# 13b_Paper1_VAS_Equivalence_Exact_Request.R
# Implements only the equivalence-analysis changes requested by VAS.
# No other analyses are added or changed.
#
# VAS requested:
# 1) TOST for each of 432 SMA–MRI pairs at alpha=.05, uncorrected,
#    using bounds +/- .05 and +/- .03 and 90% CIs.
# 2) Number and percentage equivalent per network and bound.
# 3) Network-level intersection–union claim:
#    supported only if all 72 pairs in that network pass at alpha=.05.
# 4) Largest |partial r| and its 90% CI per network.
# 5) Do not use Holm-across-432 for the new main equivalence conclusion.

options(stringsAsFactors = FALSE)

message("Choose: 13_unregularised_432_partial_correlations.csv")
p <- normalizePath(file.choose(), winslash="/", mustWork=TRUE)

if (basename(p) != "13_unregularised_432_partial_correlations.csv") {
  stop("Wrong file selected.")
}

x <- read.csv(p, check.names=FALSE, na.strings=c("NA","NaN",""))

required <- c(
  "node1","node2","network","N",
  "partial_r","fisher_z","SE_z"
)
missing <- setdiff(required, names(x))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse=", "))
}

if (nrow(x) != 432L) stop("Expected 432 SMA–MRI pairs.")
if (length(unique(x$network)) != 6L) stop("Expected six networks.")
if (any(table(x$network) != 72L)) stop("Expected 72 SMA–MRI pairs per network.")

# 90% Fisher-z CIs from the existing unregularised partial correlations and SEs.
zcrit <- qnorm(.95)
x$lower90 <- tanh(x$fisher_z - zcrit * x$SE_z)
x$upper90 <- tanh(x$fisher_z + zcrit * x$SE_z)

# Pairwise uncorrected TOST at alpha=.05 for the two requested bounds.
bounds <- c(.05, .03)
pairwise <- list()

for (bound in bounds) {
  e <- x
  delta <- atanh(bound)

  e$bound_r <- bound
  e$p_lower <- pnorm((e$fisher_z + delta) / e$SE_z, lower.tail=FALSE)
  e$p_upper <- pnorm((e$fisher_z - delta) / e$SE_z)
  e$p_TOST <- pmax(e$p_lower, e$p_upper)

  e$equivalent <- e$p_TOST <= .05

  # Consistency check with the requested 90% CI criterion.
  equiv_ci <- (e$lower90 >= -bound & e$upper90 <= bound)
  if (any(e$equivalent != equiv_ci)) {
    stop("TOST and 90% CI equivalence classifications do not match.")
  }

  pairwise[[as.character(bound)]] <- e
}

pairwise <- do.call(rbind, pairwise)
rownames(pairwise) <- NULL

# Number and percentage equivalent per network and bound,
# plus the requested intersection–union network-level claim.
summary_rows <- list()
k <- 0L

for (bound in bounds) {
  for (net in unique(x$network)) {
    e <- pairwise[
      pairwise$bound_r == bound & pairwise$network == net,
      , drop=FALSE
    ]

    k <- k + 1L
    summary_rows[[k]] <- data.frame(
      network = net,
      bound_r = bound,
      pairs = 72L,
      equivalent_n = sum(e$equivalent),
      equivalent_percent = 100 * mean(e$equivalent),
      all_72_pass_IUT = all(e$equivalent),
      stringsAsFactors = FALSE
    )
  }
}

network_summary <- do.call(rbind, summary_rows)
rownames(network_summary) <- NULL

# Largest |partial r| and corresponding 90% CI per network.
largest_rows <- list()
k <- 0L

for (net in unique(x$network)) {
  e <- x[x$network == net, , drop=FALSE]
  idx <- which.max(abs(e$partial_r))

  k <- k + 1L
  largest_rows[[k]] <- data.frame(
    network = net,
    largest_abs_partial_r = abs(e$partial_r[idx]),
    signed_partial_r = e$partial_r[idx],
    lower90 = e$lower90[idx],
    upper90 = e$upper90[idx],
    stringsAsFactors = FALSE
  )
}

largest_summary <- do.call(rbind, largest_rows)
rownames(largest_summary) <- NULL

# Save only the outputs requested by VAS.
outdir <- file.path(dirname(p), "VAS_EQUIVALENCE_EXACT_REQUEST")
dir.create(outdir, showWarnings=FALSE)

write.csv(
  pairwise,
  file.path(outdir, "VAS_pairwise_TOST_432.csv"),
  row.names=FALSE,
  na=""
)

write.csv(
  network_summary,
  file.path(outdir, "VAS_equivalence_counts_percent_IUT_by_network.csv"),
  row.names=FALSE,
  na=""
)

write.csv(
  largest_summary,
  file.path(outdir, "VAS_largest_partial_r_90CI_by_network.csv"),
  row.names=FALSE,
  na=""
)

message("")
message("COMPLETE")
message("Output folder: ", outdir)
message("")
message("Network equivalence summary:")
print(network_summary, row.names=FALSE)
message("")
message("Largest |partial r| and 90% CI:")
print(largest_summary, row.names=FALSE)
message("")
message("Reference to add/check in manuscript:")
message("Berger RL, Hsu JC. Bioequivalence trials, intersection-union tests and equivalence confidence sets. Statistical Science. 1996;11(4):283-319. doi:10.1214/ss/1032280304.")
