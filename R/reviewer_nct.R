paper1_extract_edge_p <- function(fit, nodes) {
  p <- fit$einv.pvals
  if (is.matrix(p) && identical(dim(p), c(18L, 18L))) {
    if ((!is.null(rownames(p)) && !identical(rownames(p), nodes)) ||
        (!is.null(colnames(p)) && !identical(colnames(p), nodes)))
      stop("NCT edge matrix has a different node order.")
    ij <- which(upper.tri(p), arr.ind = TRUE)
    e <- data.frame(i = ij[, 1], j = ij[, 2], p_raw = p[ij])
  } else {
    p <- as.data.frame(p, check.names = FALSE)
    if (!all(c("Var1", "Var2", "p-value") %in% names(p))) stop("Unsupported NCT edge table.")
    idx <- function(v) {
      s <- as.character(v); x <- match(s, nodes); bad <- is.na(x)
      x[bad] <- suppressWarnings(as.integer(s[bad]))
      if (anyNA(x) || any(x < 1L | x > 18L)) stop("Unknown NCT node labels.")
      x
    }
    a <- idx(p$Var1); b <- idx(p$Var2)
    e <- data.frame(i = pmin(a, b), j = pmax(a, b), p_raw = as.numeric(p[["p-value"]]))
  }
  if (nrow(e) != 153L || any(e$i >= e$j) || anyDuplicated(e[c("i", "j")]) ||
      any(!is.finite(e$p_raw)) || any(e$p_raw < 0 | e$p_raw > 1)) stop("Invalid edge-test output.")
  e[order(e$j, e$i), , drop = FALSE]
}

paper1_validate_primary <- function(nets, settings, helpers, enforce_sizes = TRUE) {
  keys <- unlist(lapply(names(helpers$MRI), function(m) paste(m, helpers$WAVES, sep = "_")))
  if (length(nets) != 6L || !setequal(names(nets), keys)) stop("Expected the six R1 networks.")
  if (!identical(settings$pooled_wave_fixed, FALSE) || settings$gamma != .5 || settings$seed != 20261003L)
    stop("These are not the R1 primary settings.")
  expected <- setNames(c(3295L, 4006L, 5348L, 3325L, 4017L, 5132L), keys)
  for (key in keys) {
    a <- nets[[key]]
    if (!is.matrix(a$x) || ncol(a$x) != 18L || nrow(a$x) != a$n ||
        a$n != length(a$ids) || anyNA(a$ids) || anyDuplicated(a$ids) ||
        anyNA(a$x) || any(!is.finite(a$x)) || !identical(colnames(a$x), helpers$NODES))
      stop("Invalid primary inputs: ", key)
    if (isTRUE(enforce_sizes) && a$n != expected[[key]]) stop("Sample differs from submitted R1: ", key)
    if (max(abs(helpers$network(a$x) - a$W)) > 1e-8)
      stop("Saved graph cannot be reconstructed with the installed software: ", key)
  }
  invisible(TRUE)
}

paper1_nct_revision <- function(cfg) {
  paper1_missing_packages(c("NetworkComparisonTest", "bootnet", "qgraph"))
  baseline <- paper1_baseline(cfg)
  # This input is participant-level. It is read locally, never exported to the repository.
  nf <- file.path(baseline$primary, "03_primary_networks_RESTRICTED.rds")
  sf <- file.path(baseline$primary, "03_settings.rds")
  of <- file.path(baseline$primary, "05_original_NCT_omnibus.csv")
  ef <- file.path(baseline$primary, "05_original_NCT_edges.csv")
  if (!all(file.exists(c(nf, sf, of, ef)))) stop("The completed R1 network inputs, settings and comparison tables are required.")
  settings <- readRDS(sf); nets <- readRDS(nf)
  h <- paper1_source_helpers(cfg$work_dir); h$CFG <- settings
  paper1_validate_primary(nets, settings, h, cfg$enforce_r1_sample_sizes)
  old <- read.csv(of, check.names = FALSE)
  old_edges <- read.csv(ef, check.names = FALSE)
  if (nrow(old) != 6L || anyDuplicated(old$comparison) || any(old$permutations != 1000L))
    stop("Expected six original 1000-permutation comparisons.")
  out <- file.path(cfg$work_dir, "R2_NCT_10000")
  private <- file.path(out, "RESTRICTED"); dir.create(private, recursive = TRUE, showWarnings = FALSE)
  K <- 10000L
  runtime_versions <- paper1_versions()
  signature <- list(schema = "primary-NCT-10000-v1", permutations = K,
    inputs = unname(tools::md5sum(c(nf, sf, of, ef))),
    R = paste(R.version$major, R.version$minor, sep = "."),
    versions = runtime_versions,
    code = unname(tools::md5sum(file.path(.paper1_repo, "R",
      c("reviewer_nct.R", "runtime.R", "engines/00_Paper1_Settings_and_Functions.R")))))
  sig <- file.path(private, "input_signature.rds")
  if (file.exists(sig)) {
    if (!identical(readRDS(sig), signature)) stop("Inputs, software or code changed. Use a new work_dir; do not mix checkpoints.")
  } else paper1_atomic_rds(signature, sig)
  completed_file <- file.path(private, "COMPLETE.rds")
  if (file.exists(completed_file)) {
    finished <- readRDS(completed_file)
    if (!identical(finished$signature, signature) ||
        !identical(unname(tools::md5sum(names(finished$checksums))), unname(finished$checksums)))
      stop("Completed R2 outputs changed. Use a new work directory.")
    message("All six R2 comparisons were already completed and verified: ", out)
    return(invisible(out))
  }
  write.csv(runtime_versions, file.path(private, "software_versions.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(private, "sessionInfo.txt"))
  omni <- edge_rows <- checks <- list()
  for (mod in names(h$MRI)) for (pair in combn(h$WAVES, 2, simplify = FALSE)) {
    a <- nets[[paste(mod, pair[1], sep = "_")]]; b <- nets[[paste(mod, pair[2], sep = "_")]]
    ids <- sort(intersect(a$ids, b$ids))
    if (length(ids) < settings$min_network_n) stop("Too few paired participants.")
    A <- a$x[match(ids, a$ids), , drop = FALSE]; B <- b$x[match(ids, b$ids), , drop = FALSE]
    label <- paste(mod, paste(pair, collapse = "_vs_"), sep = "_")
    seed <- as.integer(settings$seed + 500L + length(omni))
    previous <- old[old$comparison == label, , drop = FALSE]
    if (nrow(previous) != 1L || previous$paired_N != length(ids)) stop("R1 paired sample does not match: ", label)
    checkpoint <- file.path(private, paste0(label, "_NCT10000.rds"))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      if (!identical(saved$signature, signature) || saved$seed != seed || saved$comparison != label)
        stop("Checkpoint signature mismatch.")
      nt <- saved$fit
      message("Using verified completed comparison: ", label)
    } else {
      message("NCT 10000: ", label, "; paired N = ", length(ids))
      RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
      set.seed(seed)
      nt <- NetworkComparisonTest::NCT(A, B, gamma = .5, it = K,
        binary.data = FALSE, paired = TRUE, weighted = TRUE,
        test.edges = TRUE, edges = "all", progressbar = TRUE, p.adjust.methods = "none")
      paper1_atomic_rds(list(signature = signature, comparison = label, seed = seed, fit = nt), checkpoint)
    }
    if (length(nt$nwinv.perm) != K || length(nt$glstrinv.perm) != K ||
        any(!is.finite(nt$nwinv.perm)) || any(!is.finite(nt$glstrinv.perm)))
      stop("Wrong permutation count or invalid null distribution: ", label)
    row <- data.frame(comparison = label, modality = mod, visit_A = pair[1], visit_B = pair[2],
      paired_N = length(ids), permutations = K, seed = seed, M = as.numeric(nt$nwinv.real),
      p_structure = as.numeric(nt$nwinv.pval), global_A = as.numeric(nt$glstrinv.sep[1]),
      global_B = as.numeric(nt$glstrinv.sep[2]), p_global = as.numeric(nt$glstrinv.pval))
    delta <- max(abs(unlist(row[c("M", "global_A", "global_B")], use.names = FALSE) -
                     unlist(previous[c("M", "global_A", "global_B")], use.names = FALSE)))
    if (!is.finite(delta) || delta > 1e-8)
      stop("Observed statistics differ from R1. Check samples/software before interpreting: ", label)
    ps <- (sum(nt$nwinv.perm >= nt$nwinv.real) + 1) / (K + 1)
    pg <- (sum(nt$glstrinv.perm >= nt$glstrinv.real) + 1) / (K + 1)
    if (max(abs(c(ps - row$p_structure, pg - row$p_global))) > 1e-10)
      stop("Installed NCT p-value convention differs from the required plus-one estimator.")
    ep <- paper1_extract_edge_p(nt, h$NODES)
    EA <- h$edges(h$network(A)); EB <- h$edges(h$network(B))
    e <- merge(EA, ep, by = c("i", "j"), sort = FALSE)
    if (nrow(e) != 153L) stop("Edge join failed.")
    e$weight_A <- e$weight; e$weight <- NULL
    e$weight_B <- EB$weight[match(paste(e$i, e$j), paste(EB$i, EB$j))]
    e$comparison <- label; e$permutations <- K
    e$p_BH_within153 <- p.adjust(e$p_raw, "BH")
    e$p_Holm_within153 <- p.adjust(e$p_raw, "holm")
    if (any(e$p_raw < 1 / (K + 1) - 1e-12) ||
        any(abs(e$p_raw * (K + 1) - round(e$p_raw * (K + 1))) > 1e-7))
      stop("Edge p-values are inconsistent with a 10000-permutation plus-one test.")
    # Verify the observed edge estimates, not just the omnibus statistics.
    pe <- old_edges[old_edges$comparison == label, , drop = FALSE]
    k <- match(paste(e$i, e$j), paste(pe$i, pe$j))
    if (nrow(pe) != 153L || anyNA(k) ||
        max(abs(e$weight_A - pe$weight_A[k]), abs(e$weight_B - pe$weight_B[k])) > 1e-8)
      stop("Observed edge weights differ from the R1 paired networks.")
    row$p_structure_R1_1000 <- previous$p_structure
    row$p_global_R1_1000 <- previous$p_global
    omni[[label]] <- row; edge_rows[[label]] <- e
    checks[[label]] <- data.frame(comparison = label, paired_N = length(ids),
      observed_statistic_max_error = delta, edge_tests = nrow(e), permutations = K,
      smallest_raw_p = min(e$p_raw), BH_significant = sum(e$p_BH_within153 < .05),
      Holm_significant = sum(e$p_Holm_within153 < .05))
    # Checkpoints are private. Final tables are written only once all six tests validate.
  }
  om <- do.call(rbind, omni); ee <- do.call(rbind, edge_rows); ck <- do.call(rbind, checks)
  rownames(om) <- rownames(ee) <- rownames(ck) <- NULL
  stopifnot(nrow(om) == 6L, nrow(ee) == 918L, nrow(ck) == 6L)
  write.csv(om, file.path(out, "21_NCT10000_omnibus.csv"), row.names = FALSE)
  write.csv(ee, file.path(out, "21_NCT10000_edges.csv"), row.names = FALSE)
  write.csv(ck, file.path(out, "21_NCT10000_validation.csv"), row.names = FALSE)
  limits <- data.frame(permutations = c(1000L, 10000L), edges = 153L,
    minimum_raw_p = 1 / (c(1000, 10000) + 1),
    minimum_Holm_adjusted_p = 153 / (c(1000, 10000) + 1),
    minimum_possible_BH_rejections_at_05 = ceiling(153 / (.05 * (c(1000, 10000) + 1))))
  write.csv(limits, file.path(out, "21_permutation_resolution.csv"), row.names = FALSE)
  paper1_atomic_rds(list(signature = signature,
    checksums = tools::md5sum(file.path(out, c("21_NCT10000_omnibus.csv", "21_NCT10000_edges.csv", "21_NCT10000_validation.csv", "21_permutation_resolution.csv")))),
    file.path(private, "COMPLETE.rds"))
  message("Six comparisons completed and checked. Outputs: ", out)
  message("All outputs remain local. Run reviewer_summaries next; do not upload the work directory.")
  invisible(out)
}
