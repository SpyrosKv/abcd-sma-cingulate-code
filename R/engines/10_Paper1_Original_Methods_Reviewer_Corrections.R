# Paper 1: original methods + requested reviewer/VAS corrections.
# Uses the existing v1 functions and audited merge; keeps v1 outputs intact.
# No wave fixed effect; original bootnet and NCT; no equivalence tests.

required <- c(
  "haven", "dplyr", "tidyselect", "lme4", "qgraph", "polycor",
  "mice", "bootnet", "NetworkComparisonTest"
)

missing_pkgs <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) install.packages(missing_pkgs)

if (any(!vapply(
  required, requireNamespace, logical(1), quietly = TRUE
))) {
  stop("Required packages could not be loaded.")
}

message(
  "Select the OLD 00_Paper1_Settings_and_Functions.R ",
  "from your v1 folder."
)

old_settings <- normalizePath(
  file.choose(), winslash = "/", mustWork = TRUE
)

if (basename(old_settings) != "00_Paper1_Settings_and_Functions.R") {
  stop("Please select the existing v1 00 settings script.")
}

old_dir <- dirname(old_settings)
CODE_DIR <- old_dir

source(old_settings, local = TRUE)

old_out <- OUT

needed <- c(
  "01_merged_RESTRICTED.rds",
  "01_sma_categories.rds"
)

if (!all(file.exists(file.path(old_out, needed)))) {
  stop(
    "The audited 01 outputs must be in the old ",
    "Paper1_revision_outputs folder."
  )
}

# Restore the original analysis settings.
CFG$pooled_wave_fixed <- FALSE
CFG$permutations <- 1000L
CFG$bootstrap_B <- 1000L
CFG$case_B <- 1000L
CFG$mice_m <- 20L
CFG$domain_permutation <- FALSE
CFG$equivalence_bound <- NA_real_

# Motion retained under Kanwal's explicit later instruction.
# This does NOT certify correspondence with FA source scans.

OUT <- file.path(
  old_dir,
  "Paper1_revision_outputs_original_methods_v2"
)

if (dir.exists(OUT)) {
  stop(
    "v2 output folder already exists. Rename it before another run; ",
    "do not mix runs."
  )
}

dir.create(OUT, recursive = TRUE)

if (!all(file.copy(file.path(old_out, needed), OUT))) {
  stop("Could not copy 01 inputs.")
}

# Copy unchanged sample/QC/selection audits only.
# Do not copy old network results or old 02d/02f results.

audit_files <- list.files(
  old_out,
  pattern = "^(01_|02_|02b_|02e_).*[.](csv|txt)$",
  full.names = TRUE
)

if (length(audit_files) &&
    !all(file.copy(audit_files, OUT))) {
  stop("Audit copy failed.")
}

# Run existing stages without reloading the old configuration.

run_existing <- function(filename) {
  p <- file.path(old_dir, filename)

  if (!file.exists(p)) stop("Missing script: ", p)

  txt <- readLines(p, warn = FALSE)

  if (!grepl(
    "^source\\('00_Paper1_Settings_and_Functions.R'\\)",
    txt[1]
  )) {
    stop("Unexpected script header: ", filename)
  }

  eval(
    parse(text = txt[-1]),
    envir = environment(run_existing)
  )
}


# ============================================================
# 03: PRIMARY RESIDUALISATION AND NETWORKS
# ============================================================
# Same nodes, Pearson correlations and EBICglasso.
# VAS corrections: retain family contribution, one per family.
# No wave fixed effect.

run_existing(
  "03_Paper1_Primary_Residuals_and_Networks.R"
)

dg <- read.csv(
  file.path(OUT, "03_residualisation_diagnostics.csv")
)

if (any(grepl(
  "visit_factor|age.*:.*visit|visit.*:.*age",
  dg$formula
))) {
  stop("Unexpected wave term in residualisation.")
}

if (any(dg$extraction != "retain_family")) {
  stop("Unexpected residual extraction.")
}


# ============================================================
# 04: ORIGINAL BOOTNET BOOTSTRAP
# ============================================================

nets <- read_stage(
  "03_primary_networks_RESTRICTED.rds"
)

cs_rows <- edge_rows <- list()

for (key in names(nets)) {

  message("04 bootnet: ", key)

  x <- as.data.frame(nets[[key]]$x)

  fit <- bootnet::estimateNetwork(
    x,
    default = "EBICglasso",
    corMethod = "cor",
    tuning = CFG$gamma
  )

  if (max(abs(fit$graph - nets[[key]]$W)) > 1e-6) {
    stop(
      "bootnet and primary EBICglasso graphs differ: ",
      key
    )
  }

  set.seed(
    CFG$seed + 400L + match(key, names(nets))
  )

  be <- bootnet::bootnet(
    fit,
    nBoot = 1000L,
    type = "nonparametric",
    statistics = c(
      "edge", "strength", "expectedInfluence"
    ),
    nCores = 1L
  )

  bc <- bootnet::bootnet(
    fit,
    nBoot = 1000L,
    type = "case",
    statistics = c(
      "strength", "expectedInfluence", "edge",
      "closeness", "betweenness"
    ),
    nCores = 1L
  )

  saveRDS(
    be,
    file.path(
      OUT,
      paste0("04_", key, "_edge_bootnet.rds")
    )
  )

  saveRDS(
    bc,
    file.path(
      OUT,
      paste0("04_", key, "_case_bootnet.rds")
    )
  )

  pdf(
    file.path(
      OUT,
      paste0("04_", key, "_edge_accuracy.pdf")
    ),
    width = 10,
    height = 12
  )

  print(plot(be, labels = FALSE, order = "sample"))
  dev.off()

  pdf(
    file.path(
      OUT,
      paste0("04_", key, "_case_stability.pdf")
    ),
    width = 10,
    height = 8
  )

  print(plot(
    bc,
    statistics = c("strength", "expectedInfluence")
  ))

  dev.off()

  cs <- bootnet::corStability(
    bc,
    statistics = c("strength", "expectedInfluence")
  )

  cs_rows[[key]] <- data.frame(
    network = key,
    statistic = names(cs),
    CS = as.numeric(cs)
  )

  bt <- be$bootTable

  req <- c("type", "node1", "node2", "value")

  if (!all(req %in% names(bt))) {
    stop("Unexpected bootnet table structure.")
  }

  bt <- bt[
    bt$type == "edge",
    req,
    drop = FALSE
  ]

  ei <- bt |>
    dplyr::group_by(node1, node2) |>
    dplyr::summarise(
      lower = quantile(value, .025, na.rm = TRUE),
      upper = quantile(value, .975, na.rm = TRUE),
      successful_values = sum(is.finite(value)),
      .groups = "drop"
    )

  ei$network <- key

  ei$scope <- paste(
    "pointwise regularised-edge bootstrap;",
    "conditional on preprocessing;",
    "not a significance test"
  )

  edge_rows[[key]] <- ei
}

csv(
  dplyr::bind_rows(cs_rows),
  "04_bootnet_centrality_CS.csv"
)

csv(
  dplyr::bind_rows(edge_rows),
  "04_bootnet_edge_intervals.csv"
)

save_environment("04")


# ============================================================
# 05: ORIGINAL PAIRED NCT, 1000 PERMUTATIONS
# ============================================================
# BH/Holm within each comparison's 153 edges.

extract_nct_edges <- function(fit) {

  p <- fit$einv.pvals

  if (is.matrix(p) &&
      identical(dim(p), c(18L, 18L))) {

    ij <- which(upper.tri(p), arr.ind = TRUE)

    ans <- data.frame(
      i = ij[, 1],
      j = ij[, 2],
      p_raw = p[ij]
    )

  } else if (is.data.frame(p) || is.matrix(p)) {

    p <- as.data.frame(p, check.names = FALSE)

    if (!all(
      c("Var1", "Var2", "p-value") %in% names(p)
    )) {
      stop(
        "Unexpected NCT edge table; ",
        "send names(fit$einv.pvals)."
      )
    }

    index <- function(v) {

      s <- as.character(v)
      i <- match(s, NODES)

      missing <- is.na(i)

      i[missing] <- suppressWarnings(
        as.integer(s[missing])
      )

      if (anyNA(i) || any(i < 1 | i > 18)) {
        stop("Unrecognised NCT node names.")
      }

      i
    }

    a <- index(p$Var1)
    b <- index(p$Var2)

    ans <- data.frame(
      i = pmin(a, b),
      j = pmax(a, b),
      p_raw = as.numeric(p[["p-value"]])
    )

  } else {

    stop(
      "Unsupported NCT edge output; ",
      "do not infer its ordering."
    )
  }

  if (
    nrow(ans) != 153L ||
    anyDuplicated(ans[c("i", "j")]) ||
    any(!is.finite(ans$p_raw)) ||
    any(ans$p_raw < 0 | ans$p_raw > 1)
  ) {
    stop("Invalid NCT edge results.")
  }

  ans
}

omnibus <- edge_tests <- list()

for (mod in names(MRI)) {

  for (pair in combn(WAVES, 2, simplify = FALSE)) {

    a <- nets[[paste(mod, pair[1], sep = "_")]]
    b <- nets[[paste(mod, pair[2], sep = "_")]]

    ids <- sort(intersect(a$ids, b$ids))

    if (length(ids) < CFG$min_network_n) {
      stop("Too few paired participants.")
    }

    ia <- match(ids, a$ids)
    ib <- match(ids, b$ids)

    stopifnot(identical(a$ids[ia], b$ids[ib]))

    A <- a$x[ia, , drop = FALSE]
    B <- b$x[ib, , drop = FALSE]

    label <- paste(
      mod,
      paste(pair, collapse = "_vs_"),
      sep = "_"
    )

    message(
      "05 original NCT: ",
      label,
      "; paired N=",
      length(ids)
    )

    set.seed(
      CFG$seed + 500L + length(omnibus)
    )

    nt <- NetworkComparisonTest::NCT(
      A,
      B,
      gamma = .5,
      it = 1000L,
      binary.data = FALSE,
      paired = TRUE,
      weighted = TRUE,
      test.edges = TRUE,
      edges = "all",
      progressbar = TRUE,
      p.adjust.methods = "none"
    )

    saveRDS(
      nt,
      file.path(
        OUT,
        paste0("05_", label, "_NCT.rds")
      )
    )

    omnibus[[label]] <- data.frame(
      comparison = label,
      modality = mod,
      visit_A = pair[1],
      visit_B = pair[2],
      paired_N = length(ids),
      permutations = 1000L,
      M = as.numeric(nt$nwinv.real),
      p_structure = as.numeric(nt$nwinv.pval),
      global_A = as.numeric(nt$glstrinv.sep[1]),
      global_B = as.numeric(nt$glstrinv.sep[2]),
      p_global = as.numeric(nt$glstrinv.pval)
    )

    ep <- extract_nct_edges(nt)

    e <- merge(
      edges(network(A)),
      ep,
      by = c("i", "j"),
      sort = FALSE
    )

    if (nrow(e) != 153L) {
      stop("NCT edge matching failed.")
    }

    eb <- edges(network(B))

    e$weight_A <- e$weight

    e$weight_B <- eb$weight[
      match(
        paste(e$i, e$j),
        paste(eb$i, eb$j)
      )
    ]

    e$weight <- NULL
    e$comparison <- label

    e$p_BH_within153 <- p.adjust(
      e$p_raw, "BH"
    )

    e$p_Holm_within153 <- p.adjust(
      e$p_raw, "holm"
    )

    edge_tests[[label]] <- e
  }
}

csv(
  dplyr::bind_rows(omnibus),
  "05_original_NCT_omnibus.csv"
)

csv(
  dplyr::bind_rows(edge_tests),
  "05_original_NCT_edges.csv"
)

save_environment("05")


# ============================================================
# 06: FAMILY AND ORDINAL SENSITIVITY
# ============================================================
# Reviewer-requested sensitivities; no wave fixed effect.

run_existing(
  "06_Paper1_Family_Residualisation_Motion_and_Ordinal_Sensitivity.R"
)


# ============================================================
# 07: ORIGINAL RESIDUAL-SCALE PMM SENSITIVITY
# ============================================================
# Original m=20.
# FA Year 2: maxit=50; other networks: maxit=20.
# Original seed=500.
#
# Structural missingness excluded BEFORE PMM.
# QC failures never enter PMM.
#
# AND/OR refer to observed availability among eligible,
# QC-passing records.
#
# Original complete(imp,1) graph retained as DESCRIPTIVE.
# All 20 graphs exported; no Rubin-pooled inference added.

d <- read_stage("01_merged_RESTRICTED.rds")

res <- read_stage(
  "03_primary_residuals_RESTRICTED.rds"
)

mi_audit <- mi_edges <- mi_comparison <- mi_logged <- list()

for (mod in names(MRI)) {

  z <- res[[mod]]$source
  r <- res[[mod]]$residuals

  stopifnot(
    all(z[[GATE]] == 1),
    all(z[[QC[[mod]]]] == 1),
    all(z$selected_family)
  )

  for (w in WAVES) {

    for (rule in c("AND", "OR")) {

      has_sma <- rowSums(
        !is.na(r[, 1:6, drop = FALSE])
      ) > 0

      has_mri <- rowSums(
        !is.na(r[, 7:18, drop = FALSE])
      ) > 0

      eligible_rows <- if (rule == "AND") {
        has_sma & has_mri
      } else {
        has_sma | has_mri
      }

      take <- z$visit == w & eligible_rows

      x <- as.data.frame(
        r[take, , drop = FALSE]
      )

      label <- paste(
        mod, w, rule, sep = "_"
      )

      if (
        nrow(x) < CFG$min_network_n ||
        any(colSums(!is.na(x)) < 10L)
      ) {
        stop(
          "Insufficient observations/donors: ",
          label
        )
      }

      iterations <- if (
        mod == "FA" && w == "ses-02A"
      ) {
        50L
      } else {
        20L
      }

      message(
        "07 original residual PMM: ",
        label
      )

      imp <- if (anyNA(x)) {

        mice::mice(
          x,
          m = 20L,
          method = "pmm",
          maxit = iterations,
          seed = 500L,
          printFlag = TRUE
        )

      } else {
        NULL
      }

      saveRDS(
        list(
          imputation = imp,
          complete_data = if (is.null(imp)) x else NULL
        ),
        file.path(
          OUT,
          paste0(
            "07_", label, "_PMM_RESTRICTED.rds"
          )
        )
      )

      events <- if (is.null(imp)) {
        NULL
      } else {
        imp$loggedEvents
      }

      if (!is.null(events)) {
        events$network_rule <- label
        mi_logged[[label]] <- events
      }

      if (anyNA(x)) {

        pdf(
          file.path(
            OUT,
            paste0("07_", label, "_trace.pdf")
          )
        )

        print(plot(imp))
        dev.off()
      }

      key <- paste(mod, w, sep = "_")

      for (k in 1:20) {

        completed <- if (is.null(imp)) {
          x
        } else {
          mice::complete(imp, k)
        }

        if (
          anyNA(completed) ||
          any(!is.finite(as.matrix(completed)))
        ) {
          stop(
            "Invalid completed residual dataset: ",
            label
          )
        }

        W <- network(completed)

        ee <- edges(W)

        ee$network <- key
        ee$rule <- rule
        ee$imputation <- k
        ee$N <- nrow(completed)

        mi_edges[[
          paste(label, k, sep = "_")
        ]] <- ee

        co <- compare_networks(
          nets[[key]]$W,
          W,
          paste(label, k, sep = "_")
        )

        co$N_primary <- nets[[key]]$n
        co$N_MI <- nrow(completed)

        co$scope <- paste(
          "descriptive sensitivity;",
          "populations may differ"
        )

        mi_comparison[[
          paste(label, k, sep = "_")
        ]] <- co

        if (k == 1L) {

          saveRDS(
            W,
            file.path(
              OUT,
              paste0(
                "07_", label,
                "_legacy_first_imputation_graph.rds"
              )
            )
          )
        }
      }

      mi_audit[[label]] <- data.frame(
        modality = mod,
        visit = w,
        rule = rule,
        N = nrow(x),
        complete_18 = sum(complete.cases(x)),
        missing_cells = sum(is.na(x)),
        m = 20L,
        iterations = iterations,
        imputed_scale = "adjusted continuous values",
        logged_events = if (is.null(events)) {
          0L
        } else {
          nrow(events)
        }
      )
    }
  }
}

csv(
  dplyr::bind_rows(mi_audit),
  "07_original_PMM_eligibility_audit.csv"
)

csv(
  dplyr::bind_rows(mi_edges),
  "07_original_PMM_all20_edges_DESCRIPTIVE.csv"
)

csv(
  dplyr::bind_rows(mi_comparison),
  "07_original_PMM_vs_primary_DESCRIPTIVE.csv"
)

if (length(mi_logged)) {
  csv(
    dplyr::bind_rows(mi_logged),
    "07_original_PMM_logged_events.csv"
  )
}

save_environment("07")


# All results remain in the authorised work directory. No public export.
save_environment("08")
message("R1 rerun complete. All outputs remain restricted until reviewed.")
