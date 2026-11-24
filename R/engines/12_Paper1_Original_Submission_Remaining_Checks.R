# Remaining original-submission diagnostics and numerical VAS verification.
# Does not change the primary sample, covariates or estimated networks.
# Uses the existing 00 functions and saved v2 inputs/models.

message(
  "Select the original 00_Paper1_Settings_and_Functions.R ",
  "from the v1 folder."
)

p00 <- normalizePath(
  file.choose(), winslash = "/", mustWork = TRUE
)

if (basename(p00) != "00_Paper1_Settings_and_Functions.R") {
  stop("Wrong settings file.")
}

CODE_DIR <- dirname(p00)
source(p00, local = TRUE)

V2 <- file.path(
  CODE_DIR,
  "Paper1_revision_outputs_original_methods_v2"
)

required <- c(
  "01_merged_RESTRICTED.rds",
  "03_primary_networks_RESTRICTED.rds",
  "03_primary_residuals_RESTRICTED.rds",
  "03_settings.rds",
  "05_original_NCT_omnibus.csv"
)

if (!all(file.exists(file.path(V2, required)))) {
  stop("Completed v2 outputs not found.")
}

CFG <- readRDS(file.path(V2, "03_settings.rds"))

if (!identical(CFG$pooled_wave_fixed, FALSE) ||
    CFG$gamma != .5) {
  stop("Unexpected v2 configuration.")
}

CFG$permutations <- 1000L

OUT <- file.path(
  CODE_DIR,
  "Paper1_original_required_checks"
)

if (dir.exists(OUT)) {
  stop(
    "Output folder already exists; ",
    "rename it before a fresh run."
  )
}

dir.create(OUT)

d <- readRDS(
  file.path(V2, "01_merged_RESTRICTED.rds")
)

nets <- readRDS(
  file.path(V2, "03_primary_networks_RESTRICTED.rds")
)

res <- readRDS(
  file.path(V2, "03_primary_residuals_RESTRICTED.rds")
)

expected <- as.vector(
  outer(names(MRI), WAVES, paste, sep = "_")
)

if (!setequal(names(nets), expected)) {
  stop("Unexpected primary networks.")
}

check_keys(d, "Merged data")

selected <- unique(
  d[d$selected_family, c("participant_id", "family")]
)

if (anyNA(selected$family) ||
    anyDuplicated(selected$family)) {
  stop("More than one selected participant per family.")
}


# ============================================================
# A. NUMERICAL VERIFICATION OF VAS RESIDUALISATION
# ============================================================

audit <- list()

for (mod in names(MRI)) {

  z <- res[[mod]]$source
  vars <- c(SMA, MRI[[mod]])

  if (!all(flags_for(z, mod)$base)) {
    stop("Invalid primary adjustment population.")
  }

  for (j in seq_along(vars)) {

    v <- vars[j]

    fit <- res[[mod]]$models[[
      paste("pooled", v, sep = "__")
    ]]

    if (is.null(fit)) {
      stop("Saved model missing: ", v)
    }

    if (any(grepl(
      "visit|wave",
      names(lme4::fixef(fit)),
      ignore.case = TRUE
    ))) {
      stop("Unexpected visit/wave fixed effect.")
    }

    used <- !is.na(z[[v]])

    zz <- droplevels(
      z[used, , drop = FALSE]
    )

    zz$motion <- zz[[CFG$motion_var]]

    groups <- names(lme4::ranef(fit))

    if (!"family" %in% groups) {
      stop("Pooled family intercept missing.")
    }

    remove_re <- if ("scanner" %in% groups) {
      ~(1 | scanner)
    } else {
      NA
    }

    expected_r <- zz[[v]] - predict(
      fit,
      newdata = zz,
      re.form = remove_re
    )

    saved_r <- res[[mod]]$residuals[used, j]

    err <- max(abs(expected_r - saved_r))

    if (!is.finite(err) || err > 1e-10) {
      stop("Residual reconstruction failed: ", v)
    }

    conditional <- zz[[v]] - predict(
      fit,
      newdata = zz,
      re.form = NULL
    )

    audit[[length(audit) + 1]] <- data.frame(
      modality = mod,
      variable = v,
      N = sum(used),
      maximum_reconstruction_error = err,
      retained_family = TRUE,
      extracted_SD = sd(saved_r),
      conditional_SD = sd(conditional)
    )
  }

  for (w in WAVES) {

    key <- paste(mod, w, sep = "_")
    n <- nets[[key]]

    rows <- match(
      paste(n$ids, w),
      paste(d$participant_id, d$visit)
    )

    if (anyNA(rows) || anyDuplicated(n$ids)) {
      stop("Primary ID matching failed.")
    }

    zz <- d[rows, , drop = FALSE]

    if (!all(flags_for(zz, mod)$included)) {
      stop("Primary inclusion failed: ", key)
    }

    if (
      n$n != length(n$ids) ||
      nrow(n$x) != n$n ||
      anyDuplicated(zz$family)
    ) {
      stop("Primary sample integrity failed: ", key)
    }

    if (max(abs(network(n$x) - n$W)) > 1e-10) {
      stop("Primary graph reconstruction failed.")
    }
  }
}

csv(
  dplyr::bind_rows(audit),
  "12_primary_residual_reconstruction.csv"
)


# ============================================================
# B. ORIGINAL WAVE-SPECIFIC SENSITIVITY
# ============================================================
# Same selected families, covariates and residual extraction.
#
# With one observation per family within each wave, the family
# random-intercept variance is not separately identifiable.
# residualise() therefore omits that term automatically.

wave_comparison <- wave_diag <- wave_edges <- list()

for (mod in names(MRI)) {

  z <- res[[mod]]$source

  for (w in WAVES) {

    zz <- z[z$visit == w, , drop = FALSE]

    if (anyDuplicated(zz$family)) {
      stop("Unexpected repeated family within a wave.")
    }
  }

  rr <- residualise(
    z,
    mod,
    scope = "wave",
    family_model = TRUE,
    extraction = "retain_family"
  )

  wave_diag[[mod]] <- rr$diagnostics

  for (w in WAVES) {

    key <- paste(mod, w, sep = "_")
    p <- nets[[key]]

    take <- z$visit == w & complete.cases(rr$data)
    ids <- z$participant_id[take]

    if (!setequal(ids, p$ids)) {
      stop("Wave sensitivity changed the sample: ", key)
    }

    x <- rr$data[take, , drop = FALSE][
      match(p$ids, ids),
      ,
      drop = FALSE
    ]

    W <- network(x)

    co <- compare_networks(
      p$W,
      W,
      paste("pooled_vs_wave_specific", key)
    )

    co$N <- p$n
    co$scope <- paste(
      "same participants; same covariates;",
      "separate-wave adjustment"
    )

    wave_comparison[[key]] <- co

    e <- edges(W)
    e$network <- key
    wave_edges[[key]] <- e
  }
}

csv(
  dplyr::bind_rows(wave_diag),
  "12_wave_specific_model_diagnostics.csv"
)

csv(
  dplyr::bind_rows(wave_comparison),
  "12_wave_specific_comparisons.csv"
)

csv(
  dplyr::bind_rows(wave_edges),
  "12_wave_specific_edges.csv"
)


# ============================================================
# C. ORIGINAL MISSINGNESS DIAGNOSTICS AND PMM VERIFICATION
# ============================================================
# Structural non-applicability never enters MCAR or PMM.
# Little's test uses the original AND network definition,
# now restricted to eligible, QC-passing participants.

if (!requireNamespace("naniar", quietly = TRUE)) {
  install.packages("naniar")
}

need("naniar")

missing_tests <- mi_integrity <- list()

for (mod in names(MRI)) {

  for (w in WAVES) {

    z <- res[[mod]]$source
    r <- res[[mod]]$residuals

    has_sma <- rowSums(
      !is.na(r[, 1:6, drop = FALSE])
    ) > 0

    has_mri <- rowSums(
      !is.na(r[, 7:18, drop = FALSE])
    ) > 0

    for (rule in c("AND", "OR")) {

      rule_mask <- if (rule == "AND") {
        has_sma & has_mri
      } else {
        has_sma | has_mri
      }

      use <- z$visit == w & rule_mask

      x <- as.data.frame(
        r[use, , drop = FALSE]
      )

      label <- paste(
        mod, w, rule, sep = "_"
      )

      file <- file.path(
        V2,
        paste0(
          "07_", label, "_PMM_RESTRICTED.rds"
        )
      )

      if (!file.exists(file)) {
        stop("Saved PMM object missing: ", label)
      }

      obj <- readRDS(file)

      imp <- if (inherits(obj, "mids")) {
        obj
      } else {
        obj$imputation
      }

      saved <- if (is.null(imp)) {
        obj$complete_data
      } else {
        imp$data
      }

      if (
        is.null(saved) ||
        !isTRUE(all.equal(
          as.matrix(x),
          as.matrix(saved),
          check.attributes = FALSE,
          tolerance = 1e-10
        ))
      ) {
        stop(
          "PMM input differs from eligible adjusted sample: ",
          label
        )
      }

      expected_iter <- if (
        mod == "FA" && w == "ses-02A"
      ) {
        50L
      } else {
        20L
      }

      if (!is.null(imp) &&
          (imp$m != 20L ||
           imp$iteration != expected_iter)) {
        stop("Unexpected PMM configuration: ", label)
      }

      if (
        !all(z[[GATE]][use] == 1) ||
        !all(z[[QC[[mod]]]][use] == 1)
      ) {
        stop("Non-eligible/QC-failing record in PMM.")
      }

      mi_integrity[[label]] <- data.frame(
        network = label,
        N = nrow(x),
        missing_cells = sum(is.na(x)),
        eligible_only = TRUE,
        QC_pass_only = TRUE,
        stored_input_matches = TRUE,
        m = if (is.null(imp)) 0L else imp$m
      )

      if (rule != "AND") next

      pattern <- mice::md.pattern(
        x, plot = FALSE
      )

      write.csv(
        pattern,
        file.path(
          OUT,
          paste0(
            "12_", label, "_missing_patterns.csv"
          )
        )
      )

      flux <- mice::flux(x)
      flux$variable <- rownames(flux)

      csv(
        flux,
        paste0(
          "12_", label, "_missing_flux.csv"
        )
      )

      warnings <- character()

      test <- if (anyNA(x)) {

        tryCatch(
          withCallingHandlers(
            naniar::mcar_test(x),
            warning = function(e) {
              warnings <<- c(
                warnings,
                conditionMessage(e)
              )
              invokeRestart("muffleWarning")
            }
          ),
          error = function(e) e
        )

      } else {
        NULL
      }

      ok <- !is.null(test) &&
        !inherits(test, "error")

      missing_tests[[label]] <- data.frame(
        network = label,
        N = nrow(x),
        missing_cells = sum(is.na(x)),

        status = if (is.null(test)) {
          "no_missing_values"
        } else if (ok) {
          "completed"
        } else {
          "failed"
        },

        statistic = if (ok) {
          as.numeric(test$statistic[1])
        } else {
          NA_real_
        },

        df = if (ok) {
          as.numeric(test$df[1])
        } else {
          NA_real_
        },

        p = if (ok) {
          as.numeric(test$p.value[1])
        } else {
          NA_real_
        },

        warnings = paste(
          unique(warnings),
          collapse = "; "
        ),

        error = if (inherits(test, "error")) {
          conditionMessage(test)
        } else {
          ""
        }
      )
    }
  }
}

csv(
  dplyr::bind_rows(mi_integrity),
  "12_PMM_input_integrity.csv"
)

csv(
  dplyr::bind_rows(missing_tests),
  "12_Little_MCAR_eligible_AND.csv"
)


# ============================================================
# D. ORIGINAL POST-HOC REGRESSION SENSITIVITY
# ============================================================
# Original pipeline: top 5 SMA–MRI edges per network.
# Selected after network estimation; descriptive sensitivity.
# These coefficients are NOT network partial correlations.
#
# Family variance cannot be estimated within a wave with
# one observation per family.
# Scanner random intercept retained where estimable.
# FA motion retained under Kanwal's later instruction.

regression <- list()

for (mod in names(MRI)) {

  for (w in WAVES) {

    key <- paste(mod, w, sep = "_")
    p <- nets[[key]]

    ee <- edges(p$W)

    ee <- ee[
      ee$domain == "SMA-MRI",
      ,
      drop = FALSE
    ]

    ee <- ee[
      order(-abs(ee$weight), ee$i, ee$j),
      ,
      drop = FALSE
    ]

    ee <- head(ee, 5)

    for (k in seq_len(nrow(ee))) {

      sma <- SMA[ee$i[k]]
      mri <- MRI[[mod]][ee$j[k] - 6L]

      z <- res[[mod]]$source

      z <- droplevels(
        z[
          z$visit == w &
            !is.na(z[[sma]]) &
            !is.na(z[[mri]]),
          ,
          drop = FALSE
        ]
      )

      z$y <- z[[mri]]
      z$sma <- z[[sma]]
      z$motion <- z[[CFG$motion_var]]

      if (anyDuplicated(z$family)) {
        stop("Unexpected family repetition in regression.")
      }

      if (nrow(z) < 100L || sd(z$sma) == 0) {
        stop("Insufficient regression data: ", key)
      }

      fixed <- c("scale(sma)", "age")

      if (nlevels(z$software) > 1) {
        fixed <- c(fixed, "software")
      }

      if (mod == "FA") {
        fixed <- c(fixed, "motion")
      }

      random <- if (
        nlevels(z$scanner) > 1 &&
        nlevels(z$scanner) < nrow(z)
      ) {
        "(1|scanner)"
      } else {
        character()
      }

      formula <- as.formula(
        paste(
          "y ~",
          paste(
            c(fixed, random),
            collapse = " + "
          )
        )
      )

      warnings <- character()

      fit <- withCallingHandlers(

        if (length(random)) {

          lme4::lmer(
            formula,
            z,
            REML = TRUE,
            control = lme4::lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 200000)
            )
          )

        } else {
          lm(formula, z)
        },

        warning = function(e) {
          warnings <<- c(
            warnings,
            conditionMessage(e)
          )
          invokeRestart("muffleWarning")
        }
      )

      beta <- if (inherits(fit, "merMod")) {
        lme4::fixef(fit)["scale(sma)"]
      } else {
        coef(fit)["scale(sma)"]
      }

      ci <- suppressMessages(
        confint(
          fit,
          parm = "scale(sma)",
          method = "Wald"
        )
      )

      conv <- if (inherits(fit, "merMod")) {
        fit@optinfo$conv$lme4$messages
      } else {
        NULL
      }

      regression[[length(regression) + 1]] <- data.frame(
        network = key,
        rank = k,
        SMA_node = ee$node1[k],
        MRI_node = ee$node2[k],
        network_weight = ee$weight[k],
        N = nrow(z),
        formula = paste(
          deparse(formula),
          collapse = " "
        ),
        beta = as.numeric(beta),
        lower = as.numeric(ci[1]),
        upper = as.numeric(ci[2]),
        direction_matches =
          sign(beta) == sign(ee$weight[k]),
        converged = !length(conv),

        singular = if (inherits(fit, "merMod")) {
          lme4::isSingular(fit)
        } else {
          FALSE
        },

        family_term =
          "not_identifiable_one_observation_per_family",

        warnings = paste(
          c(unique(warnings), conv),
          collapse = "; "
        ),

        interpretation = paste(
          "original selected-edge descriptive sensitivity;",
          "not a network-edge significance test"
        )
      )
    }
  }
}

csv(
  dplyr::bind_rows(regression),
  "12_ORIGINAL_posthoc_regression_sensitivity.csv"
)


# ============================================================
# E. VALIDATE ALREADY COMPLETED REQUESTED ANALYSES
# ============================================================

checks <- list()

for (name in c(
  "04_bootnet_centrality_CS.csv",
  "04_bootnet_edge_intervals.csv",
  "05_original_NCT_omnibus.csv",
  "05_original_NCT_edges.csv",
  "06_ordinal_matrix_diagnostics.csv",
  "02b_model_sample_and_convergence.csv"
)) {

  file <- file.path(V2, name)

  if (!file.exists(file)) {
    stop("Expected result missing: ", name)
  }

  a <- read.csv(file)

  checks[[name]] <- data.frame(
    file = name,
    rows = nrow(a),
    exists = TRUE
  )
}

boot <- read.csv(
  file.path(V2, "04_bootnet_edge_intervals.csv")
)

if (
  nrow(boot) != 918L ||
  any(boot$successful_values != 1000L) ||
  any(boot$lower > boot$upper)
) {
  stop("Invalid bootstrap summaries.")
}

nct <- read.csv(
  file.path(V2, "05_original_NCT_omnibus.csv")
)

if (
  nrow(nct) != 6L ||
  any(nct$permutations != 1000L)
) {
  stop("Unexpected original NCT settings.")
}

e <- read.csv(
  file.path(V2, "05_original_NCT_edges.csv")
)

for (label in unique(e$comparison)) {

  a <- e[
    e$comparison == label,
    ,
    drop = FALSE
  ]

  if (
    nrow(a) != 153L ||
    max(abs(
      a$p_BH_within153 -
      p.adjust(a$p_raw, "BH")
    )) > 1e-10 ||
    max(abs(
      a$p_Holm_within153 -
      p.adjust(a$p_raw, "holm")
    )) > 1e-10
  ) {
    stop(
      "Original within-comparison correction failed: ",
      label
    )
  }
}

ordinal <- read.csv(
  file.path(V2, "06_ordinal_matrix_diagnostics.csv")
)

if (
  nrow(ordinal) != 6L ||
  any(ordinal$min_eigenvalue <= 1e-8)
) {
  stop("Invalid ordinal matrices.")
}

models <- read.csv(
  file.path(
    V2,
    "02b_model_sample_and_convergence.csv"
  )
)

if (
  nrow(models) != 9L ||
  !all(models$converged)
) {
  stop("Selection model convergence failed.")
}

csv(
  dplyr::bind_rows(checks),
  "12_existing_results_inventory.csv"
)

save_environment("12")


# ============================================================
# F. EXPORT AGGREGATE RESULTS ONLY
# ============================================================

share <- file.path(
  OUT,
  "SHAREABLE_REQUIRED_CHECKS"
)

dir.create(share)

files <- list.files(
  OUT,
  pattern = "[.](csv|txt)$",
  full.names = TRUE
)

files <- files[
  !grepl(
    "RESTRICTED",
    basename(files),
    ignore.case = TRUE
  )
]

if (!all(file.copy(files, share))) {
  stop("Aggregate export failed.")
}

writeLines(
  c(
    "Original-submission diagnostics and VAS numerical checks only.",
    "No new permutations, correction families, motion exclusions or equivalence tests.",
    "Original Little MCAR and missingness diagnostics restricted to eligible QC-passing AND samples.",
    "Original selected-edge regressions remain post-hoc descriptive.",
    "Regression coefficients do not test network partial correlations.",
    "Family random intercept is unidentifiable within a wave with one observation per family.",
    "Scan provenance, MAR and practical absence of associations cannot be established by these tests.",
    "Inspect MCAR status/warnings and regression diagnostics before claiming completion."
  ),
  file.path(
    share,
    "12_READ_ME.txt"
  )
)

message(
  "Required checks complete. Send: ",
  share
)