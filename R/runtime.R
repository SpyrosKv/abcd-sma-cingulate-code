paper1_packages <- function() c("haven", "dplyr", "tidyselect", "lme4", "qgraph",
  "polycor", "mice", "bootnet", "NetworkComparisonTest", "brglm2",
  "detectseparation", "naniar", "xml2", "base64enc", "ggplot2", "pdftools", "renv")

paper1_path <- function(p, must = TRUE) {
  if (length(p) != 1L || is.na(p) || !nzchar(p)) stop("A configured path is empty.")
  normalizePath(path.expand(p), winslash = "/", mustWork = must)
}
paper1_inside <- function(path, parent) {
  a <- tolower(paper1_path(path, FALSE)); b <- tolower(paper1_path(parent, FALSE))
  identical(a, b) || startsWith(a, paste0(sub("/+$", "", b), "/"))
}
paper1_load_config <- function(file = "") {
  if (!nzchar(file)) file <- file.path(.paper1_repo, "config.local.R")
  if (!file.exists(file)) stop("Copy config.example.R to config.local.R and set work_dir / existing_analysis_dir first.")
  e <- new.env(parent = baseenv()); sys.source(file, e)
  if (!exists("paper1_config", e, inherits = FALSE) || !is.list(e$paper1_config))
    stop("Configuration must define paper1_config as a list.")
  cfg <- e$paper1_config
  cfg$work_dir <- paper1_path(cfg$work_dir, FALSE)
  if (paper1_inside(cfg$work_dir, .paper1_repo)) stop("work_dir must be outside the repository.")
  if (nzchar(cfg$existing_analysis_dir)) {
    cfg$existing_analysis_dir <- paper1_path(cfg$existing_analysis_dir)
    if (paper1_inside(cfg$existing_analysis_dir, .paper1_repo)) stop("Saved data must be outside the repository.")
  }
  if (!identical(as.integer(cfg$nct_permutations), 10000L))
    stop("The reviewer-revision stage is locked to 10000 permutations.")
  cfg
}
paper1_authorised <- function(cfg) {
  if (!isTRUE(cfg$authorised_data_use))
    stop("Confirm current DUC/institutional authorisation and approved storage, then set authorised_data_use=TRUE. Do not upload ABCD data here.")
  dir.create(cfg$work_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cfg$work_dir) || paper1_inside(paper1_path(cfg$work_dir), .paper1_repo))
    stop("The actual work directory must exist outside the repository, including after resolving links.")
}
paper1_versions <- function() {
  pkgs <- paper1_packages()
  data.frame(package = pkgs, version = vapply(pkgs, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, character(1)))
}
paper1_preflight <- function(cfg) {
  fs <- list.files(.paper1_repo, "[.]R$", recursive = TRUE, full.names = TRUE)
  fs <- fs[!grepl("config[.]local[.]R$", fs)]
  for (f in fs) parse(file = f, keep.source = FALSE)
  v <- paper1_versions(); print(v, row.names = FALSE)
  message("All packaged R files parsed. No participant data were loaded.")
  message("Missing packages: ", paste(v$package[is.na(v$version)], collapse = ", "))
  if (nzchar(cfg$existing_analysis_dir)) {
    p <- file.path(cfg$existing_analysis_dir, "Paper1_revision_outputs_original_methods_v2")
    n <- c("03_settings.rds", "03_primary_networks_RESTRICTED.rds", "05_original_NCT_omnibus.csv", "05_original_NCT_edges.csv")
    print(data.frame(required_file = n, exists = file.exists(file.path(p, n))), row.names = FALSE)
  }
  invisible(v)
}
paper1_baseline <- function(cfg) {
  parent <- if (nzchar(cfg$existing_analysis_dir)) cfg$existing_analysis_dir else file.path(cfg$work_dir, "analysis")
  p <- file.path(parent, "Paper1_revision_outputs_original_methods_v2")
  if (!dir.exists(p)) stop("Completed R1 folder not found. Point existing_analysis_dir to the folder above Paper1_revision_outputs_original_methods_v2, or run r1_all.")
  list(parent = parent, primary = p)
}
paper1_scope <- function(code_dir, choices = character()) {
  e <- new.env(parent = globalenv()); e$CODE_DIR <- code_dir
  e$file.choose <- local({
    paths <- choices; counter <- 0L
    function(...) {
      counter <<- counter + 1L
      if (counter > length(paths)) stop("Unexpected file prompt; the input routing must be checked.")
      paper1_path(paths[[counter]])
    }
  })
  # Keep all historical source() calls in the isolated stage environment.
  e$source <- local({ target <- e
    function(file, local = FALSE, ...) {
      dest <- if (is.environment(local)) local else if (isTRUE(local)) parent.frame() else target
      base::sys.source(file, envir = dest, keep.source = FALSE)
      invisible(NULL)
    }
  })
  e
}
paper1_missing_packages <- function(pkgs) {
  m <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(m)) stop("Install missing packages with setup_packages.R: ", paste(m, collapse = ", "))
}
paper1_source_helpers <- function(root) {
  paper1_missing_packages(c("haven", "dplyr", "tidyselect", "lme4", "qgraph", "polycor", "mice"))
  e <- paper1_scope(root)
  sys.source(file.path(.paper1_repo, "R", "engines", "00_Paper1_Settings_and_Functions.R"), e)
  e
}
paper1_atomic_rds <- function(x, path) {
  tmp <- tempfile("checkpoint_", tmpdir = dirname(path)); on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (file.exists(path)) stop("Refusing to overwrite a completed checkpoint: ", basename(path))
  if (!file.rename(tmp, path)) stop("Could not finalise checkpoint.")
}
paper1_find <- function(parent, name, roots = NULL) {
  if (is.null(roots)) roots <- file.path(parent, c("Paper1_revision_outputs_original_methods_v2",
    "Paper1_original_required_checks", "Paper1_VAS_remaining_checks", "Paper1_Retention_Results",
    "Paper1_Account_Holder_Comparison"))
  roots <- roots[dir.exists(roots)]
  files <- unique(unlist(lapply(roots, list.files, recursive = TRUE, full.names = TRUE)))
  hits <- files[basename(files) == name]
  if (!length(hits)) stop("Required result missing: ", name)
  if (length(unique(unname(tools::md5sum(hits)))) != 1L)
    stop("Conflicting copies of ", name, "; use one verified set of outputs.")
  hits[1L]
}
paper1_stage_map <- function() c(
  merge = "01_Paper1_Merge_Using_Visit_Labels.R",
  sample_flow = "02_Paper1_Sample_Flow_Missingness_and_Selection.R",
  missingness = "02b_Paper1_Bias_Reduced_Missingness_and_Selection.R",
  primary = "10_Paper1_Original_Methods_Reviewer_Corrections.R",
  checks = "12_Paper1_Original_Submission_Remaining_Checks.R",
  sensitivities = "13_Paper1_VAS_Equivalence_and_Final_Sensitivities.R",
  equivalence = "13b_Paper1_VAS_Equivalence_Exact_Request.R",
  wave_pmm = "14_Paper1_Matched_Wave_Specific_PMM.R",
  education = "15_Paper1_Parental_Education_Verification_and_Selection.R",
  education_audit = "15b_Paper1_Parental_Education_Matching_Audit.R",
  domains = "16_Paper1_VAS_Domain_Strength_NCT.R",
  families = "17_Paper1_VAS_Five_Family_Selections.R",
  distributions = "18_Paper1_VAS_Distributions_and_Exclusion_Groups.R",
  retention = "19a_Paper1_Retention_Original_Methods.R",
  account_holders = "20_Paper1_Early_Later_Account_Holders.R",
  motion_source = "02e_Paper1_Full_MRI_Source_Metadata_and_Motion_Audit.R")

paper1_reference_stage <- function(stage, cfg) {
  # Fresh reruns have their own code copy and results, outside the repository.
  root <- file.path(cfg$work_dir, "analysis")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  engines <- list.files(file.path(.paper1_repo, "R", "engines"), "[.]R$", full.names = TRUE)
  for (f in engines) {
    dest <- file.path(root, basename(f))
    if (file.exists(dest) && unname(tools::md5sum(dest)) != unname(tools::md5sum(f))) {
      if (!isTRUE(cfg$allow_existing_code_change)) stop("Execution copy changed: ", basename(f), ". Use a new work directory to avoid mixing runs.")
    }
    if (!file.copy(f, dest, overwrite = TRUE)) stop("Could not copy execution code.")
  }
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(root)
  settings <- file.path(root, "00_Paper1_Settings_and_Functions.R")
  v2 <- file.path(root, "Paper1_revision_outputs_original_methods_v2")
  input <- function(n) {
    x <- cfg$inputs[[n]]
    if (is.null(x) || !nzchar(x)) stop("Set inputs$", n, " for this stage.")
    x <- paper1_path(x)
    if (paper1_inside(x, .paper1_repo)) stop("Participant inputs must be outside the repository.")
    x
  }
  choices <- switch(stage,
    merge = c(input("paper1"), input("mri_qc"), input("sma_eligibility")),
    education = c(settings, input("education")),
    education_audit = c(settings, input("education")),
    equivalence = paper1_find(root, "13_unregularised_432_partial_correlations.csv"),
    retention = c(settings, file.path(v2, "03_primary_networks_RESTRICTED.rds")),
    account_holders = file.path(v2, "01_merged_RESTRICTED.rds"),
    motion_source = c(input("full_mri"), input("mri_qc")),
    primary = settings, checks = settings, sensitivities = settings, wave_pmm = settings,
    domains = settings, families = settings, distributions = settings, character())
  if (stage == "r1_reports") {
    cfg2 <- cfg; cfg2$existing_analysis_dir <- root
    return(paper1_reports(cfg2, revision = FALSE))
  }
  if (!stage %in% names(paper1_stage_map())) stop("Unknown reference stage.")
  e <- paper1_scope(root, choices)
  message("Running ", stage, ". All outputs remain in the authorised work directory.")
  sys.source(file.path(root, paper1_stage_map()[[stage]]), e)
  if (stage == "education") paper1_reference_stage("education_audit", cfg)
  invisible(TRUE)
}
