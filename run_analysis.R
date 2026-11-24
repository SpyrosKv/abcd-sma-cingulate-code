# Source this file; sourcing alone does not load participant data or run analyses.
.paper1_repo <- local({
  # Use the innermost matching source file, including when called from a test.
  files <- vapply(sys.frames(), function(e) {
    x <- e$ofile
    if (is.character(x) && length(x)) x[[1L]] else ""
  }, character(1))
  files <- files[nzchar(files) & basename(files) == "run_analysis.R"]
  root <- if (length(files)) dirname(normalizePath(tail(files, 1L), winslash = "/", mustWork = TRUE)) else
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "R", "runtime.R"))) stop("Cannot locate the analysis project root.")
  root
})
source(file.path(.paper1_repo, "R", "runtime.R"), local = TRUE)
source(file.path(.paper1_repo, "R", "reviewer_nct.R"), local = TRUE)
source(file.path(.paper1_repo, "R", "reviewer_summaries.R"), local = TRUE)
source(file.path(.paper1_repo, "R", "reports.R"), local = TRUE)

paper1_run <- function(stage = "check", config_file = Sys.getenv("PAPER1_CONFIG", "")) {
  cfg <- paper1_load_config(config_file)
  valid <- c("check", "r2_nct", "reviewer_summaries", "r2_reports", "r1_all",
             "merge", "sample_flow", "missingness", "primary", "checks",
             "sensitivities", "equivalence", "wave_pmm", "education",
             "domains", "families", "distributions", "retention", "account_holders",
             "r1_reports", "motion_source")
  if (!stage %in% valid) stop("Unknown stage. Choose: ", paste(valid, collapse = ", "))
  if (stage == "check") return(invisible(paper1_preflight(cfg)))
  paper1_authorised(cfg)
  if (stage == "r2_nct") return(invisible(paper1_nct_revision(cfg)))
  if (stage == "reviewer_summaries") return(invisible(paper1_reviewer_summaries(cfg)))
  if (stage == "r2_reports") return(invisible(paper1_reports(cfg, revision = TRUE)))
  if (stage == "r1_all") {
    for (nm in c("paper1", "mri_qc", "sma_eligibility", "education")) {
      f <- cfg$inputs[[nm]]
      if (is.null(f) || !nzchar(f) || !file.exists(f)) stop("Set the prepared input path inputs$", nm, " before starting the full reconstruction.")
      if (paper1_inside(f, .paper1_repo)) stop("Prepared inputs must be outside the repository.")
    }
    if (dir.exists(file.path(cfg$work_dir, "analysis", "Paper1_revision_outputs_original_methods_v2")))
      stop("Full reconstruction needs a fresh work_dir. Resume a known unfinished stage individually instead of overwriting the baseline.")
    # Deliberately excludes obsolete direct stages 04/05/07 and the old 5000 test.
    plan <- c("merge", "sample_flow", "missingness", "primary", "checks",
              "sensitivities", "equivalence", "wave_pmm", "education", "domains",
              "families", "distributions", "retention", "account_holders", "r1_reports")
    for (s in plan) paper1_reference_stage(s, cfg)
    return(invisible(TRUE))
  }
  invisible(paper1_reference_stage(stage, cfg))
}
