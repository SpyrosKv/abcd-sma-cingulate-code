# Record an environment after the local numerical verification has succeeded.
# This does not claim that the November environment is the historical October environment.
source("run_analysis.R")
cfg <- paper1_load_config(Sys.getenv("PAPER1_CONFIG", ""))
paper1_authorised(cfg)
paper1_missing_packages(paper1_packages())
dest <- file.path(cfg$work_dir, "software_environment")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
write.csv(paper1_versions(), file.path(dest, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(dest, "sessionInfo.txt"))
renv::snapshot(project = .paper1_repo, packages = paper1_packages(),
  lockfile = file.path(dest, "renv.lock"), prompt = FALSE)
message("Environment files saved locally: ", dest)
message("Review and include the verified renv.lock in the repository release; keep config.local.R private.")
