# Run separately; never installs or upgrades packages during an analysis stage.
source("run_analysis.R")
pkgs <- paper1_packages()
absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) install.packages(absent, repos = "https://cloud.r-project.org")
if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
  stop("Some dependencies are unavailable. See docs/software_environment.md.")
message("Dependencies available. Restart R if installation updated loaded libraries.")
message("Versions are recorded at execution. Historical package versions are not invented.")
