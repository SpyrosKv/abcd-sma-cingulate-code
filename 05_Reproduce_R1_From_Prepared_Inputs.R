# Optional, expensive full R1 reconstruction. Not required to repeat just the six NCTs.
# Run only in a NEW work_dir. Existing R1 outputs are not overwritten.
source("run_analysis.R")
paper1_run("r1_all")
