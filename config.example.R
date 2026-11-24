# Copy to config.local.R (ignored by Git), or set PAPER1_CONFIG to an external file.
# Use forward slashes on Windows. All data and outputs must be outside this repository.
paper1_config <- list(
  # Set TRUE only after confirming current institutional/DUC authorisation and storage.
  authorised_data_use = FALSE,
  work_dir = "",              # New, approved directory for execution and outputs.
  existing_analysis_dir = "", # R1 scripts folder containing the completed result folders.
  # Leave inputs empty when using the saved-R1 route. Required for r1_all only.
  inputs = list(paper1 = "", mri_qc = "", sma_eligibility = "",
                education = "", full_mri = ""),
  # Optional S2/S3 templates; plain self-contained templates are supplied otherwise.
  templates = list(S2 = "", S3 = ""),
  nct_permutations = 10000L,
  enforce_r1_sample_sizes = TRUE,
  allow_existing_code_change = FALSE
)
