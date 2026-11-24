# Input contract

No data are supplied. A researcher must have lawful access to the same Release 7.0
prepared tables and a current approved analysis environment. This package starts
from tabular exports; it does not acquire or process raw MRI images.

## Saved-R1 reviewer route

`existing_analysis_dir` is the directory above these result folders, not the
network RDS itself:

- `Paper1_revision_outputs_original_methods_v2`: `03_settings.rds`,
  `03_primary_networks_RESTRICTED.rds`, `05_original_NCT_omnibus.csv`,
  `05_original_NCT_edges.csv` for the new NCT; other saved outputs for reports.
- `Paper1_original_required_checks`: script-12 diagnostics and outputs.
- `Paper1_VAS_remaining_checks`: script-13–18 outputs, including its private
  and historical `SHAREABLE_VAS_CHECKS` subfolders.
- `Paper1_Retention_Results`: completed script-19a outputs and checksum marker.
- `Paper1_Account_Holder_Comparison`: script-20 early/later summaries.

The network object contains participant-level adjusted data and IDs. Keep it
local; do not send it to the code repository, Zenodo, OSF or ChatGPT.

## Full R1 reconstruction

Configure the actual files by contents, not by filename guesses:

1. `inputs$paper1`: prepared original Paper1 SAV with original identifiers,
   labelled visit codes, six SMA items, both sets of 12 MRI variables, demographic
   variables, family identifier, age, scanner serial and scanner software.
2. `inputs$mri_qc`: MRI QC extract with identifiers/visits, separate FA/T2
   recommended inclusion flags, both dMRI motion summaries and manual audit fields.
3. `inputs$sma_eligibility`: eligibility extract with identifiers/visits and the
   social-media account ownership gate, used to validate the original gate.
4. `inputs$education`: source SAV with identifiers/visits and the verified caregiver
   education variable `ab_g_dyn__cohort_edu__cgs` (expected categories 1–5).
5. `inputs$full_mri`: optional full MRI SAV used only for the `motion_source` audit.

Keys are `participant_id` and `session_id`. Visits are decoded from each file's
labels before joining; numeric visit codes must not be assumed comparable.
Duplicate or contradictory keys/gates stop the analysis. The original gate is
not silently backfilled from the validation extract.

SMA variables: `nt_y_stq__socmed_008` through `nt_y_stq__socmed_013`, with
labelled response categories 1–6; code 777 is treated as nonresponse.
Gate: `nt_y_stq__socmed_001` (yes=1, no=0).

MRI suffixes: `gcpd`, `gcpv`, `gsca`, `gscma`, `gscmp`, `scm`, in `lh`/`rh`.
FA prefix: `mr_y_dti__is__fa__wm__dst__`.
T2 prefix: `mr_y_smri__t2__wm__dst__`.
The variable suffix is `__lh_mean` or `__rh_mean`.

QC: `mr_y_qc__incl__dmri_indicator`, `mr_y_qc__incl__smri__t2_indicator`.
Motion: `mr_y_qc__raw__dmri__mot__ok_mean`, `mr_y_qc__raw__dmri__mot_mean`.
Manual audit: `mr_y_qc__post__man__t2_score`,
`mr_y_qc__post__man__dmri_score`, `mr_y_qc__post__man__fsurf_score`.

The exact covariate and demographic names are declared in the settings engine.
Inputs must retain original value labels, missing codes and units. If a header
does not match, stop and inspect the documented export; do not rename a different
measure to satisfy the code.

## Provenance limits

The initial export/merge that created the prepared original SAV is not reconstructed
by these scripts. Record release DOI/version, source table list, extraction filters
and any upstream processing from the data provider. Exact motion-to-FA scan
correspondence is a provenance question that statistical checks cannot establish.
No new claim of scan correspondence is made by this package.

References: NBDC Data Hub https://www.nbdc-datahub.org/data-access-process ;
ABCD Release 7.0 documentation https://docs.abcdstudy.org/v/7_0_0/ .
