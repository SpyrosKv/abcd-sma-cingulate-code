# Computational verification

## Locally executed workflow — 24 November 2026

The author executed the packaged workflow in RStudio on Windows, using the four prepared input files. All four input MD5 checksums matched the recorded R1 provenance, and the required variable checks passed (42 main-table variables, 9 MRI-QC variables, 3 eligibility variables and 3 education-source variables).

The complete `r1_all` route finished, including imputation, bootstrap, alternative family selections, domain comparisons, caregiver education, retention, account-holder analyses and report generation. A local comparison against the R1 reference established:

- All 14 settings matched.
- The same participants and node order in all six primary networks.
- Zero maximum absolute difference in adjusted network inputs and graph weights.
- Exact parsed content matches for all 109 compared analysis CSV tables. The inventory excluded restricted records, generated report copies and run-specific metadata/checksum files; it is not a claim of byte equality for every private intermediate object.

Primary complete-case sample sizes were FA 3,295 / 4,006 / 5,348 and T2-weighted 3,325 / 4,017 / 5,132, in visit order ses-02A / ses-04A / ses-06A.

The `r2_nct` stage was then run from the freshly reconstructed inputs. Its six omnibus rows, 918 edge-test rows and six validation rows were identical byte for byte to the earlier verified 10,000-permutation outputs. BH/Holm corrections, the plus-one p-value grid, domain edge counts and reconstruction of omnibus statistics from edge weights were independently checked. An earlier local audit also recalculated all 918 raw p-values directly from saved full-precision permutation arrays; these arrays remain private.

Base-R smoke tests and the synthetic NCT API check passed in the local R environment. Synthetic tests contain no ABCD records.

## Source scope of the verification

All R files in this release except `R/reports.R` are byte-identical to the corresponding files in the locally executed package. The reporting wrapper subsequently received descriptive changes: S3 extraction labels, explicit corrected-edge summaries, testing-family rationale, family-selection interpretation and neutral equivalence wording. Those changes do not fit models, change p-values, or change sample selection. Static syntax/source and HTML-fixture checks accompany preparation; the updated generated reports still require local rendering and disclosure review before submission.

Successful reproduction does not establish exchangeability, absence of selection bias, clinical importance, exact motion/FA source-scan correspondence, or validity in other software environments. Dense-network warnings and the package's paired-NCT limitation remain relevant. Initial maximum-likelihood selection models can warn about separation; the subsequent bias-reduced outputs are the reported replacements.

## Private records and environment

Source/data signatures, configurations, session information, raw input checksums and participant-level objects remain local. A source-file manifest records release bytes without publishing those private objects. Direct package versions observed during the run are documented in `software_environment.md` and `validated_package_versions.csv`; no historical environment or dependency lockfile has been inferred.

The public deposit and its permanent identifier are separate publication steps. A scientific reproducibility statement should specify that the reconstruction begins with the prepared input tables, not raw MRI images or an automatically recreated source export.
