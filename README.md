# SMA symptoms and cingulate white-matter MRI networks

Analysis code for **Social Media Addiction Symptoms and Cingulate White-Matter MRI Measures in Adolescent Social Media Users: A Cross-Wave Network Analysis of the ABCD Study**.

The repository contains R code and documentation. ABCD participant-level data, private configurations, intermediate objects, generated reports and research approvals are not distributed here.

## Scope and verification

The statistical workflow was executed locally on 24 November 2026 from four prepared SAV inputs. All 14 checked settings, the participant memberships and adjusted values of all six primary networks, and the graph weights matched the R1 reference. All 109 compared analysis CSV tables matched exactly. The six subsequent 10,000-permutation comparisons produced three CSV files identical byte for byte to the independently checked R2 outputs.

This is computational verification in the recorded Windows/R environment, not validation of all statistical assumptions or of upstream data-extraction provenance. The reporting wrapper contains subsequent wording corrections; see [validation](docs/validation.md) and [changes](CHANGELOG.md).

## Requirements

- Authorised access to the matching ABCD Study Release 7.0 prepared tables and an approved local analysis environment.
- R and the packages listed in `setup_packages.R`; the versions observed during verification are documented in `docs/software_environment.md`.
- Four inputs: the prepared Paper 1 table, MRI QC extract, questionnaire-eligibility extract, and caregiver-education source.

See [input requirements](docs/input_data_requirements.md) for exact fields and coding. The pipeline starts from prepared tabular exports. It does not download ABCD data, process raw MRI scans, or reconstruct the provider's upstream export that produced these SAV files.

## Full reproduction from prepared inputs

Open `ABCD_SMA_Cingulate.Rproj`. Copy `config.example.R` to `config.local.R`. Set an approved **new external work directory**, leave `existing_analysis_dir` empty, and provide the four required `inputs` paths. Set `authorised_data_use = TRUE` only when the applicable access and storage requirements are met. The configuration remains private and is ignored by Git.

```r
source("setup_packages.R")
source("run_analysis.R")
paper1_run("check")
source("tests/smoke_test.R")
source("tests/nct_synthetic_check.R")
paper1_run("r1_all")
paper1_run("r2_nct")
paper1_run("reviewer_summaries")
paper1_run("r2_reports")
```

This route does **not** require another researcher's pre-existing R1 results. The full R1 reconstruction creates its own inputs for the R2 stage. Bootstraps, imputations and family-selection analyses can take several hours. Full reconstruction requires a fresh work directory; diagnose an interruption before resuming individual stages.

## Repeat the R2 comparisons from saved R1 inputs

For an existing completed analysis, set `existing_analysis_dir` to the directory above `Paper1_revision_outputs_original_methods_v2`. Keep `work_dir` separate and outside this repository. The six comparisons need the saved network object, settings and original comparison CSVs; reports additionally need the other completed R1 output folders. Then run `r2_nct`, `reviewer_summaries` and `r2_reports`.

The root scripts `01` through `05` provide corresponding entry points. Do not run the scripts in `R/engines` in filename order: the runner selects the final implementations, including replacements for superseded scripts. Historical internal filenames are retained to preserve the analysis provenance.

## Statistical specifications

- Visits: `ses-02A`, `ses-04A`, `ses-06A`, decoded separately from each source's visit labels before joining by participant and visit.
- Account ownership assessed at each visit; structural non-administration excluded from SMA imputation. Separate FA/T2-weighted inclusion flags are applied.
- One participant per family, selected once across visits with seed 20261003.
- Pooled primary residualisation without visit fixed effects or age-by-visit terms. Fixed-effect predictions and imaging scanner-intercept predictions are removed; the family contribution is retained. FA adjustment includes the specified OK-scan motion summary.
- Six 18-node Pearson EBICglasso networks, gamma .5; 1,000 nonparametric and case-dropping bootstrap samples.
- Original primary paired NCTs use 1,000 permutations; revised primary NCTs use 10,000. BH and Holm corrections are applied across 153 edges within each comparison. Omnibus p-values remain unadjusted and exploratory.
- Domain comparisons retain their separate family of 18 tests and their 1,000-permutation implementation. Family-selection and other R1 sensitivity comparisons also retain their original permutation counts.
- The workflow includes missingness/selection, residual-scale PMM, ordinal, family and motion sensitivities, caregiver education, retention and early/later account-holder comparisons.
- Exploratory equivalence uses pairwise unadjusted TOST at bounds of ±.05 and ±.03, and an all-72-pair intersection-union criterion. Detailed outputs also retain separately adjusted pairwise tests.

MRI measures summarise white matter adjacent to Destrieux cingulate cortical parcels; they are not tractography-defined cingulum measurements. Cross-wave comparisons concern between-person networks, not within-person developmental change or causal effects.

## Outputs and privacy

Outputs are written outside the repository. The work directory contains controlled participant-level information, including IDs, residuals, imputed values, model and bootstrap objects. **Never publish that directory**, including directories whose historical names contain `SHAREABLE`.

The code makes no network upload. Generated reports are local review outputs and require disclosure review, including small participant cells and reconstruction through percentages or totals, before journal or repository release. The manuscript, response letter and S1 Word document are edited separately.

Before committing this source tree, run:

```sh
python tools/check_public_release.py
```

The file allowlist and scanner supplement manual review; they do not certify institutional or legal compliance. See [release boundaries](docs/privacy_and_release.md).

## Citation and reuse

Author attribution is provided in `CITATION.cff`. After a versioned repository release is archived in a citable repository, add the actual repository URL and DOI to that file and cite the archived version in the manuscript. No DOI is asserted by this source package.

No software reuse licence has been specified. Public access to source code does not itself grant additional reuse rights; a licence may be added by the relevant rights holders. ABCD data access and reuse remain governed separately by the applicable data-use terms.
