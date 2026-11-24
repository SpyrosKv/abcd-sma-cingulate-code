# Execution order and result mapping

The source engines retain the authors' original calculations. An isolated runner
provides explicit file paths in place of file-selection dialogs and keeps outputs
outside the repository. Never execute all engine filenames in order.

| Stage | Engine | Main purpose |
|---|---|---|
| merge | 01 | Decode visits, validate joins, retain fixed family selection |
| sample_flow | 02 | Sequential exclusions, observed distributions, missingness summaries |
| missingness | 02b | Bias-reduced models superseding naive script-02 logistic inference |
| primary | 10, calling 03 and 06 | R1 pooled networks, bootnet, paired NCT1000, sensitivity, PMM |
| checks | 12 | Reconstruction, wave-specific adjustment, eligible MCAR, selected-edge regressions |
| sensitivities | 13 | All-pair uncertainty, alternate residualisation, wave bootstraps/NCT |
| equivalence | 13b | Final pairwise TOST and all-72 intersection–union summary |
| wave_pmm | 14 | Matched wave-specific PMM comparisons |
| education | 15 and 15b | Caregiver-education matching, coverage and selection contrasts |
| domains | 16 | 18 paired domain-strength tests, separate BH/Holm family |
| families | 17 | Five fixed family selections, networks, bootstraps and paired NCT |
| distributions | 18 | Item distributions and disjoint exclusion-group comparisons |
| retention | 19a | Three retention outcomes, original predictor construction |
| account_holders | early/later script 20 | Descriptive source-sample uptake comparisons |
| r1_reports | 19b plus final equivalence report patch | S2/S3 and S1 statistical figures |
| r2_nct | R/reviewer_nct.R | Six primary paired NCT comparisons at 10000 permutations |
| reviewer_summaries | R/reviewer_summaries.R | Numerical consistency and inference-sensitivity summaries |
| r2_reports | R/reports.R | New primary comparison tables, final equivalence wording and reviewer clarifications |

## R1 versus R2

The R1 baseline remains at 1000 primary NCT permutations and is never overwritten.
Only the six primary comparisons are repeated at 10000 for R2. Domain and
family/wave-specific sensitivities remain at their R1 settings; they are not
silently promoted to higher resolution. Edge-level inference in those sensitivities
retains the same permutation-resolution limitation and must be described accordingly.

The original script 11 combined checks now covered by script 12 with an older
5000-permutation repeat. It is not required by the supported workflow. The older
direct bootstrap/permutation/PMM scripts and optional-equivalence script are
excluded to prevent accidental use of superseded methods. The retained 13 engine
also generates earlier Holm-adjusted equivalence audit tables; final claims use
the later 13b/IUT outputs, as in the submitted R1.

New report snapshots replace only the primary omnibus and edge-comparison CSVs.
Statistical point estimates, domain tests and other sensitivity outputs are read
from the completed R1. The reports contain no automatic assurance that all aggregate
tables meet disclosure rules. S1 Word/manuscript edits require explicit author review.

## Software/source checks

The bundled execution engines are modified only to align settings defaults to
R1, keep temporary SAV extraction within approved storage and remove the blanket
copy of restricted CSVs into a folder labelled shareable. Method-specific engine
calculations are retained. Source hashes and changes are in source_provenance.csv.
The runner records execution software separately and stops rather than silently
accepting changed samples or observed edge weights in the new comparison stage.
