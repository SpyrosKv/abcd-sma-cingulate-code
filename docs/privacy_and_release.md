# Data protection and release boundaries

The code deposit and controlled ABCD data are separate resources. The NBDC DUC
prohibits redistribution of raw or derived individual-level data, including
submission to third-party systems/large language models. Sharing code does not
grant data access to another researcher.

Before execution, confirm current institutional and DUC authorisation, named
authorised users, analysis scope and approved storage against the current
institutional records. Historical approval correspondence alone does not confirm
current access. The package does not grant access or extend an approval period.
Signed approvals and access correspondence remain institutional records and are
not part of the public code deposit.

Never publish the analysis working directory, even in a private GitHub repository.
Model and bootstrap objects can retain input measurements; an RDS extension or
absence of visible IDs is not evidence of safe sharing. Individual adjusted values,
imputed measurements and lists of selected participants remain controlled data.
Legacy folder names beginning SHAREABLE are only operational names.

Public aggregate results require an independent disclosure check. The reviewed
DUC specifies a minimum public cell size of 10. Small counts, percentages and
complementary totals can reveal a small cell even after its direct count is hidden.
This package does not silently suppress values because doing so could alter the
manuscript without author review; all generated tables remain local review files.

For release, use a fresh code-only repository from this ZIP. Run
`python tools/check_public_release.py` and inspect `git diff --cached --stat`
before the first push. The scanner is a guard, not proof of compliance; manually
review source text and the staged file list. If sensitive data were ever committed,
deleting the current file is insufficient because Git history retains old blobs.

No signed DUC, passwords, API keys, participant data, execution logs, outputs or
HTML reports are included in this distribution. Archive the verified code release
in Zenodo and cite its DOI after deposition. Link the new code deposit from the
existing OSF retrospective record without changing its historical chronology.
