# 15b: έλεγχος αντιστοίχισης εκπαίδευσης.
# Δεν αλλάζει IDs, επισκέψεις, τιμές ή αναλυτικά δείγματα.

message('Choose 00_Paper1_Settings_and_Functions.R')
p00 <- normalizePath(file.choose(), winslash='/', mustWork=TRUE)

if (basename(p00) != '00_Paper1_Settings_and_Functions.R') {
  stop('Wrong settings file.')
}

CODE_DIR <- dirname(p00)
source(p00, local=TRUE)

V2 <- file.path(
  CODE_DIR, 'Paper1_revision_outputs_original_methods_v2'
)

CFG <- readRDS(file.path(V2, '03_settings.rds'))

if (!identical(CFG$pooled_wave_fixed, FALSE)) {
  stop('Wrong v2 configuration.')
}

d <- readRDS(file.path(V2, '01_merged_RESTRICTED.rds'))
check_keys(d, 'merged analysis data')

OUT <- file.path(CODE_DIR, 'Paper1_VAS_remaining_checks')
share <- file.path(OUT, 'SHAREABLE_VAS_CHECKS')
dir.create(share, recursive=TRUE, showWarnings=FALSE)

export <- function(x, name) {
  write.csv(
    x, file.path(share, name),
    row.names=FALSE, na=''
  )
}

message(
  'Choose the SAME dataset.sav used in script 15, containing education.'
)

p_source <- normalizePath(
  file.choose(), winslash='/', mustWork=TRUE
)

edvar <- 'ab_g_dyn__cohort_edu__cgs'

s <- read_sav_source(
  p_source, c('participant_id', 'session_id', edvar)
)

s <- decode_visit(s, 'education_source_15b')

export(
  attr(s, 'visit_mapping'),
  '15b_source_visit_code_mapping.csv'
)

ed <- num(s[[edvar]])

if (any(!is.na(ed) & !(ed %in% 1:5))) {
  stop('Unexpected education codes: review before continuing.')
}

key <- function(id, visit) {
  paste(id, visit, sep='\r')
}

ix <- match(
  key(d$participant_id, d$visit),
  key(s$participant_id, s$visit)
)

matched <- !is.na(ix)
observed <- matched & !is.na(ed[ix])
present_id <- d$participant_id %in% s$participant_id

# Έλεγχος για πιθανές διαφορές από κενά γύρω από τα IDs.
# Δεν διορθώνει και δεν χρησιμοποιεί αυτά τα IDs για την ανάλυση.
trim_ix <- match(
  key(trimws(d$participant_id), d$visit),
  key(trimws(s$participant_id), s$visit)
)

reason <- ifelse(
  observed,
  'matched_observed_education',
  ifelse(
    matched,
    'matched_missing_education',
    ifelse(
      present_id,
      'no_same_visit_ID_exists_at_other_visit',
      'ID_absent_from_source'
    )
  )
)

# Αναπαραγωγή της αντιστοίχισης που αποθηκεύτηκε στο 15.
verified <- file.path(
  OUT, 'RESTRICTED', '15_verified_education_RESTRICTED.rds'
)

if (!file.exists(verified)) {
  stop('Complete script 15 first: verified education file missing.')
}

v <- readRDS(verified)

if (!identical(v$variable, edvar)) {
  stop('Different verified education field.')
}

if (anyDuplicated(v$data[c('participant_id', 'visit')])) {
  stop('Duplicate verified keys.')
}

vi <- match(
  key(d$participant_id, d$visit),
  key(v$data$participant_id, v$data$visit)
)

same_values <- isTRUE(all.equal(
  as.numeric(v$data$education[vi]),
  ed[ix],
  check.attributes=FALSE
))

if (anyNA(vi) || !same_values) {
  stop(
    'Education does not reproduce script 15. Check source choice before proceeding.'
  )
}

# Διαθέσιμες εγγραφές και εκπαίδευση ανά επίσκεψη στο source.
visit_summary <- lapply(sort(unique(s$visit)), function(w) {
  take <- s$visit == w

  data.frame(
    visit=w,
    source_records=sum(take),
    observed_education=sum(!is.na(ed[take])),
    missing_education=sum(is.na(ed[take]))
  )
})

export(
  dplyr::bind_rows(visit_summary),
  '15b_source_visit_coverage.csv'
)

coverage <- groups <- other <- list()

for (w in WAVES) {
  for (mod in names(MRI)) {

    take <- which(d$visit == w)
    z <- d[take, , drop=FALSE]
    f <- flags_for(z, mod)

    for (pop in c(
      'all_account_holders',
      'selected_family_account_holders'
    )) {

      use <- f$eligible & if (
        pop == 'selected_family_account_holders'
      ) f$family_ok else TRUE

      # Χωριστά για συμπεριληφθέντες και αποκλεισμένους.
      for (group in c('included', 'excluded')) {

        rows <- take[
          use & if (group == 'included') {
            f$included
          } else {
            !f$included
          }
        ]

        for (r in c(
          'matched_observed_education',
          'matched_missing_education',
          'no_same_visit_ID_exists_at_other_visit',
          'ID_absent_from_source'
        )) {

          groups[[length(groups) + 1]] <- data.frame(
            visit=w,
            modality=mod,
            population=pop,
            group=group,
            reason=r,
            N=sum(reason[rows] == r),
            group_N=length(rows)
          )
        }
      }

      rows <- take[use]

      coverage[[length(coverage) + 1]] <- data.frame(
        visit=w,
        modality=mod,
        population=pop,
        eligible_N=length(rows),
        matched_source=sum(matched[rows]),
        observed_education=sum(observed[rows]),
        included_N=sum(f$included[use]),
        excluded_N=sum(!f$included[use]),
        unmatched_with_trimmed_exact_match=sum(
          !matched[rows] & !is.na(trim_ix[rows])
        )
      )

      # Αν δεν υπάρχει ίδια επίσκεψη, σε ποιες άλλες υπάρχει το ID;
      absent <- rows[!matched[rows] & present_id[rows]]

      for (sw in sort(unique(s$visit))) {

        ids <- s$participant_id[s$visit == sw]

        other[[length(other) + 1]] <- data.frame(
          analysis_visit=w,
          modality=mod,
          population=pop,
          source_visit=sw,
          unmatched_IDs_present_at_source_visit=sum(
            d$participant_id[absent] %in% ids
          )
        )
      }
    }
  }
}

# Επαλήθευση ότι αναπαράγονται τα συνοπτικά counts του 15.
a <- dplyr::bind_rows(coverage)

old_path <- file.path(
  share, '15_education_join_and_observed_counts.csv'
)

if (!file.exists(old_path)) {
  stop('Script 15 aggregate counts missing.')
}

old <- read.csv(old_path)

cols <- c(
  'visit', 'modality', 'population',
  'eligible_N', 'matched_source', 'observed_education',
  'included_N', 'excluded_N'
)

sort_rows <- function(x) {
  x[
    order(x$visit, x$modality, x$population),
    cols, drop=FALSE
  ]
}

if (!isTRUE(all.equal(
  sort_rows(a), sort_rows(old),
  check.attributes=FALSE
))) {
  stop('Counts do not reproduce 15: check source and settings.')
}

export(
  a,
  '15b_matching_reproduction_and_ID_format.csv'
)

export(
  dplyr::bind_rows(groups),
  '15b_education_missingness_by_inclusion_group.csv'
)

export(
  dplyr::bind_rows(other),
  '15b_unmatched_IDs_other_visit_presence.csv'
)

export(
  data.frame(
    source_file=basename(p_source),
    source_md5=unname(tools::md5sum(p_source)),
    analysis_md5=unname(tools::md5sum(
      file.path(V2, '01_merged_RESTRICTED.rds')
    )),
    reproduced_15_values=TRUE,
    reproduced_15_counts=TRUE
  ),
  '15b_input_audit.csv'
)

writeLines(
  c(
    'Diagnostic only: no ID cleaning, visit recoding, carry-forward or imputation.',
    'Visits decoded from each source SPSS labels; this checks internal consistency, not original extraction provenance.',
    'Other-visit counts overlap: the same ID can appear at several source visits.',
    'Presence at another visit does not authorise copying its education value.',
    'A matching source/result does not prove absent records were never collected.'
  ),
  file.path(share, '15b_READ_ME.txt')
)

writeLines(
  capture.output(sessionInfo()),
  file.path(share, '15b_sessionInfo.txt')
)

message('15b complete. Send all 15b_ files from: ', share)