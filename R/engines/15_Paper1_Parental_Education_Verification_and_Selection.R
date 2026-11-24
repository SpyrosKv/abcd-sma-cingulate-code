# First run: leave the three settings below unchanged; export metadata only.
# Then confirm the exact variable and nonresponse codes before calculating SMDs.
# No guessed mapping, baseline carry-forward, or new covariate adjustment.
EDUCATION_VARIABLE <- 'ab_g_dyn__cohort_edu__cgs'
EDUCATION_MISSING_CODES <- numeric(0)
EDUCATION_CODING_CONFIRMED <- TRUE

message('Choose the existing 00_Paper1_Settings_and_Functions.R.')
p00 <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R') stop('Wrong settings file.')
CODE_DIR <- dirname(p00);source(p00,local=TRUE)
V2 <- file.path(CODE_DIR,'Paper1_revision_outputs_original_methods_v2')
CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(CFG$pooled_wave_fixed,FALSE)) stop('Wrong v2 configuration.')
d <- readRDS(file.path(V2,'01_merged_RESTRICTED.rds'))
OUT <- file.path(CODE_DIR,'Paper1_VAS_remaining_checks')
share <- file.path(OUT,'SHAREABLE_VAS_CHECKS')
dir.create(share,recursive=TRUE,showWarnings=FALSE)
message('Choose a SOURCE .sav that may contain parental education, not an MRI-only QC extract.')
p_source <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
header <- haven::read_sav(p_source,n_max=0)
labels <- vapply(header,function(x){a<-attr(x,'label');if(is.null(a))'' else paste(a,collapse=' ')},character(1))
dictionary <- data.frame(variable=names(header),description=unname(labels))
candidate <- grepl('educ|degree|qualification',paste(dictionary$variable,dictionary$description),ignore.case=TRUE)
csv(dictionary[candidate,,drop=FALSE],'15_education_candidate_dictionary.csv')
print(dictionary[candidate,,drop=FALSE],row.names=FALSE)

if(!nzchar(EDUCATION_VARIABLE)) {
 message('15 metadata only. Send 15_education_candidate_dictionary.csv before choosing a field.')
} else {
 if(!EDUCATION_VARIABLE %in% names(header)) stop('Chosen education field absent.')
 required <- c('participant_id','session_id',EDUCATION_VARIABLE)
 if(!all(required %in% names(header))) stop('Source needs participant_id and labelled session_id.')
 raw <- read_sav_source(p_source,required)
 v <- raw[[EDUCATION_VARIABLE]]
 lab <- attr(v,'labels')
 levels <- data.frame(value=if(is.null(lab))numeric(0)else unname(lab),
   label=if(is.null(lab))character(0)else names(lab))
 csv(levels,'15_education_value_labels.csv')
 counts <- as.data.frame(table(num(v),useNA='ifany'));names(counts) <- c('value','N')
 csv(counts,'15_education_source_value_counts.csv')
 if(!EDUCATION_CODING_CONFIRMED) {
  message('15 coding review only. Confirm value labels and missing codes before setting TRUE.')
 } else {
  raw <- decode_visit(raw,'parental_education_source')
  education <- num(raw[[EDUCATION_VARIABLE]])
  education[education %in% EDUCATION_MISSING_CODES] <- NA_real_
  if(any(!is.finite(education[!is.na(education)])) || any(education<0,na.rm=TRUE))
   stop('Invalid/negative codes remain; check missing-value coding.')
  key <- paste(d$participant_id,d$visit);source_key <- paste(raw$participant_id,raw$visit)
  if(anyDuplicated(source_key)) stop('Duplicate education source keys.')
  ix <- match(key,source_key);ed <- education[ix]
  # Restricted join for the separate exclusion-group contrasts in script 18.
  private <- file.path(OUT,'RESTRICTED')
  dir.create(private,showWarnings=FALSE)
  saveRDS(list(variable=EDUCATION_VARIABLE,
    data=data.frame(participant_id=d$participant_id,visit=d$visit,education=ed)),
    file.path(private,'15_verified_education_RESTRICTED.rds'))
  coverage <- result <- list()
  for(w in WAVES) for(mod in names(MRI)) {
   take <- d$visit==w;z <- d[take,,drop=FALSE];xx <- ed[take]
   f <- flags_for(z,mod)
   for(pop in c('all_account_holders','selected_family_account_holders')) {
    use <- f$eligible & if(pop=='selected_family_account_holders')f$family_ok else TRUE
    x <- xx[use];inc <- f$included[use]
    coverage[[length(coverage)+1]] <- data.frame(visit=w,modality=mod,population=pop,
     eligible_N=sum(use),matched_source=sum(!is.na(ix[take][use])),
     observed_education=sum(!is.na(x)),included_N=sum(inc),excluded_N=sum(!inc))
    for(level in sort(unique(x[!is.na(x)]))) {
     a <- mean(x[inc]==level,na.rm=TRUE);b <- mean(x[!inc]==level,na.rm=TRUE)
     denominator <- sqrt((a*(1-a)+b*(1-b))/2)
     result[[length(result)+1]] <- data.frame(visit=w,modality=mod,population=pop,
      variable=EDUCATION_VARIABLE,category=level,included_N=sum(inc),excluded_N=sum(!inc),
      observed_included=sum(!is.na(x[inc])),observed_excluded=sum(!is.na(x[!inc])),
      included_proportion=a,excluded_proportion=b,
      SMD=if(is.finite(denominator)&&denominator>0)(a-b)/denominator else NA_real_)
    }
   }
  }
  csv(dplyr::bind_rows(coverage),'15_education_join_and_observed_counts.csv')
  if(!length(result))stop('No observed education categories at target visits. No carry-forward assumed.')
  csv(dplyr::bind_rows(result),'15_parental_education_selection_SMD.csv')
  message('15 education comparisons complete; inspect visit coverage before reporting.')
 }
}
save_environment('15')
files <- list.files(OUT,pattern='^15_.*[.](csv|txt)$',full.names=TRUE)
if(!all(file.copy(files,share,overwrite=TRUE))) stop('Aggregate export failed.')
message('Send ONLY aggregate 15 files from: ',share)
