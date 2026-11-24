# VAS item-distribution addition and disjoint exclusion-group selection contrasts.
# Does not fit new missingness models or change primary inclusion.
message('Choose the existing 00_Paper1_Settings_and_Functions.R.')
p00 <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R')stop('Wrong settings file.')
CODE_DIR <- dirname(p00);source(p00,local=TRUE)
V2 <- file.path(CODE_DIR,'Paper1_revision_outputs_original_methods_v2')
CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(CFG$pooled_wave_fixed,FALSE))stop('Unexpected primary configuration.')
d <- readRDS(file.path(V2,'01_merged_RESTRICTED.rds'))
cats <- readRDS(file.path(V2,'01_sma_categories.rds'))
OUT <- file.path(CODE_DIR,'Paper1_VAS_remaining_checks')
private <- file.path(OUT,'RESTRICTED');share <- file.path(OUT,'SHAREABLE_VAS_CHECKS')
dir.create(private,recursive=TRUE,showWarnings=FALSE);dir.create(share,showWarnings=FALSE)
education_file <- file.path(private,'15_verified_education_RESTRICTED.rds')
education_available <- file.exists(education_file)
if(education_available) {
 ed <- readRDS(education_file)
 check_keys(ed$data,'Verified education')
 ix <- match(paste(d$participant_id,d$visit),paste(ed$data$participant_id,ed$data$visit))
 if(anyNA(ix))stop('Education join no longer matches primary data.')
 d$verified_parental_education <- ed$data$education[ix]
}
mean_ok <- function(x)if(length(x))mean(x)else NA_real_
dist <- frequency <- counts <- contrasts <- list()
for(w in WAVES) {
 z <- d[d$visit==w,,drop=FALSE]
 for(v in SMA) {
  x <- z[[v]][z[[GATE]] %in% 1];a <- x[!is.na(x)]
  m2 <- if(length(a))mean((a-mean(a))^2)else NA_real_
  skew <- if(is.finite(m2)&&m2>0)mean((a-mean(a))^3)/m2^1.5 else NA_real_
  kurt <- if(is.finite(m2)&&m2>0)mean((a-mean(a))^4)/m2^2 else NA_real_
  dist[[paste(w,v)]] <- data.frame(visit=w,item=v,N_eligible=length(x),N_observed=length(a),
   N_missing=sum(is.na(x)),mean=mean_ok(a),SD=if(length(a)>1)sd(a)else NA_real_,
   skewness_moment=skew,kurtosis_moment=kurt,excess_kurtosis_moment=kurt-3,
   floor_percent=100*mean_ok(a==min(cats[[v]])),ceiling_percent=100*mean_ok(a==max(cats[[v]])))
  frequency[[paste(w,v)]] <- data.frame(visit=w,item=v,category=cats[[v]],
   N=vapply(cats[[v]],function(k)sum(a==k),numeric(1)),denominator_observed=length(a))
 }
 for(mod in names(MRI)) {
  f <- flags_for(z,mod)
  sma_problem <- !f$sma_complete
  mri_problem <- !f$qc_ok | !f$mri_complete
  group <- ifelse(f$included,'included',ifelse(sma_problem & mri_problem,'SMA_and_MRI_QC',
   ifelse(sma_problem,'SMA_only',ifelse(mri_problem,'MRI_QC_only','other_design_or_covariate_exclusion'))))
  for(pop in c('all_account_holders','selected_family_account_holders')) {
   use <- f$eligible & if(pop=='selected_family_account_holders')f$family_ok else TRUE
   zz <- z[use,,drop=FALSE];g <- group[use]
   levels <- c('included','SMA_only','MRI_QC_only','SMA_and_MRI_QC','other_design_or_covariate_exclusion')
   for(gr in levels) {
    take <- use & group==gr
    counts[[paste(w,mod,pop,gr)]] <- data.frame(visit=w,modality=mod,population=pop,group=gr,
     N=sum(take),N_total=sum(use),QC_failure=sum(take & z[[QC[[mod]]]] %in% 0),
     QC_unknown=sum(take & is.na(z[[QC[[mod]]]])),missing_MRI=sum(take & !f$mri_complete),
     missing_SMA=sum(take & !f$sma_complete),unselected_or_unknown_family=sum(take & !f$family_ok),
     invalid_adjustment_covariates=sum(take & !f$cov_ok),invalid_or_missing_motion=sum(take & !f$mot_ok))
   }
   variables <- c('age',DEMOS,SMA,MRI[[mod]],MOTION,
     if(education_available)'verified_parental_education'else character())
   for(gr in setdiff(levels,'included')) {
    inc <- g=='included';exc <- g==gr
    if(!sum(inc) || !sum(exc))next
    for(v in variables) {
     categorical <- v %in% c(DEMOS,'verified_parental_education')
     if(categorical) {
      x <- as.character(haven::zap_labels(zz[[v]]))
      for(k in sort(unique(x[!is.na(x)]))) {
       a <- x[inc];b <- x[exc];a <- a[!is.na(a)];b <- b[!is.na(b)]
       pa <- mean_ok(a==k);pb <- mean_ok(b==k)
       denominator <- sqrt((pa*(1-pa)+pb*(1-pb))/2)
       contrasts[[length(contrasts)+1]] <- data.frame(visit=w,modality=mod,population=pop,
        excluded_group=gr,variable=v,category=k,included_N=sum(inc),excluded_N=sum(exc),
        observed_included=length(a),observed_excluded=length(b),included_value=pa,excluded_value=pb,
        SMD=if(is.finite(denominator)&&denominator>0)(pa-pb)/denominator else NA_real_)
      }
     }else {
      x <- num(zz[[v]]);a <- x[inc];b <- x[exc]
      a <- a[is.finite(a)];b <- b[is.finite(b)]
      denominator <- if(length(a)>1 && length(b)>1)sqrt((var(a)+var(b))/2)else NA_real_
      contrasts[[length(contrasts)+1]] <- data.frame(visit=w,modality=mod,population=pop,
       excluded_group=gr,variable=v,category='continuous',included_N=sum(inc),excluded_N=sum(exc),
       observed_included=length(a),observed_excluded=length(b),included_value=mean_ok(a),excluded_value=mean_ok(b),
       SMD=if(is.finite(denominator)&&denominator>0)(mean(a)-mean(b))/denominator else NA_real_)
     }
    }
   }
  }
 }
}
count_table <- dplyr::bind_rows(counts)
totals <- count_table |> dplyr::group_by(visit,modality,population) |>
 dplyr::summarise(sum_groups=sum(N),denominator=unique(N_total),.groups='drop')
if(any(totals$sum_groups!=totals$denominator))stop('Exclusion groups do not partition the population.')
csv(dplyr::bind_rows(dist),'18_SMA_distributions_with_kurtosis.csv')
csv(dplyr::bind_rows(frequency),'18_SMA_category_frequencies.csv')
csv(count_table,'18_disjoint_exclusion_group_counts.csv')
csv(dplyr::bind_rows(contrasts),'18_included_vs_each_exclusion_group_SMD.csv')
csv(data.frame(education_available=education_available,
 source_variable=if(education_available)ed$variable else NA_character_), '18_education_status.csv')
writeLines(c('Group counts partition each population; detailed reasons may overlap within groups.',
 'MRI problem includes unavailable MRI and failed/unknown modality QC.',
 'Other group retains design/covariate/motion exclusions without SMA/MRI missingness.',
 'SMDs use observed values only; report observed group Ns.',
 'QC-failed MRI values are audit descriptions, not validated biological measurements.',
 'Moment kurtosis has normal reference 3; excess kurtosis has normal reference 0.',
 if(education_available)'Verified education included.'else 'Education pending: complete coding verification in 15, then rerun 18.'
),file.path(OUT,'18_READ_ME.txt'))
save_environment('18')
files <- list.files(OUT,pattern='^18_.*[.](csv|txt)$',full.names=TRUE)
if(!all(file.copy(files,share,overwrite=TRUE)))stop('Aggregate export failed.')
message('18 complete. Send ONLY: ',share)
