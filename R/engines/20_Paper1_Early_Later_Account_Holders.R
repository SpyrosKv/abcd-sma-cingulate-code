# Paper 1: descriptive early/later account-holder comparison.
# Uses existing merged data before account/QC exclusions. No networks are refitted.
local({
  if (!requireNamespace('haven', quietly=TRUE)) stop('Install haven first.')
  message('Choose 01_merged_RESTRICTED.rds from Paper1_revision_outputs_original_methods_v2.')
  p <- normalizePath(file.choose(), winslash='/', mustWork=TRUE)
  if (basename(p) != '01_merged_RESTRICTED.rds' ||
      basename(dirname(p)) != 'Paper1_revision_outputs_original_methods_v2')
    stop('Select the v2 merged RDS, not a SAV, residual file or network object.')
  d <- readRDS(p); code_dir <- dirname(dirname(p))
  out <- file.path(code_dir, 'Paper1_Account_Holder_Comparison',
                   'SHAREABLE_SUMMARIES')
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  put <- function(x, name) write.csv(x, file.path(out, name), row.names=FALSE, na='')
  gate <- 'nt_y_stq__socmed_001'
  demos <- c('ab_g_stc__cohort_sex', 'ab_g_stc__cohort_race__nih',
             'ab_g_stc__cohort_ethn', 'ab_g_dyn__cohort_income__hhold__3lvl',
             'ab_g_dyn__design_site')
  waves <- c('ses-02A','ses-04A','ses-06A')
  required <- c('participant_id','visit','selected_family','family','age',gate,demos)
  if (!is.data.frame(d) || !all(required %in% names(d)))
    stop('Merged data lack required fields: ', paste(setdiff(required,names(d)),collapse=', '))
  d$participant_id <- trimws(as.character(d$participant_id))
  d$visit <- as.character(d$visit)
  if (anyNA(d$participant_id) || any(d$participant_id=='') || anyNA(d$visit) ||
      any(!d$visit %in% waves) || anyDuplicated(d[c('participant_id','visit')]))
    stop('Invalid participant/visit keys. Resolve before running this comparison.')
  g <- as.numeric(haven::zap_labels(d[[gate]]))
  if (any(!is.na(g) & !g %in% c(0,1))) stop('Account status must be 0, 1 or NA.')
  if (!any(g==0,na.rm=TRUE)) stop('No non-account holders: this may be an eligibility-filtered file.')
  if (!is.logical(d$selected_family) || anyNA(d$selected_family))
    stop('selected_family must be a complete logical field.')
  ids <- sort(unique(d$participant_id))
  selected <- vapply(ids,function(id) {
    x <- unique(d$selected_family[d$participant_id==id])
    if(length(x)!=1L) stop('Family selection changes across visits.')
    x
  },logical(1))
  sf <- unique(d[d$selected_family,c('participant_id','family')])
  if(anyNA(sf$family) || anyDuplicated(sf$participant_id) || anyDuplicated(sf$family))
    stop('Invalid fixed one-participant-per-family selection.')
  keys <- paste(d$participant_id,d$visit,sep='|')
  ix <- lapply(waves,function(w) match(paste(ids,w,sep='|'),keys))
  acc <- sapply(ix,function(i) g[i]); colnames(acc) <- waves
  a <- acc[,1]; b <- acc[,2]; c <- acc[,3]
  later_yes <- b %in% 1 | c %in% 1
  later_both_no <- b %in% 0 & c %in% 0
  group <- rep('Unclassified_missing_status',length(ids))
  group[is.na(a)] <- 'Unknown_initial_account_status'
  group[a %in% 0 & later_both_no] <- 'No_account_at_all_three_observed_visits'
  group[a %in% 0 & later_yes] <- 'Later_observed_account_holder'
  group[a %in% 1] <- 'Early_observed_account_holder'
  group_levels <- c('Early_observed_account_holder','Later_observed_account_holder',
    'No_account_at_all_three_observed_visits','Unknown_initial_account_status',
    'Unclassified_missing_status')

  ed_path <- file.path(code_dir,'Paper1_VAS_remaining_checks','RESTRICTED',
                       '15_verified_education_RESTRICTED.rds')
  if(!file.exists(ed_path)) {
    message('Choose the existing 15_verified_education_RESTRICTED.rds (from script 15).')
    ed_path <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
  }
  if(basename(ed_path)!='15_verified_education_RESTRICTED.rds')
    stop('Wrong verified education file.')
  ed <- readRDS(ed_path)
  if(!is.list(ed) || !is.data.frame(ed$data) ||
     !all(c('participant_id','visit','education') %in% names(ed$data)))
    stop('Unexpected verified education object.')
  ek <- paste(ed$data$participant_id,ed$data$visit,sep='|')
  if(anyNA(ed$data$participant_id) || anyNA(ed$data$visit) || anyDuplicated(ek))
    stop('Invalid education keys.')
  ex <- match(keys,ek)
  if(anyNA(ex)) stop('Education audit does not cover all merged participant/visit keys.')
  d$verified_parental_education <- ed$data$education[ex]
  z <- d[ix[[1]],,drop=FALSE] # all values refer to ses-02A; no carry-forward
  variables <- c('age',demos,'verified_parental_education')
  dictionary <- handling <- value_labels <- list()
  for(v in variables) {
    x <- z[[v]]; labs <- attr(x,'labels'); desc <- attr(x,'label')
    dictionary[[length(dictionary)+1]] <- data.frame(variable=v,
      label=if(is.null(desc)) v else paste(desc,collapse=' '),
      education_source=if(v=='verified_parental_education') ed$variable else '')
    codes <- numeric(0)
    if(!is.null(labs)) {
      value_labels[[length(value_labels)+1]] <- data.frame(variable=v,
        code=as.character(unname(labs)),label=names(labs))
      # Only explicitly labelled nonresponse categories are removed; export the decisions.
      drop <- grepl('declin|refus|don.?t know|do not know|prefer not|not applicable|missing|not asked|not administered',
                    names(labs),ignore.case=TRUE)
      codes <- unname(labs[drop])
      if(any(drop)) handling[[length(handling)+1]] <- data.frame(variable=v,
        code=as.character(codes),label=names(labs)[drop])
    }
    if(v=='age') {
      val <- as.numeric(haven::zap_labels(x)); val[!is.finite(val)] <- NA_real_
      if(any(val<=0,na.rm=TRUE)) stop('Invalid age. Check source units/coding.')
    } else val <- as.character(haven::zap_labels(x))
    val[as.character(val) %in% as.character(codes)] <- NA
    z[[v]] <- val
  }
  counts <- patterns <- comparisons <- list()
  avg <- function(x) if(length(x)) mean(x) else NA_real_
  pattern <- apply(acc,1,function(x) paste(ifelse(is.na(x),'unknown',x),collapse='/'))
  for(pop in c('all_source_participants','fixed_family_selection')) {
    use <- if(pop=='fixed_family_selection') selected else rep(TRUE,length(ids))
    for(gr in group_levels) counts[[length(counts)+1]] <- data.frame(
      population=pop,group=gr,N=sum(use & group==gr),total_N=sum(use),
      initial_visit_record_N=sum(use & group==gr & !is.na(ix[[1]])))
    pt <- as.data.frame(table(group=group[use],account_pattern=pattern[use]),stringsAsFactors=FALSE)
    pt <- pt[pt$Freq>0,,drop=FALSE]; names(pt)[3] <- 'N'; pt$population <- pop
    patterns[[length(patterns)+1]] <- pt
    early <- use & group=='Early_observed_account_holder'
    later <- use & group=='Later_observed_account_holder'
    if(!any(early) || !any(later)) stop('One comparison group is empty: ',pop)
    for(v in variables) {
      x <- z[[v]]; aa <- x[early]; bb <- x[later]
      aa <- aa[!is.na(aa)]; bb <- bb[!is.na(bb)]
      categories <- if(v=='age') 'continuous' else sort(unique(c(aa,bb)))
      if(!length(categories)) categories <- 'no_observed_values'
      for(k in categories) {
        if(v=='age') {
          va <- avg(aa); vb <- avg(bb)
          sa <- if(length(aa)>1) sd(aa) else NA_real_
          sb <- if(length(bb)>1) sd(bb) else NA_real_
          denom <- sqrt((sa^2+sb^2)/2)
        } else {
          va <- avg(aa==k); vb <- avg(bb==k); sa <- sb <- NA_real_
          denom <- sqrt((va*(1-va)+vb*(1-vb))/2)
        }
        comparisons[[length(comparisons)+1]] <- data.frame(population=pop,
          reference_visit='ses-02A',variable=v,category=k,
          early_N=sum(early),later_N=sum(later),
          early_observed_N=length(aa),later_observed_N=length(bb),
          early_missing_N=sum(early)-length(aa),later_missing_N=sum(later)-length(bb),
          early_value=va,later_value=vb,early_SD=sa,later_SD=sb,
          SMD=if(is.finite(denom)&&denom>0) (va-vb)/denom else NA_real_)
      }
    }
  }
  ct <- do.call(rbind,counts)
  if(any(tapply(ct$N,ct$population,sum)!=tapply(ct$total_N,ct$population,unique)))
    stop('Group counts do not partition each source population.')
  put(ct,'20_group_counts.csv')
  put(do.call(rbind,patterns),'20_account_status_patterns.csv')
  put(do.call(rbind,comparisons),'20_early_later_demographic_comparisons.csv')
  put(do.call(rbind,dictionary),'20_variable_dictionary.csv')
  put(if(length(value_labels)) do.call(rbind,value_labels) else
      data.frame(variable=character(),code=character(),label=character()),
      '20_variable_value_labels.csv')
  put(if(length(handling)) do.call(rbind,handling) else
      data.frame(variable=character(),code=character(),label=character()),
      '20_nonresponse_code_handling.csv')
  put(data.frame(input=c(basename(p),basename(ed_path)),
      MD5=unname(tools::md5sum(c(p,ed_path)))), '20_input_provenance.csv')
  writeLines(c(
    'Early: account reported at ses-02A, irrespective of later account status.',
    'Later: no account reported at ses-02A, account reported at ses-04A or ses-06A.',
    'These labels refer to observed account status, not the true date of adoption.',
    'Unknown initial status and incompletely observed non-holder trajectories are separate groups.',
    'Both source-cohort and fixed-family comparisons are exported; fixed-family is the primary descriptive contrast.',
    'Comparison values are from ses-02A only, with observed denominators; no carry-forward or imputation.',
    'Age retains the existing source units (months); categorical values are proportions.',
    'SMD = early minus later, divided by the unweighted pooled group SD.',
    'Categorical SMDs use each category versus other observed categories.',
    'Undefined SMDs remain missing; no p-values, new network models or QC exclusions are added.',
    'Future account response availability affects classification; this is a descriptive selection comparison.',
    'Review variable labels/nonresponse handling before using these summaries in the manuscript.'
  ),file.path(out,'20_definitions.txt'))
  writeLines(capture.output(sessionInfo()),file.path(out,'20_sessionInfo.txt'))
  print(ct,row.names=FALSE)
  message('20 complete. Send all 20_ files from: ',out)
})
