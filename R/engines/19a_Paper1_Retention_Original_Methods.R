# Update the three retention outcomes in the original Paper 1 analysis.
# Uses saved revised inputs and the locked family selection. Does not fit networks.
options(stringsAsFactors=FALSE)
needed <- c('haven','dplyr','brglm2','detectseparation')
absent <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(absent)) stop('Install these packages first: ',paste(absent,collapse=', '))
message('Choose 00_Paper1_Settings_and_Functions.R from your analysis scripts.')
p00 <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R') stop('Wrong settings file.')
E <- new.env(parent=globalenv()); E$CODE_DIR <- dirname(p00)
sys.source(p00,envir=E)
message('Choose 03_primary_networks_RESTRICTED.rds from Paper1_revision_outputs_original_methods_v2.')
pnet <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(pnet)!='03_primary_networks_RESTRICTED.rds') stop('Wrong network input.')
V2 <- dirname(pnet)
input_names <- c('01_merged_RESTRICTED.rds','03_settings.rds','03_primary_networks_RESTRICTED.rds','03_all_network_edges.csv')
if(!all(file.exists(file.path(V2,input_names)))) stop('The completed revised input files must be in the selected folder.')
E$CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(E$CFG$pooled_wave_fixed,FALSE) || E$CFG$gamma!=.5) stop('Wrong primary analysis configuration.')
d <- readRDS(file.path(V2,'01_merged_RESTRICTED.rds'))
nets <- readRDS(pnet)
E$check_keys(d,'Revised source')
WAVES <- E$WAVES; GATE <- E$GATE; SMA <- E$SMA; MRI <- E$MRI
if(!all(d$visit %in% WAVES)) stop('Unexpected visits.')
if(anyNA(d$selected_family)) stop('Missing family-selection flag.')
if(!isTRUE(all.equal(as.numeric(d$age),as.numeric(d[[E$AGE]]))))
 stop('Saved age differs from the original age variable; check the unit before proceeding.')
# The original Rmd used ab_g_dyn__visit_age * 12 for age in months.
if(any(d$age[is.finite(d$age)]<0)) stop('Invalid age.')
OUT <- file.path(dirname(V2),'Paper1_Retention_Results')
dir.create(OUT,showWarnings=FALSE)
PRIVATE <- file.path(OUT,'RESTRICTED'); dir.create(PRIVATE,showWarnings=FALSE)
SHARE <- file.path(OUT,'SHAREABLE_RETENTION'); dir.create(SHARE,showWarnings=FALSE)
if(file.exists(file.path(SHARE,'19a_complete.csv'))) unlink(file.path(SHARE,'19a_complete.csv'))
write_share <- function(x,n) write.csv(x,file.path(SHARE,n),row.names=FALSE,na='')
provenance <- data.frame(file=input_names,md5=unname(tools::md5sum(file.path(V2,input_names))))
write_share(provenance,'19a_input_provenance.csv')

# Validate complete-case membership against the saved primary networks.
fFA <- E$flags_for(d,'FA'); fT2 <- E$flags_for(d,'T2w')
for(mod in c('FA','T2w')) for(w in WAVES) {
 f <- if(mod=='FA') fFA else fT2
 ids <- d$participant_id[d$visit==w & f$included]
 key <- paste(mod,w,sep='_')
 if(is.null(nets[[key]]) || !setequal(ids,nets[[key]]$ids))
  stop('Reconstructed sample differs from the saved primary network: ',key)
}
# Participant population preserves the original union across target visits,
# restricted to the locked family sample and account holders at >=1 target visit.
# This is sustained eligible data availability, not pure study dropout.
ids <- sort(unique(d$participant_id[d$selected_family & d[[GATE]] %in% 1]))
if(!length(ids)) stop('No eligible participants.')
x <- d[d$participant_id %in% ids,,drop=FALSE]
x <- x[order(x$participant_id,match(x$visit,WAVES)),]
rows <- split(seq_len(nrow(x)),x$participant_id)
status <- data.frame(participant_id=d$participant_id,visit=d$visit,
 any_network=fFA$base & rowSums(!is.na(d[SMA]))>0 & rowSums(!is.na(d[MRI$FA]))>0 |
             fT2$base & rowSums(!is.na(d[SMA]))>0 & rowSums(!is.na(d[MRI$T2w]))>0,
 complete_FA=fFA$included,complete_T2w=fT2$included)
lookup <- paste(status$participant_id,status$visit,sep='|')
first_observed <- function(z,col) {
 order_rows <- c(which(z$visit==WAVES[1]),which(z$visit!=WAVES[1]))
 values <- z[[col]][order_rows]; take <- which(!is.na(values))[1]
 if(is.na(take)) return(list(value=NA,visit=NA_character_))
 raw <- haven::zap_labels(values[take]); if(is.factor(raw)) raw <- as.character(raw)
 value <- suppressWarnings(as.numeric(raw))
 if(is.na(value)) stop('Non-numeric original predictor coding: ',col)
 list(value=value,visit=z$visit[order_rows[take]])
}
variables <- c(baseline_age_months='age',sex='ab_g_stc__cohort_sex',site='ab_g_dyn__design_site',
 race='ab_g_stc__cohort_race__nih',ethnicity='ab_g_stc__cohort_ethn',income='ab_g_dyn__cohort_income__hhold__3lvl')
if(!all(variables %in% names(x))) stop('Required original predictor missing.')
people <- vector('list',length(rows)); predictor_source <- list(); person_states <- list()
for(i in seq_along(rows)) {
 z <- x[rows[[i]],,drop=FALSE]; id <- z$participant_id[1]
 rec <- list(participant_id=id)
 for(v in names(variables)) {
  got <- first_observed(z,variables[[v]])
  rec[[v]] <- got$value * if(v=='baseline_age_months') 12 else 1
  predictor_source[[length(predictor_source)+1L]] <- data.frame(variable=v,source_visit=if(is.na(got$visit)) 'missing' else got$visit)
 }
 for(j in seq_along(WAVES)) {
  w <- WAVES[j]; r <- match(paste(id,w,sep='|'),lookup); rz <- which(z$visit==w)
  rec[[paste0('present_',j)]] <- length(rz)==1L
  gate <- if(length(rz)) as.numeric(z[[GATE]][rz]) else NA_real_
  rec[[paste0('state_',j)]] <- if(!length(rz)) 'visit_absent' else if(is.na(gate)) 'account_unknown' else if(gate==0) 'no_account' else if(gate==1) 'account_holder' else stop('Invalid account coding.')
  for(outcome in c('any_network','complete_FA','complete_T2w'))
   rec[[paste0(outcome,'_',j)]] <- if(is.na(r)) FALSE else status[[outcome]][r]
 }
 for(outcome in c('any_network','complete_FA','complete_T2w'))
  rec[[paste0('retained_',outcome)]] <- as.integer(all(unlist(rec[paste0(outcome,'_',1:3)])))
 people[[i]] <- as.data.frame(rec)
}
person <- dplyr::bind_rows(people)
if(anyDuplicated(person$participant_id)) stop('Duplicate person records.')
for(v in setdiff(names(variables),'baseline_age_months')) person[[v]] <- factor(person[[v]])
saveRDS(person,file.path(PRIVATE,'19a_retention_person_level_RESTRICTED.rds'))
source_tab <- as.data.frame(table(dplyr::bind_rows(predictor_source)))
names(source_tab)[ncol(source_tab)] <- 'N';source_tab <- source_tab[source_tab$N>0,]
write_share(source_tab,'19a_predictor_source_visits.csv')

wave_counts <- lapply(seq_along(WAVES),function(j) {
 z <- d[d$visit==WAVES[j],]; state <- person[[paste0('state_',j)]]
 data.frame(visit=WAVES[j],source_records=nrow(z),source_account_holders=sum(z[[GATE]] %in% 1),
 retention_population=nrow(person),visit_present=sum(person[[paste0('present_',j)]]),
 account_holders=sum(state=='account_holder'),no_account=sum(state=='no_account'),
 account_unknown=sum(state=='account_unknown'),visit_absent=sum(state=='visit_absent'),
 any_network_data=sum(person[[paste0('any_network_',j)]]),
 complete_FA=sum(person[[paste0('complete_FA_',j)]]),complete_T2w=sum(person[[paste0('complete_T2w_',j)]]))
})
write_share(dplyr::bind_rows(wave_counts),'19a_retention_wave_availability.csv')
patterns <- data.frame(participation=apply(person[paste0('present_',1:3)],1,function(a)paste(as.integer(a),collapse='')),
 eligibility=apply(person[paste0('state_',1:3)],1,paste,collapse=' | '))
patterns <- as.data.frame(table(patterns));names(patterns)[3]<-'N';patterns<-patterns[patterns$N>0,]
write_share(patterns,'19a_retention_patterns.csv')
outcomes <- c(retained_any_network='Eligible SMA and at least one imaging domain at all three visits',
 retained_complete_FA='Complete eligible FA network data at all three visits',
 retained_complete_T2w='Complete eligible T2-weighted network data at all three visits')
summary_counts <- data.frame(outcome=names(outcomes),definition=unname(outcomes),N_population=nrow(person),
 N_retained=vapply(names(outcomes),function(o)sum(person[[o]]),numeric(1)))
summary_counts$N_not_retained <- summary_counts$N_population-summary_counts$N_retained
write_share(summary_counts,'19a_retention_outcome_counts.csv')
# Original descriptive comparison, with observed denominators made explicit.
descriptions <- list()
for(o in names(outcomes)) for(g in 0:1) for(v in names(variables)) {
 a <- person[[v]][person[[o]]==g]; obs <- a[!is.na(a)]
 if(v=='baseline_age_months') {
  descriptions[[length(descriptions)+1]] <- data.frame(outcome=o,retained=g,variable=v,category='continuous',
   N_group=length(a),N_observed=length(obs),N_category=NA,percent=NA,
   mean=if(length(obs))mean(obs) else NA,SD=if(length(obs)>1)sd(obs) else NA)
 } else for(lev in levels(person[[v]])) {
  descriptions[[length(descriptions)+1]] <- data.frame(outcome=o,retained=g,variable=v,category=lev,
   N_group=length(a),N_observed=length(obs),N_category=sum(obs==lev),
   percent=if(length(obs))100*mean(obs==lev) else NA,mean=NA,SD=NA)
 }
}
write_share(dplyr::bind_rows(descriptions),'19a_retained_vs_not_retained_descriptives.csv')

model_summaries <- coefficient_rows <- fit_attempts <- references <- list()
for(o in names(outcomes)) {
 dat <- person[c(o,names(variables))];names(dat)[1] <- 'y'
 dat <- droplevels(dat[complete.cases(dat),])
 if(nrow(dat)<30 || length(unique(dat$y))!=2L) stop('Insufficient outcome variation: ',o)
 usable <- names(variables)[vapply(dat[names(variables)],function(a)length(unique(a))>1,logical(1))]
 form <- reformulate(usable,response='y')
 # Centering improves numerical conditioning; age slope remains per month.
 age_center <- mean(dat$baseline_age_months);dat$baseline_age_months <- dat$baseline_age_months-age_center
 X <- model.matrix(form,dat)
 if(qr(X)$rank!=ncol(X)) stop('Rank-deficient retention model: ',o,'. No categories were silently removed.')
 for(v in usable) if(is.factor(dat[[v]])) references[[length(references)+1]] <-
  data.frame(outcome=o,variable=v,reference_category=levels(dat[[v]])[1])
 sep <- glm(form,data=dat,family=binomial(),method=detectseparation::detect_separation)
 if(length(coef(sep))!=ncol(X) || anyNA(coef(sep))) stop('Separation diagnosis failed: ',o)
 separated <- any(is.infinite(coef(sep)))
 fit <- NULL; method <- 'Maximum likelihood'; warns <- character()
 capture_fit <- function(expr) tryCatch(withCallingHandlers(expr,
  warning=function(w){warns <<- c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e){warns <<- c(warns,conditionMessage(e));NULL})
 adequate <- function(f) {
  if(is.null(f) || !isTRUE(f$converged)) return(FALSE)
  b <- coef(f);se <- tryCatch(sqrt(diag(vcov(f))),error=function(e)rep(NA_real_,length(b)))
  length(b)==ncol(X) && all(is.finite(b)) && all(is.finite(se) & se>0) &&
   all(is.finite(exp(b))) && all(is.finite(exp(b+1.96*se)))
 }
 if(!separated) fit <- capture_fit(glm(form,data=dat,family=binomial(),control=glm.control(maxit=200)))
 fit_attempts[[length(fit_attempts)+1]] <- data.frame(outcome=o,method='Maximum likelihood',
  attempted=!separated,separation=separated,converged=adequate(fit),warnings=paste(warns,collapse='; '))
 if(!adequate(fit)) {
  method <- 'Mean bias-reduced logistic regression'
  starts <- setNames(rep(0,ncol(X)),colnames(X));starts[1] <- qlogis(mean(dat$y))
  for(j in 1:3) {
   warns <- character()
   fit <- capture_fit(glm(form,data=dat,family=binomial(),start=starts,method=brglm2::brglmFit,type='AS_mean',
    control=list(maxit=c(1000L,2000L,4000L)[j],epsilon=1e-8,slowit=c(.5,.2,.05)[j],max_step_factor=20)))
   fit_attempts[[length(fit_attempts)+1]] <- data.frame(outcome=o,method=method,attempted=TRUE,
    separation=separated,converged=adequate(fit),warnings=paste(warns,collapse='; '))
   if(adequate(fit)) break
  }
 }
 write_share(dplyr::bind_rows(fit_attempts),'19a_model_fit_attempts.csv')
 if(!adequate(fit)) stop('Retention model failed numerical checks: ',o,'. Send 19a_model_fit_attempts.csv.')
 b <- coef(fit);se <- sqrt(diag(vcov(fit)));pv <- 2*pnorm(abs(b/se),lower.tail=FALSE)
 coefficient_rows[[o]] <- data.frame(outcome=o,term=names(b),beta=unname(b),SE=unname(se),
  OR=exp(unname(b)),lower=exp(unname(b-1.96*se)),upper=exp(unname(b+1.96*se)),p=unname(pv),
  method=method,interval='Approximate Wald',row.names=NULL)
 model_summaries[[o]] <- data.frame(outcome=o,definition=outcomes[[o]],N_population=nrow(person),
  N_model=nrow(dat),N_retained=sum(dat$y),N_not_retained=sum(dat$y==0),
  N_missing_predictors=nrow(person)-nrow(dat),parameters=ncol(X),
  predictors=paste(usable,collapse=' + '),method=method,separation=separated,converged=fit$converged,
  age_center_months=age_center,AIC=if(method=='Maximum likelihood')AIC(fit) else NA_real_,
  null_deviance=if(method=='Maximum likelihood')fit$null.deviance else NA_real_,
  residual_deviance=if(method=='Maximum likelihood')fit$deviance else NA_real_)
 saveRDS(fit,file.path(PRIVATE,paste0('19a_',o,'_model_RESTRICTED.rds')))
}
co <- dplyr::bind_rows(coefficient_rows)
co$p_BH_retention <- co$p_Holm_retention <- NA_real_
ix <- co$term!='(Intercept)'
co$p_BH_retention[ix] <- p.adjust(co$p[ix],'BH');co$p_Holm_retention[ix] <- p.adjust(co$p[ix],'holm')
write_share(co,'19a_retention_logistic_coefficients.csv')
write_share(dplyr::bind_rows(model_summaries),'19a_retention_model_summary.csv')
write_share(dplyr::bind_rows(references),'19a_retention_reference_categories.csv')
write_share(data.frame(component=c('Population','Outcomes','Eligibility','Complete cases','Predictors','Age','Missing predictors','Estimation','Multiplicity','Interpretation'),
 definition=c('Locked one-participant-per-family sample; account holder at one or more of the three target visits. This retains the original union-of-visits design while excluding never-eligible records.',
 'The original three outcomes: any SMA-imaging data, complete FA data and complete T2-weighted data at all three visits.',
 'Each visit requires account ownership and modality-specific QC. Missing account status is not treated as eligible.',
 'Complete-case visit membership is checked against the saved primary network participant IDs.',
 'Age, sex, site, race, ethnicity and income; ses-02A value when observed, otherwise first observed target-visit value, following the original script.',
 'Age in months, as in the original retention script; centered within the model sample. The age odds ratio is per month.',
 'Complete-case predictor analysis; no predictor imputation. Source visits and missing counts are reported.',
 'Original maximum-likelihood logistic regression when estimable; mean bias reduction only for separation or failed numerical checks. Method is identified for every model.',
 'Raw p-values plus BH and Holm across all non-intercept coefficients of the three retention models, as one reporting family.',
 'Sustained eligible analytic-data availability combines visit participation, account status, imaging quality and completeness; it is not a measure of study dropout alone.')),
 '19a_retention_definitions.csv')
writeLines(capture.output(sessionInfo()),file.path(SHARE,'19a_sessionInfo.txt'))
files_final <- list.files(SHARE,pattern='[.]csv$',full.names=TRUE)
write_share(data.frame(file=basename(files_final),md5=unname(tools::md5sum(files_final))), '19a_complete.csv')
message('19a complete. Aggregated outputs: ',normalizePath(SHARE,winslash='/'))
message('Next run 19b; select 19a_retention_model_summary.csv when requested. Do not share RESTRICTED files.')
