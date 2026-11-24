source('00_Paper1_Settings_and_Functions.R')
d <- read_stage('01_merged_RESTRICTED.rds'); cats <- read_stage('01_sma_categories.rds')
flow <- reasons <- dist <- freq <- selection <- motion_audit <- list()
for(w in WAVES) {
 z <- d[d$visit==w,,drop=FALSE]
 for(v in SMA) {
  x <- z[[v]][z[[GATE]] %in% 1]; observed <- x[!is.na(x)]
  sk <- if(length(observed)>2 && sd(observed)>0) mean((observed-mean(observed))^3)/mean((observed-mean(observed))^2)^1.5 else NA_real_
  dist[[length(dist)+1]] <- data.frame(visit=w,item=v,eligible=length(x),observed=length(observed),
    missing=sum(is.na(x)),mean=mean(observed),sd=sd(observed),skewness_moment=sk,
    floor_percent=100*mean(observed==min(cats[[v]])),ceiling_percent=100*mean(observed==max(cats[[v]])))
  freq[[length(freq)+1]] <- data.frame(visit=w,item=v,category=cats[[v]],
    count=vapply(cats[[v]],function(k)sum(observed==k),numeric(1)),denominator_observed=length(observed))
 }
 for(mod in names(MRI)) {
  f <- flags_for(z,mod); stage <- rep(TRUE,nrow(z)); counts <- c(source=nrow(z))
  reason <- rep('included',nrow(z))
  checks <- c('eligible','family_ok','qc_ok','cov_ok','mot_ok','mri_complete','sma_complete')
  for(k in checks) {
   newly <- stage & !f[[k]]; reason[newly] <- paste0('excluded_',k)
   stage <- stage & f[[k]]; counts[k] <- sum(stage)
  }
  flow[[length(flow)+1]] <- data.frame(visit=w,modality=mod,step=names(counts),N=as.integer(counts))
  reasons[[length(reasons)+1]] <- data.frame(participant_id=z$participant_id,visit=w,modality=mod,
    first_exclusion_reason=reason,account_no=z[[GATE]] %in% 0,account_unknown=is.na(z[[GATE]]),
    qc_record_absent=!(z$in_mri_extract %in% TRUE),qc_value_missing=is.na(z[[QC[[mod]]]]),
    qc_failed=z[[QC[[mod]]]] %in% 0,missing_family=is.na(z$family),
    eligible_item_nonresponse=f$eligible & !f$sma_complete)
  # Selection contrasts within account holders; also within the locked family-selected group.
  for(pop in c('all_account_holders','selected_family_account_holders')) {
   idx <- f$eligible & if(pop=='selected_family_account_holders') f$family_ok else TRUE
   zz <- z[idx,,drop=FALSE]; inc <- f$included[idx]
   for(v in c('age',DEMOS,SMA,MRI[[mod]],MOTION)) {
    if(v %in% DEMOS) {
      xx <- as.character(haven::zap_labels(zz[[v]])); lev <- sort(unique(xx[!is.na(xx)]))
      for(k in lev) {
       a<-mean(xx[inc]==k,na.rm=TRUE); b<-mean(xx[!inc]==k,na.rm=TRUE)
       denom<-sqrt((a*(1-a)+b*(1-b))/2)
       selection[[length(selection)+1]]<-data.frame(visit=w,modality=mod,population=pop,variable=v,category=k,
        included_n=sum(inc),excluded_n=sum(!inc),observed_included=sum(!is.na(xx[inc])),
        observed_excluded=sum(!is.na(xx[!inc])),included_value=a,excluded_value=b,
        SMD=if(is.finite(denom)&&denom>0)(a-b)/denom else NA_real_)
      }
    } else {
      xx<-num(zz[[v]]); a<-xx[inc]; b<-xx[!inc]; denom<-sqrt((var(a,na.rm=TRUE)+var(b,na.rm=TRUE))/2)
      selection[[length(selection)+1]]<-data.frame(visit=w,modality=mod,population=pop,variable=v,category='continuous',
       included_n=sum(inc),excluded_n=sum(!inc),observed_included=sum(!is.na(a)),observed_excluded=sum(!is.na(b)),
       included_value=mean(a,na.rm=TRUE),excluded_value=mean(b,na.rm=TRUE),
       SMD=if(is.finite(denom)&&denom>0)(mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE))/denom else NA_real_)
    }
   }
  }
 }
 for(pop in c('all_dmri_pass','eligible_selected_dmri_pass')) {
  take <- z[[QC[['FA']]]] %in% 1
  if(pop=='eligible_selected_dmri_pass')take<-take & z[[GATE]] %in% 1 & z$selected_family
  a<-z[[MOTION[1]]][take]; b<-z[[MOTION[2]]][take]; both<-is.finite(a)&is.finite(b)
  motion_audit[[length(motion_audit)+1]]<-data.frame(visit=w,population=pop,N=sum(take),
   missing_ok=sum(is.na(a)),missing_all=sum(is.na(b)),negative_ok=sum(a<0,na.rm=TRUE),
   zero_ok=sum(a==0,na.rm=TRUE),max_ok=if(any(is.finite(a)))max(a[is.finite(a)])else NA,
   correlation=if(sum(both)>2 && sd(a[both])>0 && sd(b[both])>0)cor(a[both],b[both])else NA,
   mean_absolute_difference=if(any(both))mean(abs(a[both]-b[both]))else NA)
 }
}
csv(dplyr::bind_rows(flow),'02_sequential_sample_flow.csv')
csv(dplyr::bind_rows(reasons),'02_record_exclusion_reasons_RESTRICTED.csv')
csv(dplyr::bind_rows(dist),'02_SMA_distributions_eligible.csv')
csv(dplyr::bind_rows(freq),'02_SMA_category_frequencies.csv')
csv(dplyr::bind_rows(selection),'02_included_vs_excluded_SMD.csv')
csv(dplyr::bind_rows(motion_audit),'02_motion_audit.csv')
# Structural non-administration and eligible-item missingness explicitly separated.
missingness <- list(); patterns <- list(); logit_results <- list()
for(w in WAVES) {
 z<-d[d$visit==w,]; eligible<-z[[GATE]] %in% 1
 for(status in c('yes','no','unknown')) {
  take<-switch(status,yes=eligible,no=z[[GATE]] %in% 0,unknown=is.na(z[[GATE]]))
  nobs<-rowSums(!is.na(z[take,SMA,drop=FALSE]))
  missingness[[length(missingness)+1]]<-data.frame(visit=w,account=status,N=sum(take),
   six_items=sum(nobs==6),partial_items=sum(nobs>0 & nobs<6),no_items=sum(nobs==0))
 }
 e<-z[eligible & z$selected_family,,drop=FALSE]
 p<-apply(is.na(e[SMA]),1,paste0,collapse='')
 patterns[[w]]<-data.frame(visit=w,pattern=names(table(p)),N=as.integer(table(p)))
 # Selection/missingness models: descriptive, single-family sample within each wave.
 for(target in c('SMA_nonresponse','FA_excluded','T2w_excluded')) {
  y<-switch(target,SMA_nonresponse=!complete.cases(e[SMA]),FA_excluded=!flags_for(e,'FA')$included,
            T2w_excluded=!flags_for(e,'T2w')$included)
  zmod<-data.frame(y=as.integer(y),age=e$age)
  for(v in DEMOS)zmod[[v]]<-clean_factor(e[[v]])
  zmod<-droplevels(zmod[complete.cases(zmod),])
  usable<-names(zmod)[-1];usable<-usable[vapply(zmod[usable],function(x)length(unique(x))>1,logical(1))]
  if(nrow(zmod)<30 || length(unique(zmod$y))<2)next
  fit<-glm(reformulate(usable,response='y'),data=zmod,family=binomial())
  co<-coef(summary(fit));
  logit_results[[length(logit_results)+1]]<-data.frame(visit=w,target=target,N=nobs(fit),
    term=rownames(co),beta=co[,1],SE=co[,2],OR=exp(co[,1]),
    lower=exp(co[,1]-1.96*co[,2]),upper=exp(co[,1]+1.96*co[,2]),p=co[,4],converged=fit$converged,
    separation_warning=any(abs(co[,1])>10 | co[,2]>10))
 }
}
csv(dplyr::bind_rows(missingness),'02_structural_vs_item_missingness.csv')
csv(dplyr::bind_rows(patterns),'02_eligible_missing_patterns.csv')
lg<-dplyr::bind_rows(logit_results)
if(nrow(lg)){lg$p_BH_all_coefficients<-p.adjust(lg$p,'BH');csv(lg,'02_missingness_selection_logistic.csv')}
# Availability transitions are NOT study withdrawal; gate transitions explicitly tabulated.
trans<-list()
for(pair in combn(WAVES,2,simplify=FALSE)) {
 a<-d[d$visit==pair[1],c('participant_id',GATE)]; b<-d[d$visit==pair[2],c('participant_id',GATE)]
 a$visit_present_A<-TRUE;b$visit_present_B<-TRUE
 t<-merge(a,b,by='participant_id',all=TRUE,suffixes=c('_A','_B'))
 status<-function(x,present)ifelse(is.na(present),'no_source_visit',ifelse(is.na(x),'gate_missing',ifelse(x==1,'account_yes','account_no')))
 tab<-as.data.frame(table(status(t[[paste0(GATE,'_A')]],t$visit_present_A),
                          status(t[[paste0(GATE,'_B')]],t$visit_present_B)))
 names(tab)<-c('status_A','status_B','N');tab$visit_A<-pair[1];tab$visit_B<-pair[2];trans[[length(trans)+1]]<-tab
}
csv(dplyr::bind_rows(trans),'02_account_and_visit_transitions.csv')
# Supporting manual-QC audit: missing manual review is NOT a failure.
manual_audit <- list()
for(w in WAVES) for(mod in names(MRI)) {
 z <- d[d$visit==w,,drop=FALSE]
 fields <- c(if(mod=='FA') MANUAL[2] else MANUAL[1], MANUAL[3])
 for(v in fields) {
  tab <- as.data.frame(table(inclusion=z[[QC[[mod]]]],manual=z[[v]],useNA='ifany'))
  tab$visit <- w; tab$modality <- mod; tab$manual_field <- v
  manual_audit[[length(manual_audit)+1]] <- tab
 }
}
csv(dplyr::bind_rows(manual_audit),'02_manual_QC_supporting_audit.csv')

save_environment('02');message('02 complete. Review flow, gate inconsistencies and motion audit before interpreting models.')
