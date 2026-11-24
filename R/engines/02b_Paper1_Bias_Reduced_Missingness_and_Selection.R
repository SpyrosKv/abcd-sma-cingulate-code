# Run in the ORIGINAL Paper1 RStudio Project after 01-08.
# Additional diagnostic models; does not overwrite 01-07 or select a new sample.
source('00_Paper1_Settings_and_Functions.R')
need('brglm2') # install.packages('brglm2') once if missing
d <- read_stage('01_merged_RESTRICTED.rds')
results <- audits <- counts <- attempts <- list()
for(w in WAVES) {
 z <- d[d$visit==w & d[[GATE]] %in% 1 & d$selected_family,,drop=FALSE]
 for(target in c('SMA_nonresponse','FA_excluded','T2w_excluded')) {
  y <- switch(target,SMA_nonresponse=!complete.cases(z[SMA]),
              FA_excluded=!flags_for(z,'FA')$included,
              T2w_excluded=!flags_for(z,'T2w')$included)
  model_data <- data.frame(y=as.integer(y),age=z$age)
  for(v in DEMOS) model_data[[v]] <- clean_factor(z[[v]])
  # Same demographic complete-case rule as 02; no invented values/categories.
  dat <- droplevels(model_data[complete.cases(model_data),,drop=FALSE])
  usable <- setdiff(names(dat),'y')
  usable <- usable[vapply(dat[usable],function(x)length(unique(x))>1,logical(1))]
  for(v in usable) if(is.factor(dat[[v]])) {
   tab <- table(dat[[v]],factor(dat$y,levels=0:1))
   counts[[length(counts)+1]] <- data.frame(visit=w,target=target,variable=v,
    category=rownames(tab),non_events=tab[,1],events=tab[,2],
    zero_outcome_cell=tab[,1]==0 | tab[,2]==0)
  }
  form <- reformulate(usable,response='y')
  if(nrow(dat)<30 || length(unique(dat$y))<2) {
   audits[[length(audits)+1]] <- data.frame(visit=w,target=target,N_eligible=nrow(z),
    N_model=nrow(dat),events=sum(dat$y),parameters=NA_integer_,rank=NA_integer_,
    converged=FALSE,status='not_estimable',warnings='Insufficient observations/outcome variation')
   next
  }
  # Center age for numerical stability; the age OR remains per one year.
  age_center <- mean(dat$age)
  dat$age <- dat$age - age_center
  X <- model.matrix(form,dat)
  if(qr(X)$rank<ncol(X)) stop('Rank deficient design: ',w,' ',target,
    '. Do not silently delete demographic categories.')
  # Avoid adjusted-data ML starting values in sparse/separated designs.
  start_safe <- setNames(rep(0,ncol(X)),colnames(X))
  start_safe['(Intercept)'] <- qlogis(mean(dat$y))
  ok <- FALSE
  for(attempt in seq_len(3)) {
   warns <- character()
   slow <- c(0.5,0.2,0.05)[attempt]
   limit <- c(1000L,2000L,4000L)[attempt]
   fit <- tryCatch(withCallingHandlers(glm(form,data=dat,family=binomial(),
     start=start_safe,method=brglm2::brglmFit,type='AS_mean',
     control=list(maxit=limit,epsilon=1e-8,slowit=slow,max_step_factor=20)),
     warning=function(cnd){warns <<- c(warns,conditionMessage(cnd));invokeRestart('muffleWarning')}),
     error=function(e){warns <<- c(warns,conditionMessage(e));NULL})
   beta <- if(!is.null(fit)) coef(fit) else rep(NA_real_,ncol(X))
   se <- if(!is.null(fit)) tryCatch(sqrt(diag(vcov(fit))),
     error=function(e)rep(NA_real_,ncol(X))) else rep(NA_real_,ncol(X))
   ok <- !is.null(fit) && isTRUE(fit$converged) &&
     all(is.finite(beta)) && all(is.finite(se)) && all(se>0)
   attempts[[length(attempts)+1]] <- data.frame(visit=w,target=target,
     attempt=attempt,slowit=slow,maxit=limit,iterations=if(is.null(fit)) NA_integer_ else fit$iter,
     converged=ok,age_center_years=age_center,
     brglm2_version=as.character(packageVersion('brglm2')),
     warnings=paste(unique(warns),collapse='; '))
   csv(dplyr::bind_rows(attempts),'02b_numerical_fit_attempts.csv')
   if(ok) break
  }
  if(!ok) {
   csv(dplyr::bind_rows(counts),'02b_category_outcome_cells.csv')
   stop('Bias-reduced model still failed: ',w,' ',target,
     '. Send 02b_numerical_fit_attempts.csv and 02b_category_outcome_cells.csv; do not use failed estimates.')
  }
  p <- 2*pnorm(abs(beta/se),lower.tail=FALSE)
  results[[length(results)+1]] <- data.frame(visit=w,target=target,N=nrow(dat),
    events=sum(dat$y),term=names(beta),beta=unname(beta),SE=unname(se),
    OR=exp(beta),lower=exp(beta-1.96*se),upper=exp(beta+1.96*se),p=unname(p),
    method='mean_bias_reduced_AS_mean',interval='approximate_Wald',
    sparse_event_model=sum(dat$y)/ncol(X)<10,converged=ok,row.names=NULL)
  audits[[length(audits)+1]] <- data.frame(visit=w,target=target,N_eligible=nrow(z),
    N_model=nrow(dat),events=sum(dat$y),parameters=ncol(X),rank=qr(X)$rank,
    converged=ok,status='estimated',warnings=paste(unique(warns),collapse='; '))
 }
}
out <- dplyr::bind_rows(results)
if(nrow(out)) {
 out$p_Holm_all_nonintercept <- out$p_BH_all_nonintercept <- NA_real_
 ix <- out$term!='(Intercept)'
 out$p_Holm_all_nonintercept[ix] <- p.adjust(out$p[ix],'holm')
 out$p_BH_all_nonintercept[ix] <- p.adjust(out$p[ix],'BH')
}
csv(out,'02b_bias_reduced_logistic.csv')
csv(dplyr::bind_rows(audits),'02b_model_sample_and_convergence.csv')
csv(dplyr::bind_rows(counts),'02b_category_outcome_cells.csv')
writeLines(c('Supersedes naive OR inference from 02_missingness_selection_logistic.csv.',
 'Same participants, predictors and outcome definitions as 02. Mean bias reduction addresses separation.',
 'Age centered at the model-sample mean; age OR is still per year. Intercept refers to centered age.',
 'Stable prevalence-based start and smaller optimization steps; all attempts recorded. Nonconvergence stops export of final estimates.',
 'Finite estimates do not eliminate sparse-data uncertainty. Intervals/p-values are approximate Wald.',
 'Missing demographic covariates still exclude observations. Account ownership=0 is never modelled as SMA nonresponse.',
 'Selection combines absent imaging, failed/missing QC and incomplete measures; it is not a pure QC-failure model.',
 'These descriptive models neither establish MAR/MCAR nor correct selection bias in the networks.'),
 file.path(OUT,'02b_interpretation_notes.txt'))
save_environment('02b')
message('02b complete. Review sparse-event flags and uncertainty; finite ORs are not proof of an adequate missingness model.')
