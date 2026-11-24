# Paper 1 revision v1.0 -- 2026-10-03
# Internal engine. Use run_analysis.R; do not source numbered engines in order.
# Set PAPER1_CODE_DIR only if the scripts are in a different directory.
if (!exists('CODE_DIR')) CODE_DIR <- normalizePath(getwd(), winslash='/', mustWork=TRUE)
CFG <- list(
  seed=20261003L, gamma=0.5, bootstrap_B=1000L, permutations=1000L,
  case_B=1000L, case_drop=c(.1,.2,.3,.4,.5,.6,.7), mice_m=20L, mice_maxit=30L,
  min_network_n=100L,
  # Provisional candidate from Kanwal's email. Source-scan correspondence still needs confirmation.
  motion_var='mr_y_qc__raw__dmri__mot__ok_mean',
  motion_source_confirmed=FALSE,
  domain_permutation=FALSE, # optional; cautious joint-network interpretation is the default
  # R1 pooled adjustment excludes visit indicators and age-by-visit terms.
  pooled_wave_fixed=FALSE,
  # Optional equivalence is disabled until a scientifically justified bound is supplied.
  equivalence_bound=NA_real_
)
RNGkind(kind='Mersenne-Twister', normal.kind='Inversion', sample.kind='Rejection')
OUT <- file.path(CODE_DIR, 'Paper1_revision_outputs')
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
need <- function(pkgs) {
  absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(absent)) stop('Install required packages first: install.packages(c(',
                         paste(sprintf('"%s"', absent), collapse=', '), '))')
}
need(c('haven','dplyr','tidyselect','lme4','qgraph','polycor','mice'))
if(packageVersion('dplyr') < '1.1.1') stop('Update dplyr to >= 1.1.1 with INSTALL_Paper1_Required_Packages.R')
SMA <- sprintf('nt_y_stq__socmed_%03d',8:13)
ROI <- c('gcpd','gcpv','gsca','gscma','gscmp','scm')
roi_vars <- function(prefix) unlist(lapply(c('lh','rh'),function(h)
  paste0(prefix,ROI,'__',h,'_mean')),use.names=FALSE)
MRI <- list(FA=roi_vars('mr_y_dti__is__fa__wm__dst__'),
            T2w=roi_vars('mr_y_smri__t2__wm__dst__'))
WAVES <- c('ses-02A','ses-04A','ses-06A')
AGE <- 'ab_g_dyn__visit_age'; FAMILY <- 'ab_g_stc__design_id__fam'
SOFTWARE <- 'ab_g_dyn__design_mr__software'; SCANNER <- 'ab_g_dyn__design_mr__serial'
GATE <- 'nt_y_stq__socmed_001'
QC <- c(FA='mr_y_qc__incl__dmri_indicator', T2w='mr_y_qc__incl__smri__t2_indicator')
MOTION <- c('mr_y_qc__raw__dmri__mot__ok_mean','mr_y_qc__raw__dmri__mot_mean')
MANUAL <- c('mr_y_qc__post__man__t2_score','mr_y_qc__post__man__dmri_score',
            'mr_y_qc__post__man__fsurf_score')
DEMOS <- c('ab_g_stc__cohort_sex','ab_g_stc__cohort_race__nih',
           'ab_g_stc__cohort_ethn','ab_g_dyn__cohort_income__hhold__3lvl',
           'ab_g_dyn__design_site')
NODES <- c('Salience','Tolerance','Mood_modification','Relapse','Withdrawal','Conflict',
           paste0(rep(c('L_','R_'),each=6),rep(ROI,2)))
csv <- function(x,name) write.csv(x,file.path(OUT,name),row.names=FALSE,na='')
read_stage <- function(name) {
  p <- file.path(OUT,name); if(!file.exists(p)) stop('Missing prerequisite: ',p)
  readRDS(p)
}
num <- function(x) as.numeric(haven::zap_labels(x))
clean_factor <- function(x) {
  v <- trimws(as.character(haven::zap_labels(x))); v[v==''] <- NA_character_; factor(v)
}
valid_num <- function(x) !is.na(x) & is.finite(x)
check_keys <- function(d,label) {
  if(anyNA(d$participant_id) || any(!nzchar(trimws(d$participant_id))) || anyNA(d$visit))
    stop(label,': missing key')
  if(anyDuplicated(d[c('participant_id','visit')])) stop(label,': duplicate participant/visit')
}
# Decode each source's labels BEFORE joining. Never assume numeric codes agree.
decode_visit <- function(d,label) {
  x <- d$session_id; raw <- haven::zap_labels(x)
  if(is.character(raw) && all(grepl('^ses-[0-9]{2}A$',raw[!is.na(raw)]))) {
    visit <- raw
  } else {
    labs <- attr(x,'labels')
    if(is.null(labs) || anyDuplicated(unname(labs)) || anyDuplicated(names(labs)))
      stop(label,': absent/ambiguous visit labels; do not infer by row order')
    visit <- names(labs)[match(raw,unname(labs))]
  }
  if(anyNA(visit) || any(!grepl('^ses-[0-9]{2}A$',visit))) stop(label,': unrecognised visit')
  mapping <- unique(data.frame(source=label,stored_code=as.character(raw),visit=visit))
  d$visit <- visit; d$participant_id <- as.character(d$participant_id)
  check_keys(d,label)
  attr(d,'visit_mapping') <- mapping
  d
}
read_sav_source <- function(path,vars=NULL) {
  # Historical paper1.sav may actually be a ZIP containing dataset.sav.
  con <- file(path,'rb'); magic <- readBin(con,'raw',n=2); close(con)
  if(identical(magic,charToRaw('PK'))) {
    zi <- utils::unzip(path,list=TRUE)
    candidates <- zi$Name[grepl('(^|/)dataset[.]sav$',zi$Name,ignore.case=TRUE)]
    if(length(candidates)!=1L) stop('ZIP must contain one identifiable dataset.sav')
    if(grepl('(^/|[.][.])',candidates)) stop('Unsafe ZIP member')
    private_tmp <- file.path(CODE_DIR, 'PRIVATE_TEMP')
    dir.create(private_tmp, recursive=TRUE, showWarnings=FALSE)
    td <- tempfile('paper1_unzip_', tmpdir=private_tmp); dir.create(td)
    on.exit(unlink(td,recursive=TRUE),add=TRUE)
    utils::unzip(path,files=candidates,exdir=td)
    path <- file.path(td,candidates)
  }
  h <- haven::read_sav(path,n_max=0)
  if(is.null(vars)) return(haven::read_sav(path))
  miss <- setdiff(vars,names(h)); if(length(miss)) stop('Missing columns: ',paste(miss,collapse=', '))
  haven::read_sav(path,col_select=tidyselect::all_of(vars))
}
flags_for <- function(d,mod,motion=CFG$motion_var,one_family=TRUE) {
  eligible <- !is.na(d[[GATE]]) & d[[GATE]]==1
  family_ok <- !is.na(d$family) & (!one_family | d$selected_family)
  qc_ok <- !is.na(d[[QC[[mod]]]]) & d[[QC[[mod]]]]==1
  cov_ok <- valid_num(d$age) & !is.na(d$software) & !is.na(d$scanner)
  mot_ok <- if(mod=='FA') valid_num(d[[motion]]) & d[[motion]]>=0 else rep(TRUE,nrow(d))
  mri_complete <- complete.cases(d[MRI[[mod]]])
  sma_complete <- complete.cases(d[SMA])
  data.frame(eligible,family_ok,qc_ok,cov_ok,mot_ok,mri_complete,sma_complete,
             base=eligible & family_ok & qc_ok & cov_ok & mot_ok,
             included=eligible & family_ok & qc_ok & cov_ok & mot_ok & mri_complete & sma_complete)
}
# Return scanner-adjusted marginal residuals: retain family/person BLUP.
# For SMA remove fixed effects only. 'conditional' is an explicit diagnostic sensitivity.
residualise <- function(d,mod,scope='pooled',family_model=TRUE,
                        extraction='retain_family',motion=CFG$motion_var) {
  vars <- c(SMA,MRI[[mod]]); ans <- matrix(NA_real_,nrow(d),length(vars),dimnames=list(NULL,NODES))
  diagnostics <- list(); models <- list()
  groups <- if(scope=='pooled') list(pooled=seq_len(nrow(d))) else split(seq_len(nrow(d)),d$visit)
  for(g in names(groups)) for(j in seq_along(vars)) {
    v <- vars[j]; ix <- groups[[g]]; z <- droplevels(d[ix,,drop=FALSE])
    use <- !is.na(z[[v]]); z <- droplevels(z[use,,drop=FALSE]); rows <- ix[use]
    if(nrow(z)<30 || sd(z[[v]])==0) stop('Insufficient variation: ',g,' ',v)
    z$y <- z[[v]]; z$motion <- z[[motion]]; z$visit_factor <- factor(z$visit)
    fixed <- 'age'
    if(scope=='pooled' && CFG$pooled_wave_fixed && nlevels(z$visit_factor)>1) fixed <- c(fixed,'visit_factor')
    if(j>6 && nlevels(z$software)>1) fixed <- c(fixed,'software')
    if(j>6 && mod=='FA') fixed <- c(fixed,'motion')
    random <- character()
    if(family_model && nlevels(z$family)>1 && nlevels(z$family)<nrow(z)) random <- c(random,'(1|family)')
    if(j>6 && nlevels(z$scanner)>1 && nlevels(z$scanner)<nrow(z)) random <- c(random,'(1|scanner)')
    f <- as.formula(paste('y ~',paste(c(fixed,random),collapse=' + ')))
    warnings <- character()
    fit <- withCallingHandlers(if(length(random)) lme4::lmer(f,z,REML=TRUE,
      control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))) else lm(f,z),
      warning=function(w){warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')})
    if(inherits(fit,'merMod')) {
      cm <- fit@optinfo$conv$lme4$messages
      if(length(cm)) stop('Model convergence problem: ',g,' ',v,': ',paste(cm,collapse='; '))
      re <- if(extraction=='conditional') NULL else if('(1|scanner)' %in% random) ~(1|scanner) else NA
      fitted_remove <- predict(fit,newdata=z,re.form=re,allow.new.levels=FALSE)
      fixed_only <- predict(fit,newdata=z,re.form=NA,allow.new.levels=FALSE)
      vc <- as.data.frame(lme4::VarCorr(fit))
      fam_v <- sum(vc$vcov[vc$grp=='family']); scan_v <- sum(vc$vcov[vc$grp=='scanner'])
      singular <- lme4::isSingular(fit)
    } else {
      fitted_remove <- fitted(fit); fixed_only <- fitted_remove
      fam_v <- scan_v <- 0; singular <- FALSE
    }
    r <- z$y-fitted_remove
    ans[rows,j] <- r
    diagnostics[[length(diagnostics)+1]] <- data.frame(modality=mod,scope=scope,group=g,
      variable=v,formula=paste(deparse(f),collapse=' '),n=nrow(z),extraction=extraction,
      family_variance=fam_v,scanner_variance=scan_v,residual_variance=sigma(fit)^2,
      extracted_sd=sd(r),conditional_sd=sd(residuals(fit)),fixed_only_residual_sd=sd(z$y-fixed_only),
      singular=singular,warnings=paste(unique(warnings),collapse='; '))
    models[[paste(g,v,sep='__')]] <- fit
  }
  list(data=ans,diagnostics=dplyr::bind_rows(diagnostics),models=models)
}
network <- function(x) {
  x <- as.matrix(x)
  if(nrow(x)<CFG$min_network_n || anyNA(x) || any(!is.finite(x))) stop('Invalid network data or N too small')
  if(any(apply(x,2,sd)<=0)) stop('Constant network node')
  R <- cor(x); if(min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)<1e-8) stop('Singular correlation matrix')
  W <- qgraph::EBICglasso(R,n=nrow(x),gamma=CFG$gamma,verbose=FALSE)
  dimnames(W) <- list(NODES,NODES); W
}
edges <- function(W) {
  ij <- which(upper.tri(W),arr.ind=TRUE)
  data.frame(i=ij[,1],j=ij[,2],node1=NODES[ij[,1]],node2=NODES[ij[,2]],
    domain=ifelse(ij[,2]<=6,'SMA-SMA',ifelse(ij[,1]>6,'MRI-MRI','SMA-MRI')),
    weight=W[ij])
}
central <- function(W) {
  mask <- outer(seq_len(18)<=6,seq_len(18)<=6,'!=')
  data.frame(node=NODES,strength=rowSums(abs(W)),expected_influence=rowSums(W),
    bridge_strength=rowSums(abs(W)*mask),bridge_expected_influence=rowSums(W*mask))
}
strengths <- function(W) {
  e <- edges(W)
  c(global=sum(abs(e$weight)),vapply(c('SMA-SMA','MRI-MRI','SMA-MRI'),
    function(g) sum(abs(e$weight[e$domain==g])),numeric(1)))
}
partial_unregularised <- function(x) {
  R <- cor(x); P <- solve(R); W <- -cov2cor(P); diag(W) <- 0; dimnames(W)<-list(NODES,NODES);W
}
compare_networks <- function(A,B,label) {
  a <- edges(A); b <- edges(B); use <- a$domain=='SMA-MRI'
  safe_cor <- function(x,y) if(sd(x)==0 || sd(y)==0) NA_real_ else cor(x,y)
  data.frame(comparison=label,edge_correlation=safe_cor(a$weight,b$weight),
    cross_edge_correlation=safe_cor(a$weight[use],b$weight[use]),
    max_abs_difference=max(abs(a$weight-b$weight)),
    mean_abs_cross_difference=mean(abs(a$weight[use]-b$weight[use])),
    global_A=sum(abs(a$weight)),global_B=sum(abs(b$weight)))
}
save_environment <- function(stage) {
  writeLines(capture.output(sessionInfo()),file.path(OUT,paste0(stage,'_sessionInfo.txt')))
  saveRDS(CFG,file.path(OUT,paste0(stage,'_settings.rds')))
}
