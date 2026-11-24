# Supplementary revision analyses requested/proposed by VAS.
# Does not replace or overwrite the completed primary v2 analysis.
message('Choose the existing 00_Paper1_Settings_and_Functions.R.')
p00 <- normalizePath(file.choose(), winslash='/', mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R') stop('Wrong settings file.')
CODE_DIR <- dirname(p00)
source(p00, local=TRUE)
need(c('NetworkComparisonTest','bootnet'))
V2 <- file.path(CODE_DIR,'Paper1_revision_outputs_original_methods_v2')
required <- c('03_settings.rds','03_primary_networks_RESTRICTED.rds',
              '03_primary_residuals_RESTRICTED.rds','06_sensitivity_edges.csv')
if(!all(file.exists(file.path(V2,required)))) stop('Completed v2 inputs missing.')
CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(CFG$pooled_wave_fixed,FALSE) || CFG$gamma!=.5)
  stop('Unexpected primary configuration.')
CFG$permutations <- 1000L
CFG$bootstrap_B <- CFG$case_B <- 1000L
OUT <- file.path(CODE_DIR,'Paper1_VAS_remaining_checks')
dir.create(OUT,showWarnings=FALSE)
private <- file.path(OUT,'RESTRICTED')
share <- file.path(OUT,'SHAREABLE_VAS_CHECKS')
dir.create(private,showWarnings=FALSE);dir.create(share,showWarnings=FALSE)
# Resume only with identical source files and the same script revision.
signature <- list(version='13-v1',source_md5=tools::md5sum(file.path(V2,required)))
sigfile <- file.path(private,'13_input_signature.rds')
if(file.exists(sigfile)) {
  if(!identical(readRDS(sigfile),signature)) stop('Inputs changed; use a new output folder.')
} else saveRDS(signature,sigfile)
nets <- readRDS(file.path(V2,'03_primary_networks_RESTRICTED.rds'))
res <- readRDS(file.path(V2,'03_primary_residuals_RESTRICTED.rds'))
if(length(nets)!=6L) stop('Expected six primary networks.')

# A. All 432 UNREGULARISED cross-domain partial correlations.
# Conditional on the saved preprocessing; approximate Gaussian Fisher-z inference.
# Conditioning set = other 16 network nodes. SE = 1/sqrt(N-16-3).
# These are NOT significance/equivalence tests of lasso estimates.
all_pairs <- matrix_audit <- list()
for(key in names(nets)) {
  p <- nets[[key]];x <- as.matrix(p$x)
  if(nrow(x)!=p$n || anyNA(x) || any(!is.finite(x)) || ncol(x)!=18L)
    stop('Invalid primary input: ',key)
  R <- cor(x);ev <- eigen(R,symmetric=TRUE,only.values=TRUE)$values
  if(min(ev)<=1e-8) stop('Cannot invert correlation matrix: ',key)
  P <- solve(R);W <- -cov2cor(P);diag(W) <- 0
  dimnames(W) <- list(NODES,NODES)
  e <- edges(W);e <- e[e$domain=='SMA-MRI',,drop=FALSE]
  if(nrow(e)!=72L || any(abs(e$weight)>=1)) stop('Invalid partial correlations.')
  # Cross-check precision-matrix estimates using pairwise residual regression.
  errors <- numeric(72)
  for(k in seq_len(72)) {
    others <- setdiff(seq_len(18),c(e$i[k],e$j[k]))
    X <- cbind(1,x[,others,drop=FALSE])
    a <- lm.fit(X,x[,e$i[k]])$residuals
    b <- lm.fit(X,x[,e$j[k]])$residuals
    errors[k] <- abs(cor(a,b)-e$weight[k])
  }
  if(max(errors)>1e-8) stop('Partial-correlation cross-check failed.')
  e$network <- key;e$N <- nrow(x);e$conditioning_nodes <- 16L
  e$partial_r <- e$weight;e$weight <- NULL
  e$fisher_z <- atanh(e$partial_r);e$SE_z <- 1/sqrt(e$N-19)
  e$lower95 <- tanh(e$fisher_z-qnorm(.975)*e$SE_z)
  e$upper95 <- tanh(e$fisher_z+qnorm(.975)*e$SE_z)
  e$lower90 <- tanh(e$fisher_z-qnorm(.95)*e$SE_z)
  e$upper90 <- tanh(e$fisher_z+qnorm(.95)*e$SE_z)
  e$scope <- 'unregularised; conditional on saved preprocessing; approximate Gaussian inference'
  all_pairs[[key]] <- e
  matrix_audit[[key]] <- data.frame(network=key,N=nrow(x),min_eigenvalue=min(ev),
    condition_number=max(ev)/min(ev),maximum_regression_check_error=max(errors))
}
est <- dplyr::bind_rows(all_pairs)
if(nrow(est)!=432L) stop('Expected 432 cross-domain estimates.')
csv(est,'13_unregularised_432_partial_correlations.csv')
csv(dplyr::bind_rows(matrix_audit),'13_unregularised_matrix_audit.csv')

# VAS proposed .05, with stricter .03 sensitivity. Neither is a validated clinical cutoff.
# Both bounds are reported. They were selected during revision, after primary results.
# Holm across 432 protects pair-specific claims/counts across the six networks,
# separately at each bound. Also export raw and within-network values transparently.
tost_rows <- list()
for(bound in c(.05,.03)) {
  e <- est;delta <- atanh(bound);e$bound_r <- bound
  e$p_lower <- pnorm((e$fisher_z+delta)/e$SE_z,lower.tail=FALSE)
  e$p_upper <- pnorm((e$fisher_z-delta)/e$SE_z)
  e$p_TOST <- pmax(e$p_lower,e$p_upper)
  e$p_Holm432 <- p.adjust(e$p_TOST,'holm')
  e$p_Holm72 <- NA_real_
  for(key in unique(e$network)) {
    ix <- e$network==key;e$p_Holm72[ix] <- p.adjust(e$p_TOST[ix],'holm')
  }
  e$equivalent_raw <- e$p_TOST<.05
  e$equivalent_Holm432 <- e$p_Holm432<.05
  e$equivalent_Holm72 <- e$p_Holm72<.05
  if(any(e$equivalent_raw!=(e$lower90 > -bound & e$upper90 < bound)))
    stop('TOST/90-percent interval consistency check failed.')
  tost_rows[[as.character(bound)]] <- e
}
tost <- dplyr::bind_rows(tost_rows)
csv(tost,'13_equivalence_TOST_both_VAS_bounds.csv')
summary_rows <- list()
for(bound in c(.05,.03)) for(key in c(unique(est$network),'ALL_SIX')) {
  e <- tost[tost$bound_r==bound & (key=='ALL_SIX' | tost$network==key),]
  summary_rows[[length(summary_rows)+1]] <- data.frame(network=key,bound_r=bound,
    pairs=nrow(e),max_abs_partial_r=max(abs(e$partial_r)),
    equivalent_raw=sum(e$equivalent_raw),equivalent_Holm432=sum(e$equivalent_Holm432),
    equivalent_Holm72=sum(e$equivalent_Holm72),
    positive=sum(e$partial_r>0),negative=sum(e$partial_r<0),
    interpretation='within stated bounds only; conditional approximate inference')
}
csv(dplyr::bind_rows(summary_rows),'13_equivalence_summary.csv')
pdf(file.path(OUT,'13_unregularised_cross_domain_intervals.pdf'),width=10,height=14)
for(key in names(nets)) {
  e <- est[est$network==key,];e <- e[order(e$partial_r),];y <- seq_len(72)
  par(mar=c(5,15,3,1))
  plot(e$partial_r,y,xlim=range(c(e$lower95,e$upper95,-.055,.055)),
    yaxt='n',ylab='',xlab='Unregularised partial correlation; pointwise 95% interval',
    main=key,pch=16,cex=.6)
  segments(e$lower95,y,e$upper95,y,col='grey40')
  axis(2,at=y,labels=paste(e$node1,e$node2,sep=' - '),las=2,cex.axis=.48)
  abline(v=0,col='grey60');abline(v=c(-.05,.05),col='blue',lty=2)
  abline(v=c(-.03,.03),col='darkorange',lty=3)
}
dev.off()

# B. Reconstruct wave-specific and no-family pooled sensitivities on primary IDs.
# Within one visit, one participant/family makes a family random intercept unidentifiable.
# The existing residualise() function therefore omits that within-wave term.
sets <- list(primary=nets,wave_specific=list(),no_family_pooled=list())
diagnostics <- list()
for(variant in c('wave_specific','no_family_pooled')) for(mod in names(MRI)) {
  z <- res[[mod]]$source
  if(!all(flags_for(z,mod)$base)) stop('Invalid eligible/QC adjustment sample.')
  rr <- residualise(z,mod,scope=if(variant=='wave_specific')'wave' else 'pooled',
    family_model=variant=='wave_specific',extraction='retain_family')
  dg <- rr$diagnostics;dg$variant <- variant;diagnostics[[paste(variant,mod)]] <- dg
  saveRDS(list(source=z,residuals=rr$data,models=rr$models),
    file.path(private,paste0('13_',variant,'_',mod,'_residuals_RESTRICTED.rds')))
  for(w in WAVES) {
    key <- paste(mod,w,sep='_');p <- nets[[key]]
    take <- z$visit==w & complete.cases(rr$data);ids <- z$participant_id[take]
    if(anyDuplicated(ids) || !setequal(ids,p$ids)) stop('Sensitivity sample changed: ',key)
    x <- rr$data[take,,drop=FALSE][match(p$ids,ids),,drop=FALSE]
    sets[[variant]][[key]] <- list(x=x,W=network(x),ids=p$ids,n=nrow(x))
  }
}
csv(dplyr::bind_rows(diagnostics),'13_sensitivity_residualisation_diagnostics.csv')
saveRDS(sets,file.path(private,'13_network_sets_RESTRICTED.rds'))

# Centrality comparisons need only saved edge matrices; no extra model is fitted.
cent <- comparison <- list()
add_central <- function(W,key,variant,reference=NULL) {
  c <- central(W);c$network <- key;c$variant <- variant
  metrics <- c('strength','expected_influence','bridge_strength','bridge_expected_influence')
  for(m in metrics)c[[paste0(m,'_rank')]] <- rank(-c[[m]],ties.method='average')
  cent[[length(cent)+1]] <<- c
  if(!is.null(reference)) for(m in metrics) {
    r <- central(reference)[[m]];v <- c[[m]]
    comparison[[length(comparison)+1]] <<- data.frame(network=key,variant=variant,
      metric=m,reference='primary',spearman=if(sd(r)>0 && sd(v)>0)cor(r,v,method='spearman')else NA_real_,
      max_abs_difference=max(abs(r-v)))
  }
}
for(variant in names(sets)) for(key in names(nets))
  add_central(sets[[variant]][[key]]$W,key,variant,
    if(variant=='primary')NULL else nets[[key]]$W)
sens <- read.csv(file.path(V2,'06_sensitivity_edges.csv'))
for(key in names(nets)) {
  e <- sens[sens$variant=='ordinal_latent_age_adjusted' & sens$network==key,]
  if(nrow(e)!=153L || anyDuplicated(e[c('i','j')])) stop('Ordinal edges missing.')
  W <- matrix(0,18,18,dimnames=list(NODES,NODES));W[cbind(e$i,e$j)] <- e$weight
  W <- W+t(W);add_central(W,key,'ordinal_latent_age_adjusted',nets[[key]]$W)
}
csv(dplyr::bind_rows(cent),'13_primary_and_sensitivity_centrality.csv')
csv(dplyr::bind_rows(comparison),'13_centrality_comparisons.csv')

# C. Original paired NCT: 1000 permutations; BH/Holm within 153 edges/comparison.
extract_p <- function(nt) {
  p <- nt$einv.pvals
  if(is.matrix(p) && identical(dim(p),c(18L,18L))) {
    ij <- which(upper.tri(p),arr.ind=TRUE)
    a <- data.frame(i=ij[,1],j=ij[,2],p_raw=p[ij])
  } else {
    p <- as.data.frame(p,check.names=FALSE)
    if(!all(c('Var1','Var2','p-value') %in% names(p))) stop('Unexpected NCT edge table.')
    index <- function(v) {
      s <- as.character(v);a <- match(s,NODES);miss <- is.na(a)
      a[miss] <- suppressWarnings(as.integer(s[miss]))
      if(anyNA(a) || any(a<1 | a>18)) stop('Unknown NCT nodes.');a
    }
    i <- index(p$Var1);j <- index(p$Var2)
    a <- data.frame(i=pmin(i,j),j=pmax(i,j),p_raw=as.numeric(p[['p-value']]))
  }
  if(nrow(a)!=153L || anyDuplicated(a[c('i','j')]) ||
    any(!is.finite(a$p_raw)) || any(a$p_raw<0 | a$p_raw>1)) stop('Invalid edge p-values.')
  a
}
omnibus <- edge_tests <- list()
for(variant in c('wave_specific','no_family_pooled')) for(mod in names(MRI))
for(pair in combn(WAVES,2,simplify=FALSE)) {
  a <- sets[[variant]][[paste(mod,pair[1],sep='_')]]
  b <- sets[[variant]][[paste(mod,pair[2],sep='_')]]
  ids <- sort(intersect(a$ids,b$ids))
  if(length(ids)<CFG$min_network_n) stop('Insufficient paired sample.')
  A <- a$x[match(ids,a$ids),,drop=FALSE];B <- b$x[match(ids,b$ids),,drop=FALSE]
  label <- paste(variant,mod,paste(pair,collapse='_vs_'),sep='_')
  checkpoint <- file.path(private,paste0('13_',label,'_NCT.rds'))
  if(file.exists(checkpoint)) nt <- readRDS(checkpoint) else {
    message('13 NCT: ',label,'; N=',length(ids))
    set.seed(CFG$seed+1300L+length(omnibus))
    nt <- NetworkComparisonTest::NCT(A,B,gamma=.5,it=1000L,binary.data=FALSE,
      paired=TRUE,weighted=TRUE,test.edges=TRUE,edges='all',
      progressbar=TRUE,p.adjust.methods='none')
    saveRDS(nt,checkpoint)
  }
  omnibus[[label]] <- data.frame(variant=variant,modality=mod,visit_A=pair[1],
    visit_B=pair[2],paired_N=length(ids),permutations=1000L,
    M=as.numeric(nt$nwinv.real),p_structure=as.numeric(nt$nwinv.pval),
    global_A=as.numeric(nt$glstrinv.sep[1]),global_B=as.numeric(nt$glstrinv.sep[2]),
    p_global=as.numeric(nt$glstrinv.pval),scope='joint network; omnibus p-values unadjusted')
  e <- merge(edges(network(A)),extract_p(nt),by=c('i','j'),sort=FALSE)
  eb <- edges(network(B));e$weight_A <- e$weight
  e$weight_B <- eb$weight[match(paste(e$i,e$j),paste(eb$i,eb$j))];e$weight <- NULL
  e$variant <- variant;e$comparison <- label
  e$p_BH_within153 <- p.adjust(e$p_raw,'BH')
  e$p_Holm_within153 <- p.adjust(e$p_raw,'holm');edge_tests[[label]] <- e
}
csv(dplyr::bind_rows(omnibus),'13_sensitivity_NCT_omnibus.csv')
csv(dplyr::bind_rows(edge_tests),'13_sensitivity_NCT_edges.csv')

# D. Original wave-specific nonparametric and case-dropping bootstraps (1000).
# Primary bootstraps have already been completed; do not repeat them.
cs_rows <- interval_rows <- list()
for(key in names(nets)) {
  edge_checkpoint <- file.path(private,paste0('13_',key,'_wave_edge_bootnet.rds'))
  if(file.exists(edge_checkpoint)) be <- readRDS(edge_checkpoint) else {
    message('13 wave-specific edge bootstrap: ',key)
    set.seed(CFG$seed+1450L+length(interval_rows))
    en <- bootnet::estimateNetwork(as.data.frame(sets$wave_specific[[key]]$x),
      default='EBICglasso',tuning=.5,corMethod='cor')
    be <- bootnet::bootnet(en,nBoots=1000L,type='nonparametric',nCores=1L,
      statistics=c('edge','strength','expectedInfluence'))
    saveRDS(be,edge_checkpoint)
  }
  bt <- be$bootTable
  if(!all(c('type','node1','node2','value') %in% names(bt))) stop('Unexpected bootstrap table.')
  bt <- bt[bt$type=='edge',c('node1','node2','value'),drop=FALSE]
  ei <- bt |> dplyr::group_by(node1,node2) |>
    dplyr::summarise(lower=quantile(value,.025,na.rm=TRUE),
      upper=quantile(value,.975,na.rm=TRUE),successful_values=sum(is.finite(value)),.groups='drop')
  if(nrow(ei)!=153L || any(ei$successful_values!=1000L)) stop('Incomplete wave bootstrap.')
  ei$network <- key;ei$scope <- 'regularised; pointwise; conditional on preprocessing; not significance tests'
  interval_rows[[key]] <- ei
  pdf(file.path(OUT,paste0('13_',key,'_wave_edge_accuracy.pdf')))
  print(plot(be,labels=FALSE,order='sample'));dev.off()
  checkpoint <- file.path(private,paste0('13_',key,'_wave_case_bootnet.rds'))
  if(file.exists(checkpoint)) bc <- readRDS(checkpoint) else {
    message('13 wave-specific case bootstrap: ',key)
    set.seed(CFG$seed+1400L+length(cs_rows))
    en <- bootnet::estimateNetwork(as.data.frame(sets$wave_specific[[key]]$x),
      default='EBICglasso',tuning=.5,corMethod='cor')
    bc <- bootnet::bootnet(en,nBoots=1000L,type='case',nCores=1L,
      statistics=c('strength','expectedInfluence','edge','closeness','betweenness'))
    saveRDS(bc,checkpoint)
  }
  cs <- bootnet::corStability(bc,statistics=c('strength','expectedInfluence'))
  cs_rows[[key]] <- data.frame(network=key,statistic=names(cs),CS=as.numeric(cs))
  pdf(file.path(OUT,paste0('13_',key,'_wave_case_stability.pdf')))
  print(plot(bc,statistics=c('strength','expectedInfluence')));dev.off()
}
csv(dplyr::bind_rows(cs_rows),'13_wave_specific_centrality_CS.csv')
csv(dplyr::bind_rows(interval_rows),'13_wave_specific_edge_intervals.csv')
writeLines(c(
 'Supplementary revision analyses; completed primary v2 estimates remain unchanged.',
 'TOST targets unregularised partial correlations conditional on other 16 nodes.',
 'Fisher-z normal approximation; conditional on fitted preprocessing, whose uncertainty is not propagated.',
 'Discrete SMA scores and distributional assumptions limit this approximation.',
 'VAS proposed bounds .05 and .03; both reported, selected during revision; not validated clinical cutoffs.',
 'Use Holm432 columns for individual equivalence claims/counts across six networks at each bound.',
 'Raw 90% intervals correspond to raw TOST only; 95% intervals are pointwise, not simultaneous.',
 'Do not claim zero association, clinical irrelevance or equivalence of lasso weights.',
 'Sensitivity NCTs retain 1000 permutations and BH/Holm within 153 edges.',
 '1000 permutations limit attainable p-values and Holm resolution.',
 'Ordinal sensitivity changes SMA adjustment/estimand; it is not a replacement primary model.',
 'Motion source-scan provenance remains external verification.'
),file.path(OUT,'13_READ_ME.txt'))
save_environment('13')
files <- list.files(OUT,pattern='^13_.*[.](csv|pdf|txt)$',full.names=TRUE)
if(!all(file.copy(files,share,overwrite=TRUE))) stop('Aggregate export failed.')
message('13 complete. Send ONLY: ',share)
