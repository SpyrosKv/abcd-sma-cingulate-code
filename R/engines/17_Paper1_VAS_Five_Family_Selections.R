# Five additional random selections, fixed across visits within each repetition.
# Re-estimate the PRIMARY pipeline, centrality, bootstraps and original paired NCTs.
# Repetitions are sensitivity analyses, not five independent study samples.
message('Choose the existing 00_Paper1_Settings_and_Functions.R.')
p00 <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R')stop('Wrong settings file.')
CODE_DIR <- dirname(p00);source(p00,local=TRUE);need(c('NetworkComparisonTest','bootnet'))
V2 <- file.path(CODE_DIR,'Paper1_revision_outputs_original_methods_v2')
CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(CFG$pooled_wave_fixed,FALSE) || CFG$gamma!=.5)stop('Unexpected primary configuration.')
OUT <- file.path(CODE_DIR,'Paper1_VAS_remaining_checks')
private <- file.path(OUT,'RESTRICTED');share <- file.path(OUT,'SHAREABLE_VAS_CHECKS')
dir.create(private,recursive=TRUE,showWarnings=FALSE);dir.create(share,showWarnings=FALSE)
d <- readRDS(file.path(V2,'01_merged_RESTRICTED.rds'))
primary <- readRDS(file.path(V2,'03_primary_networks_RESTRICTED.rds'))
family_seeds <- CFG$seed+1700L+seq_len(5)
signature <- list(version='17-v1',seeds=family_seeds,inputs=tools::md5sum(file.path(V2,
 c('01_merged_RESTRICTED.rds','03_primary_networks_RESTRICTED.rds','03_settings.rds'))),
 helper=tools::md5sum(p00),NCT=as.character(packageVersion('NetworkComparisonTest')),
 bootnet=as.character(packageVersion('bootnet')))
sigfile <- file.path(private,'17_input_signature.rds')
if(file.exists(sigfile)) {
 if(!identical(readRDS(sigfile),signature))stop('Inputs/settings changed; use a new output folder.')
}else saveRDS(signature,sigfile)
fam <- unique(d[c('participant_id','family')]);fam <- fam[!is.na(fam$family),]
if(anyDuplicated(fam$participant_id))stop('Participant has multiple family IDs.')
fam <- fam[order(as.character(fam$family),fam$participant_id),]
groups <- split(fam$participant_id,fam$family,drop=TRUE)
original_selected <- unique(d$participant_id[d$selected_family])

extract_p <- function(nt) {
 p <- nt$einv.pvals
 if(is.matrix(p) && identical(dim(p),c(18L,18L))) {
  ij <- which(upper.tri(p),arr.ind=TRUE);e <- data.frame(i=ij[,1],j=ij[,2],p_raw=p[ij])
 }else {
  p <- as.data.frame(p,check.names=FALSE)
  if(!all(c('Var1','Var2','p-value') %in% names(p)))stop('Unexpected NCT edge table.')
  index <- function(v) {
   s <- as.character(v);i <- match(s,NODES);missing <- is.na(i)
   i[missing] <- suppressWarnings(as.integer(s[missing]))
   if(anyNA(i) || any(i<1 | i>18))stop('Unknown NCT nodes.');i
  }
  i <- index(p$Var1);j <- index(p$Var2)
  e <- data.frame(i=pmin(i,j),j=pmax(i,j),p_raw=as.numeric(p[['p-value']]))
 }
 if(nrow(e)!=153L || anyDuplicated(e[c('i','j')]) || any(!is.finite(e$p_raw)) ||
   any(e$p_raw<0 | e$p_raw>1))stop('Invalid edge p-values.')
 e
}
family_audit <- flow <- diag_rows <- comparisons <- edge_rows <- central_rows <- list()
omnibus <- nct_edges <- cs_rows <- intervals <- list()
for(rep in seq_len(5)) {
 seed <- family_seeds[rep];label <- paste0('selection_',rep)
 run_dir <- file.path(private,paste0('17_',label));dir.create(run_dir,showWarnings=FALSE)
 set.seed(seed)
 selected <- unlist(lapply(groups,function(ids)sample(sort(ids),1)),use.names=FALSE)
 selected_families <- fam$family[match(selected,fam$participant_id)]
 if(anyNA(selected_families) || anyDuplicated(selected) || anyDuplicated(selected_families))
   stop('Family selection integrity failed.')
 dd <- d;dd$selected_family <- dd$participant_id %in% selected
 saveRDS(data.frame(participant_id=selected,family=selected_families),
   file.path(run_dir,'selected_participants_RESTRICTED.rds'))
 family_audit[[label]] <- data.frame(repetition=rep,seed=seed,N_selected=length(selected),
   N_families=length(unique(selected_families)),N_multichild_families=sum(lengths(groups)>1),
   different_selected_from_primary=sum(!selected %in% original_selected))
 checkpoint <- file.path(run_dir,'networks_and_diagnostics_RESTRICTED.rds')
 if(file.exists(checkpoint))obj <- readRDS(checkpoint)else {
  nets <- diagnostics <- list()
  for(mod in names(MRI)) {
   f <- flags_for(dd,mod);z <- droplevels(dd[f$base,,drop=FALSE])
   rr <- residualise(z,mod,scope='pooled',family_model=TRUE,extraction='retain_family')
   diagnostics[[mod]] <- rr$diagnostics
   for(w in WAVES) {
    use <- z$visit==w & complete.cases(rr$data);x <- rr$data[use,,drop=FALSE]
    ids <- z$participant_id[use]
    if(anyDuplicated(ids) || anyDuplicated(z$family[use]))stop('Non-independent network sample.')
    key <- paste(mod,w,sep='_')
    nets[[key]] <- list(x=x,W=network(x),ids=ids,n=nrow(x))
   }
   saveRDS(list(source=z,residuals=rr$data,models=rr$models),
     file.path(run_dir,paste0(mod,'_residuals_RESTRICTED.rds')))
  }
  obj <- list(networks=nets,diagnostics=dplyr::bind_rows(diagnostics));saveRDS(obj,checkpoint)
 }
 nets <- obj$networks;dg <- obj$diagnostics;dg$repetition <- rep
 if(any(grepl('visit_factor',dg$formula)) || any(dg$extraction!='retain_family'))
   stop('Unexpected residualisation.')
 diag_rows[[label]] <- dg
 for(mod in names(MRI))for(w in WAVES) {
  zz <- dd[dd$visit==w,,drop=FALSE];f <- flags_for(zz,mod)
  stage <- rep(TRUE,nrow(zz));counts <- c(source=nrow(zz))
  for(k in c('eligible','family_ok','qc_ok','cov_ok','mot_ok','mri_complete','sma_complete')) {
   stage <- stage & f[[k]];counts[k] <- sum(stage)
  }
  flow[[paste(label,mod,w)]] <- data.frame(repetition=rep,modality=mod,visit=w,
   step=names(counts),N=as.integer(counts))
 }
 for(key in names(nets)) {
  n <- nets[[key]];p <- primary[[key]]
  e <- edges(n$W);e$network <- key;e$repetition <- rep;e$N <- n$n
  edge_rows[[paste(label,key)]] <- e
  ct <- central(n$W);ct$network <- key;ct$repetition <- rep
  central_rows[[paste(label,key)]] <- ct
  co <- compare_networks(p$W,n$W,paste(label,key))
  co$network <- key;co$repetition <- rep;co$N_primary <- p$n;co$N_repetition <- n$n
  co$common_participants <- length(intersect(p$ids,n$ids))
  co$scope <- 'different family selections; descriptive comparison of complete primary samples'
  comparisons[[paste(label,key)]] <- co

  # Nonparametric edge bootstrap, original settings.
  file_be <- file.path(run_dir,paste0(key,'_edge_bootnet.rds'))
  if(file.exists(file_be))be <- readRDS(file_be)else {
   message('17 selection ',rep,' edge bootstrap: ',key)
   set.seed(seed+100L+match(key,names(nets)))
   en <- bootnet::estimateNetwork(as.data.frame(n$x),default='EBICglasso',tuning=.5,corMethod='cor')
   be <- bootnet::bootnet(en,nBoots=1000L,type='nonparametric',nCores=1L,
     statistics=c('edge','strength','expectedInfluence'));saveRDS(be,file_be)
  }
  bt <- be$bootTable
  if(!all(c('type','node1','node2','value') %in% names(bt)))stop('Unexpected bootstrap table.')
  bt <- bt[bt$type=='edge',c('node1','node2','value'),drop=FALSE]
  ei <- bt |> dplyr::group_by(node1,node2) |>
    dplyr::summarise(lower=quantile(value,.025,na.rm=TRUE),upper=quantile(value,.975,na.rm=TRUE),
      successful_values=sum(is.finite(value)),.groups='drop')
  if(nrow(ei)!=153L || any(ei$successful_values!=1000L))stop('Incomplete edge bootstrap.')
  ei$network <- key;ei$repetition <- rep
  ei$scope <- 'regularised; pointwise; conditional on preprocessing; not significance tests'
  intervals[[paste(label,key)]] <- ei
  file_bc <- file.path(run_dir,paste0(key,'_case_bootnet.rds'))
  if(file.exists(file_bc))bc <- readRDS(file_bc)else {
   message('17 selection ',rep,' case bootstrap: ',key)
   set.seed(seed+200L+match(key,names(nets)))
   en <- bootnet::estimateNetwork(as.data.frame(n$x),default='EBICglasso',tuning=.5,corMethod='cor')
   bc <- bootnet::bootnet(en,nBoots=1000L,type='case',nCores=1L,
     statistics=c('strength','expectedInfluence','edge','closeness','betweenness'))
   saveRDS(bc,file_bc)
  }
  cs <- bootnet::corStability(bc,statistics=c('strength','expectedInfluence'))
  cs_rows[[paste(label,key)]] <- data.frame(repetition=rep,network=key,
   statistic=names(cs),CS=as.numeric(cs))
  rm(be,bc);gc()
 }
 # Original paired NCTs within each repetition; no tests between overlapping repetitions.
 for(mod in names(MRI))for(pair in combn(WAVES,2,simplify=FALSE)) {
  a <- nets[[paste(mod,pair[1],sep='_')]];b <- nets[[paste(mod,pair[2],sep='_')]]
  ids <- sort(intersect(a$ids,b$ids))
  if(length(ids)<CFG$min_network_n)stop('Too few paired participants.')
  A <- a$x[match(ids,a$ids),,drop=FALSE];B <- b$x[match(ids,b$ids),,drop=FALSE]
  nm <- paste(mod,paste(pair,collapse='_vs_'),sep='_')
  checkpoint <- file.path(run_dir,paste0(nm,'_NCT.rds'))
  if(file.exists(checkpoint))nt <- readRDS(checkpoint)else {
   message('17 selection ',rep,' NCT: ',nm)
   set.seed(seed+300L+length(omnibus))
   nt <- NetworkComparisonTest::NCT(A,B,gamma=.5,it=1000L,binary.data=FALSE,
    paired=TRUE,weighted=TRUE,test.edges=TRUE,edges='all',progressbar=TRUE,p.adjust.methods='none')
   saveRDS(nt,checkpoint)
  }
  omnibus[[paste(label,nm)]] <- data.frame(repetition=rep,modality=mod,
   visit_A=pair[1],visit_B=pair[2],paired_N=length(ids),permutations=1000L,
   M=as.numeric(nt$nwinv.real),p_structure=as.numeric(nt$nwinv.pval),
   global_A=as.numeric(nt$glstrinv.sep[1]),global_B=as.numeric(nt$glstrinv.sep[2]),
   p_global=as.numeric(nt$glstrinv.pval))
  e <- merge(edges(network(A)),extract_p(nt),by=c('i','j'),sort=FALSE)
  eb <- edges(network(B));e$weight_A <- e$weight
  e$weight_B <- eb$weight[match(paste(e$i,e$j),paste(eb$i,eb$j))];e$weight <- NULL
  e$repetition <- rep;e$comparison <- nm
  e$p_BH_within153 <- p.adjust(e$p_raw,'BH');e$p_Holm_within153 <- p.adjust(e$p_raw,'holm')
  nct_edges[[paste(label,nm)]] <- e
 }
 # Export completed repetitions progressively; restart resumes checkpoints.
 csv(dplyr::bind_rows(family_audit),'17_family_selection_audit.csv')
 csv(dplyr::bind_rows(flow),'17_family_selection_sample_flow.csv')
 csv(dplyr::bind_rows(diag_rows),'17_family_selection_residualisation.csv')
 csv(dplyr::bind_rows(comparisons),'17_family_selection_network_comparisons.csv')
 csv(dplyr::bind_rows(edge_rows),'17_family_selection_edges.csv')
 csv(dplyr::bind_rows(central_rows),'17_family_selection_centrality.csv')
 csv(dplyr::bind_rows(intervals),'17_family_selection_edge_intervals.csv')
 csv(dplyr::bind_rows(cs_rows),'17_family_selection_centrality_CS.csv')
 csv(dplyr::bind_rows(omnibus),'17_family_selection_NCT_omnibus.csv')
 csv(dplyr::bind_rows(nct_edges),'17_family_selection_NCT_edges.csv')
}
if(length(family_audit)!=5L || length(omnibus)!=30L || length(cs_rows)!=30L)
 stop('Five repetitions did not complete.')
writeLines(c('Five additional family selections; selection fixed across visits within each repetition.',
 'Selection precedes eligibility/QC and is not chosen to maximise usable data.',
 'Primary pooled residualisation retains family contribution; no wave fixed effect.',
 'Each repetition: six networks, original nonparametric/case bootstraps and six paired NCTs.',
 '1000 bootstrap samples/permutations; BH/Holm within 153 edges/comparison.',
 'Repetitions overlap and are not independent samples; no pooling or between-repetition significance tests.',
 'Sensitivity results cannot be selected according to favourable p-values.',
 'Report all five selections and uncertainty; MI/ordinal/domain analyses are not repeated here.'
),file.path(OUT,'17_READ_ME.txt'))
save_environment('17')
files <- list.files(OUT,pattern='^17_.*[.](csv|txt)$',full.names=TRUE)
if(!all(file.copy(files,share,overwrite=TRUE)))stop('Aggregate export failed.')
message('17 complete. Send ONLY: ',share)
