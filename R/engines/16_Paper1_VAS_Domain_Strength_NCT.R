# VAS domain decomposition: fit all 18 nodes, then sum only the requested edges.
# Uses the installed NCT package's own paired permutation procedure.
message('Choose the existing 00_Paper1_Settings_and_Functions.R.')
p00 <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R')stop('Wrong settings file.')
CODE_DIR <- dirname(p00);source(p00,local=TRUE);need('NetworkComparisonTest')
V2 <- file.path(CODE_DIR,'Paper1_revision_outputs_original_methods_v2')
CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(CFG$pooled_wave_fixed,FALSE) || CFG$gamma!=.5)stop('Unexpected configuration.')
OUT <- file.path(CODE_DIR,'Paper1_VAS_remaining_checks')
private <- file.path(OUT,'RESTRICTED');share <- file.path(OUT,'SHAREABLE_VAS_CHECKS')
dir.create(private,recursive=TRUE,showWarnings=FALSE);dir.create(share,showWarnings=FALSE)
nets <- readRDS(file.path(V2,'03_primary_networks_RESTRICTED.rds'))
signature <- list(version='16-v1',source=tools::md5sum(file.path(V2,
 c('03_primary_networks_RESTRICTED.rds','03_settings.rds'))),NCT=as.character(packageVersion('NetworkComparisonTest')))
sigfile <- file.path(private,'16_input_signature.rds')
if(file.exists(sigfile)) {
 if(!identical(readRDS(sigfile),signature))stop('Inputs/package changed; use a new output folder.')
}else saveRDS(signature,sigfile)
native <- getFromNamespace('NCT_estimator_GGM','NetworkComparisonTest')
domain_mask <- function(domain) {
 e <- edges(matrix(0,18,18));e <- e[e$domain==domain,]
 m <- matrix(FALSE,18,18);m[cbind(e$i,e$j)] <- TRUE;m | t(m)
}
domain_estimator <- function(x,domain) {
 # Same full GGM estimator/settings used by the original NCT.
 W <- native(x,gamma=.5,make.positive.definite=TRUE,corMethod='cor',verbose=FALSE)
 W[!domain_mask(domain)] <- 0;diag(W) <- 0;W
}
domains <- c('SMA-SMA','MRI-MRI','SMA-MRI');rows <- list()
for(mod in names(MRI))for(pair in combn(WAVES,2,simplify=FALSE)) {
 a <- nets[[paste(mod,pair[1],sep='_')]];b <- nets[[paste(mod,pair[2],sep='_')]]
 ids <- sort(intersect(a$ids,b$ids))
 if(length(ids)<CFG$min_network_n)stop('Too few paired participants.')
 A <- a$x[match(ids,a$ids),,drop=FALSE];B <- b$x[match(ids,b$ids),,drop=FALSE]
 full_A <- native(A,gamma=.5,make.positive.definite=TRUE,corMethod='cor',verbose=FALSE)
 full_B <- native(B,gamma=.5,make.positive.definite=TRUE,corMethod='cor',verbose=FALSE)
 strengths_A <- strengths(full_A);strengths_B <- strengths(full_B)
 if(abs(sum(strengths_A[domains])-strengths_A['global'])>1e-10 ||
    abs(sum(strengths_B[domains])-strengths_B['global'])>1e-10)
   stop('Domain/global decomposition failed.')
 # Same seed per paired comparison gives the same permutation sequence to all domains.
 pair_seed <- CFG$seed+1600L+match(mod,names(MRI))*10L+
   match(paste(pair,collapse='_'),vapply(combn(WAVES,2,simplify=FALSE),paste,character(1),collapse='_'))
 for(domain in domains) {
  label <- paste(mod,paste(pair,collapse='_vs_'),domain,sep='_')
  checkpoint <- file.path(private,paste0('16_',label,'_NCT.rds'))
  if(file.exists(checkpoint))nt <- readRDS(checkpoint)else {
   message('16 domain NCT: ',label);set.seed(pair_seed)
   nt <- NetworkComparisonTest::NCT(A,B,it=1000L,paired=TRUE,weighted=TRUE,
    abs=TRUE,test.edges=FALSE,test.centrality=FALSE,progressbar=TRUE,
    estimator=domain_estimator,estimatorArgs=list(domain=domain),p.adjust.methods='none')
   saveRDS(nt,checkpoint)
  }
  if(max(abs(as.numeric(nt$glstrinv.sep)-c(strengths_A[domain],strengths_B[domain])))>1e-8)
    stop('Domain NCT strength does not match full-network decomposition.')
  perm <- as.numeric(nt$glstrinv.perm)
  if(length(perm)!=1000L || any(!is.finite(perm)))stop('Invalid domain permutation draws.')
  rows[[label]] <- data.frame(modality=mod,visit_A=pair[1],visit_B=pair[2],domain=domain,
   paired_N=length(ids),permutations=1000L,
   possible_edges=switch(domain,'SMA-SMA'=15L,'MRI-MRI'=66L,'SMA-MRI'=72L),
   strength_A=as.numeric(strengths_A[domain]),strength_B=as.numeric(strengths_B[domain]),
   signed_difference=as.numeric(strengths_B[domain]-strengths_A[domain]),
   absolute_difference=as.numeric(nt$glstrinv.real),p_raw=as.numeric(nt$glstrinv.pval),
   estimand='sum of absolute domain weights from full 18-node GGM; conditional on preprocessing')
 }
}
result <- dplyr::bind_rows(rows)
if(nrow(result)!=18L || any(!is.finite(result$p_raw)))stop('Expected 18 valid domain tests.')
# New domain family: all 18 aggregate tests, distinct from original edge-test families.
result$p_BH18 <- p.adjust(result$p_raw,'BH')
result$p_Holm18 <- p.adjust(result$p_raw,'holm')
csv(result,'16_domain_strength_NCT_18_tests.csv')
writeLines(c('Full 18-node networks estimated before selecting domain weights; zero edges included.',
 'Native NCT estimator and paired permutation procedure; 1000 permutations.',
 'All 18 aggregate tests form a separate family; report raw, BH18 and Holm18 transparently.',
 'Does not change original within-153 edge corrections or unadjusted omnibus tests.',
 'A nonsignificant domain test does not establish equivalence or absence of difference.',
 'Dependent-data NCT and preprocessing have the same limitations as the primary paired analysis.'
),file.path(OUT,'16_READ_ME.txt'))
save_environment('16')
files <- list.files(OUT,pattern='^16_.*[.](csv|txt)$',full.names=TRUE)
if(!all(file.copy(files,share,overwrite=TRUE)))stop('Aggregate export failed.')
message('16 complete. Send ONLY: ',share)
