# Correct replacement for the inherited mixed pooled-PMM / wave-complete-case comparison.
# Run after 13. Original residual-scale PMM: m=20, seed=500, 50/20 iterations.
message('Choose the existing 00_Paper1_Settings_and_Functions.R.')
p00 <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p00)!='00_Paper1_Settings_and_Functions.R') stop('Wrong settings file.')
CODE_DIR <- dirname(p00);source(p00,local=TRUE)
V2 <- file.path(CODE_DIR,'Paper1_revision_outputs_original_methods_v2')
CFG <- readRDS(file.path(V2,'03_settings.rds'))
if(!identical(CFG$pooled_wave_fixed,FALSE) || CFG$gamma!=.5) stop('Unexpected configuration.')
OUT <- file.path(CODE_DIR,'Paper1_VAS_remaining_checks')
private <- file.path(OUT,'RESTRICTED');share <- file.path(OUT,'SHAREABLE_VAS_CHECKS')
required <- file.path(private,c('13_network_sets_RESTRICTED.rds',
  '13_wave_specific_FA_residuals_RESTRICTED.rds',
  '13_wave_specific_T2w_residuals_RESTRICTED.rds'))
if(!all(file.exists(required))) stop('Complete 13 first.')
sets <- readRDS(required[1]);audit <- comparisons <- all_edges <- all_central <- logged <- list()
signature <- list(version='14-v1',source_md5=tools::md5sum(required))
sigfile <- file.path(private,'14_input_signature.rds')
if(file.exists(sigfile)) {
 if(!identical(readRDS(sigfile),signature)) stop('Inputs changed; do not mix PMM runs.')
} else saveRDS(signature,sigfile)
for(mod in names(MRI)) {
 rs <- readRDS(file.path(private,paste0('13_wave_specific_',mod,'_residuals_RESTRICTED.rds')))
 z <- rs$source;r <- rs$residuals
 if(nrow(z)!=nrow(r) || !all(flags_for(z,mod)$base)) stop('Invalid adjustment population.')
 for(w in WAVES) for(rule in c('AND','OR')) {
  has_sma <- rowSums(!is.na(r[,1:6,drop=FALSE]))>0
  has_mri <- rowSums(!is.na(r[,7:18,drop=FALSE]))>0
  available <- if(rule=='AND')has_sma & has_mri else has_sma | has_mri
  take <- z$visit==w & available;x <- as.data.frame(r[take,,drop=FALSE])
  key <- paste(mod,w,sep='_');label <- paste(key,rule,sep='_')
  if(nrow(x)<CFG$min_network_n || any(colSums(!is.na(x))<10L)) stop('Too few donors: ',label)
  # Complete-case reference uses the SAME wave-specific adjusted values.
  cc <- complete.cases(x);ids_cc <- z$participant_id[take][cc]
  reference <- sets$wave_specific[[key]]
  if(!setequal(ids_cc,reference$ids)) stop('Complete-case sample changed: ',label)
  cc_data <- x[cc,,drop=FALSE][match(reference$ids,ids_cc),,drop=FALSE]
  if(max(abs(network(cc_data)-reference$W))>1e-10) stop('Reference reconstruction failed.')
  iterations <- if(mod=='FA' && w=='ses-02A')50L else 20L
  checkpoint <- file.path(private,paste0('14_',label,'_PMM_RESTRICTED.rds'))
  if(file.exists(checkpoint)) obj <- readRDS(checkpoint) else {
   message('14 wave-specific PMM: ',label)
   imp <- if(anyNA(x))mice::mice(x,m=20L,method='pmm',maxit=iterations,
     seed=500L,printFlag=TRUE)else NULL
   obj <- list(imputation=imp,complete_data=if(is.null(imp))x else NULL)
   saveRDS(obj,checkpoint)
  }
  imp <- obj$imputation
  saved_x <- if(is.null(imp))obj$complete_data else imp$data
  if(!isTRUE(all.equal(as.matrix(saved_x),as.matrix(x),check.attributes=FALSE)))
   stop('Saved PMM input mismatch: ',label)
  if(!is.null(imp) && (imp$m!=20L || imp$iteration!=iterations)) stop('Wrong PMM settings.')
  events <- if(is.null(imp))NULL else imp$loggedEvents
  if(!is.null(events)) {events$network_rule <- label;logged[[label]] <- events}
  if(!is.null(imp)) {
   pdf(file.path(OUT,paste0('14_',label,'_trace.pdf')))
   print(plot(imp));dev.off()
  }
  for(k in 1:20) {
   completed <- if(is.null(imp))x else mice::complete(imp,k)
   if(anyNA(completed) || any(!is.finite(as.matrix(completed)))) stop('Invalid PMM completion.')
   W <- network(completed);e <- edges(W);e$network <- key;e$rule <- rule
   e$imputation <- k;e$N <- nrow(completed);all_edges[[paste(label,k)]] <- e
   co <- compare_networks(reference$W,W,paste(label,k,sep='_'))
   co$network <- key;co$rule <- rule;co$imputation <- k
   co$N_CC <- reference$n;co$N_PMM <- nrow(completed)
   co$scope <- 'same wave-specific preprocessing; different missing-data handling/sample'
   comparisons[[paste(label,k)]] <- co
   ct <- central(W);ct$network <- key;ct$rule <- rule;ct$imputation <- k
   all_central[[paste(label,k)]] <- ct
   if(k==1L)saveRDS(W,file.path(private,paste0('14_',label,'_first_imputation_graph.rds')))
  }
  audit[[label]] <- data.frame(modality=mod,visit=w,rule=rule,N=nrow(x),
   N_complete=sum(cc),missing_cells=sum(is.na(x)),m=20L,iterations=iterations,
   logged_events=if(is.null(events))0L else nrow(events),
   adjustment='wave-specific',scale='continuous adjusted values',
   eligible_and_QC_pass=TRUE,CC_reference_matches=TRUE)
 }
}
csv(dplyr::bind_rows(audit),'14_wave_PMM_input_audit.csv')
csv(dplyr::bind_rows(comparisons),'14_wave_PMM_vs_wave_CC_DESCRIPTIVE.csv')
csv(dplyr::bind_rows(all_edges),'14_wave_PMM_all20_edges_DESCRIPTIVE.csv')
csv(dplyr::bind_rows(all_central),'14_wave_PMM_all20_centrality_DESCRIPTIVE.csv')
if(length(logged))csv(dplyr::bind_rows(logged),'14_wave_PMM_logged_events.csv')
writeLines(c('Replaces the inherited pooled-PMM vs wave-specific complete-case comparison.',
 'All inputs are account-eligible, modality-QC-passing, selected-family records.',
 'PMM uses wave-specific continuous adjusted values, m=20, seed=500, 50/20 iterations.',
 'All 20 graphs reported descriptively; no Rubin-pooled network or NCT inference.',
 'First-imputation graph preserves the original reporting convention; it is descriptive.',
 'Review trace plots and logged events. PMM does not establish MAR or remove selection bias.'
),file.path(OUT,'14_READ_ME.txt'))
save_environment('14')
files <- list.files(OUT,pattern='^14_.*[.](csv|pdf|txt)$',full.names=TRUE)
if(!all(file.copy(files,share,overwrite=TRUE))) stop('Aggregate export failed.')
message('14 complete. Send ONLY: ',share)
