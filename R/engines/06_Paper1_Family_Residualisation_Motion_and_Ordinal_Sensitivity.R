source('00_Paper1_Settings_and_Functions.R')
d<-read_stage('01_merged_RESTRICTED.rds');primary<-read_stage('03_primary_networks_RESTRICTED.rds')
res_primary<-read_stage('03_primary_residuals_RESTRICTED.rds');cats<-read_stage('01_sma_categories.rds')
comparison<-diagnostics<-all_edges<-ordinal_pd<-list()
# Change one analytical component at a time. Complete-case comparisons use common IDs.
variants<-list(
 no_family_random_intercept=list(scope='pooled',family_model=FALSE,extraction='retain_family',one=TRUE,motion=CFG$motion_var),
 wave_specific_no_family=list(scope='wave',family_model=FALSE,extraction='retain_family',one=TRUE,motion=CFG$motion_var),
 conditional_family_removed_DIAGNOSTIC=list(scope='pooled',family_model=TRUE,extraction='conditional',one=TRUE,motion=CFG$motion_var),
 all_families_retained_DESCRIPTIVE=list(scope='pooled',family_model=TRUE,extraction='retain_family',one=FALSE,motion=CFG$motion_var),
 alternate_motion_common_sample=list(scope='pooled',family_model=TRUE,extraction='retain_family',one=TRUE,motion=MOTION[2])
)
for(name in names(variants)) for(mod in names(MRI)) {
 if(name=='alternate_motion_common_sample' && mod!='FA')next
 v<-variants[[name]];f<-flags_for(d,mod,motion=v$motion,one_family=v$one)
 take<-f$base
 if(name=='alternate_motion_common_sample')take<-take & flags_for(d,mod,motion=CFG$motion_var)$base
 z<-droplevels(d[take,,drop=FALSE]);message(name,': ',mod)
 rr<-residualise(z,mod,scope=v$scope,family_model=v$family_model,extraction=v$extraction,motion=v$motion)
 dg<-rr$diagnostics;dg$variant<-name;diagnostics[[length(diagnostics)+1]]<-dg
 # For alternate motion, refit the reference on the same base population too.
 reference<-if(name=='alternate_motion_common_sample') residualise(z,mod,motion=CFG$motion_var) else NULL
 for(w in WAVES) {
  key<-paste(mod,w,sep='_');n<-primary[[key]]
  use<-z$visit==w & complete.cases(rr$data);x<-rr$data[use,,drop=FALSE];ids<-z$participant_id[use]
  if(name=='all_families_retained_DESCRIPTIVE') {
   # Different sample by design. No independent-observation CI or p-value is calculated here.
   A<-n$W;B<-network(x);scope_label<-'different_sample_family_sensitivity'
   nref<-n$n
  } else {
   if(!is.null(reference)) {
    ref_use<-z$visit==w & complete.cases(reference$data)
    ref_x<-reference$data[ref_use,,drop=FALSE];ref_ids<-z$participant_id[ref_use]
   } else {ref_x<-n$x;ref_ids<-n$ids}
   common<-sort(intersect(ids,ref_ids));x<-x[match(common,ids),,drop=FALSE]
   A<-network(ref_x[match(common,ref_ids),,drop=FALSE]);B<-network(x)
   scope_label<-'matched_participants';nref<-length(common)
  }
  co<-compare_networks(A,B,paste(name,key));co$N_reference<-nref;co$N_sensitivity<-nrow(x);co$scope<-scope_label
  comparison[[length(comparison)+1]]<-co
  ee<-edges(B);ee$variant<-name;ee$network<-key;all_edges[[length(all_edges)+1]]<-ee
 }
}
# Ordinal sensitivity: original ordered SMA items (NEVER ordinalise continuous residuals).
# MRI values are the primary adjusted values; age is partialled from the latent correlation
# matrix via its Schur complement. The matched Pearson benchmark uses the same construction.
for(mod in names(MRI)) for(w in WAVES) {
 key<-paste(mod,w,sep='_');p<-primary[[key]];rs<-res_primary[[mod]]
 ix<-match(p$ids,rs$source$participant_id[rs$source$visit==w])
 rows<-which(rs$source$visit==w)[ix];z<-rs$source[rows,,drop=FALSE]
 X<-data.frame(z[SMA],rs$residuals[rows,7:18,drop=FALSE],age=z$age,check.names=FALSE)
 names(X)[1:18]<-NODES
 continuous<-X
 for(k in 1:6)X[[k]]<-ordered(X[[k]],levels=cats[[SMA[k]]])
 hc<-polycor::hetcor(X,ML=FALSE,std.err=FALSE,use='complete.obs',pd=FALSE)
 R<-hc$correlations
 if(any(!is.finite(R)))stop('Ordinal correlation failed: ',key)
 min_eig<-min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)
 # No silent nearest-PD repair: an inadmissible matrix requires investigation.
 if(min_eig<=1e-8)stop('Ordinal matrix is not positive definite: ',key,'. Inspect sparse categories; do not silently smooth.')
 adjust_age<-function(R) {
  S<-R[1:18,1:18]-R[1:18,19,drop=FALSE]%*%solve(R[19,19,drop=FALSE])%*%R[19,1:18,drop=FALSE]
  cov2cor(S)
 }
 W_ord<-qgraph::EBICglasso(adjust_age(R),nrow(X),gamma=CFG$gamma,verbose=FALSE)
 W_cont<-qgraph::EBICglasso(adjust_age(cor(continuous)),nrow(X),gamma=CFG$gamma,verbose=FALSE)
 ee<-edges(W_ord);ee$variant<-'ordinal_latent_age_adjusted';ee$network<-key;all_edges[[length(all_edges)+1]]<-ee
 co<-compare_networks(W_cont,W_ord,paste('ordinal_vs_matched_Pearson',key));co$N_reference<-nrow(X);co$N_sensitivity<-nrow(X);co$scope<-'same_sample_same_age_adjustment_different_correlation'
 comparison[[length(comparison)+1]]<-co
 co<-compare_networks(p$W,W_ord,paste('ordinal_vs_primary',key));co$N_reference<-nrow(X);co$N_sensitivity<-nrow(X);co$scope<-'same_sample_different_SMA_adjustment_estimand'
 comparison[[length(comparison)+1]]<-co
 ordinal_pd[[key]]<-data.frame(network=key,N=nrow(X),min_eigenvalue=min_eig)
 saveRDS(list(correlation=R,ordinal=W_ord,matched_Pearson=W_cont),file.path(OUT,paste0('06_',key,'_ordinal.rds')))
}
csv(dplyr::bind_rows(comparison),'06_sensitivity_network_comparisons.csv')
csv(dplyr::bind_rows(diagnostics),'06_sensitivity_residualisation_diagnostics.csv')
csv(dplyr::bind_rows(all_edges),'06_sensitivity_edges.csv')
csv(dplyr::bind_rows(ordinal_pd),'06_ordinal_matrix_diagnostics.csv')
save_environment('06');message('06 complete: ordinal and all-family analyses are sensitivity analyses, not replacements for primary inference.')
