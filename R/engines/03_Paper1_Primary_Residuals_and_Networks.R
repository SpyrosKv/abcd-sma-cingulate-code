source('00_Paper1_Settings_and_Functions.R')
d<-read_stage('01_merged_RESTRICTED.rds')
if(!CFG$motion_source_confirmed)message('Motion covariate is provisional: confirm correspondence to FA source scans before final submission.')
res<-list(); nets<-list(); diag<-list();et<-list();ct<-list()
for(mod in names(MRI)) {
 f<-flags_for(d,mod); z<-droplevels(d[f$base,,drop=FALSE])
 # Fit on available outcomes within the QC/eligibility/family/covariate-defined population.
 rr<-residualise(z,mod,scope='pooled',family_model=TRUE,extraction='retain_family')
 diag[[mod]]<-rr$diagnostics
 res[[mod]]<-list(source=z,residuals=rr$data,models=rr$models)
 for(w in WAVES) {
  use<-z$visit==w & complete.cases(rr$data)
  x<-rr$data[use,,drop=FALSE]; key<-paste(mod,w,sep='_'); W<-network(x)
  nets[[key]]<-list(W=W,x=x,ids=z$participant_id[use],mod=mod,visit=w,n=nrow(x))
  e<-edges(W);e$network<-key;e$N<-nrow(x);et[[key]]<-e
  ctab<-central(W);ctab$network<-key;ct[[key]]<-ctab
 }
}
saveRDS(res,file.path(OUT,'03_primary_residuals_RESTRICTED.rds'))
saveRDS(nets,file.path(OUT,'03_primary_networks_RESTRICTED.rds'))
csv(dplyr::bind_rows(diag),'03_residualisation_diagnostics.csv')
csv(dplyr::bind_rows(et),'03_all_network_edges.csv');csv(dplyr::bind_rows(ct),'03_descriptive_centrality.csv')
# Fixed layout and edge scale across six networks; graph is statistical, not anatomical.
coords<-rbind(cbind(rep(-1,6),seq(1,-1,length.out=6)),
              cbind(rep(1,12),seq(1,-1,length.out=12)))
maxedge<-max(vapply(nets,function(n)max(abs(n$W)),numeric(1)))
pdf(file.path(OUT,'03_network_figures.pdf'),width=12,height=9)
for(key in names(nets)) {
 n<-nets[[key]]
 qgraph::qgraph(n$W,layout=coords,labels=NODES,groups=list(SMA=1:6,MRI=7:18),
   color=c('#E7AE55','#70A5CC'),maximum=maxedge,minimum=0,cut=0,
   vsize=4,label.cex=.6,title=paste(key,'N =',n$n),legend=FALSE)
}
dev.off()
csv(data.frame(node=NODES,source_FA=c(SMA,MRI$FA),source_T2w=c(SMA,MRI$T2w)), '03_node_dictionary.csv')
save_environment('03');message('03 complete. These are revised cross-sectional networks; no within-person change is estimated.')
