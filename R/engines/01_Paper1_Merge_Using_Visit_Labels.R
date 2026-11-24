source('00_Paper1_Settings_and_Functions.R')
message('Choose ORIGINAL Paper1 dataset (paper1.sav or extracted dataset.sav), NOT the MRI QC source.')
p_raw <- file.choose()
message('Choose MRI_QC_original_codes.sav')
p_mri <- file.choose()
message('Choose SMA_eligibility_original_codes.sav')
p_sma <- file.choose()
if(length(unique(normalizePath(c(p_raw,p_mri,p_sma))))!=3) stop('Choose three different files')
raw_vars <- c('participant_id','session_id',AGE,FAMILY,SOFTWARE,SCANNER,GATE,DEMOS,SMA,unlist(MRI))
r <- decode_visit(read_sav_source(p_raw,unique(raw_vars)),'Paper1')
m <- decode_visit(read_sav_source(p_mri,c('participant_id','session_id',unname(QC),MOTION,MANUAL)),'MRI_QC')
s <- decode_visit(read_sav_source(p_sma,c('participant_id','session_id',GATE)),'SMA_eligibility')
csv(dplyr::bind_rows(attr(r,'visit_mapping'),attr(m,'visit_mapping'),attr(s,'visit_mapping')),'01_visit_mapping.csv')
# Keep metadata for scale-category verification and coding audits.
metadata <- lapply(list(raw=r,mri=m,sma=s),function(d) lapply(d,attributes))
saveRDS(metadata,file.path(OUT,'01_variable_metadata.rds'))
r <- r[r$visit %in% WAVES,,drop=FALSE]
for(v in c(AGE,GATE,SMA,unlist(MRI))) r[[v]] <- num(r[[v]])
for(v in c(unname(QC),MOTION,MANUAL)) m[[v]] <- num(m[[v]])
s[[GATE]] <- num(s[[GATE]])
for(d in list(r,s)) if(any(!is.na(d[[GATE]]) & !d[[GATE]] %in% c(0,1))) stop('Unexpected gate code')
for(v in c(unname(QC),MANUAL)) if(any(!is.na(m[[v]]) & !m[[v]] %in% c(0,1))) stop('Unexpected QC code: ',v)
# No inner join: preserve the original source cohort and audit unmatched records.
names(s)[names(s)==GATE] <- 'account_new'; s$in_sma_extract <- TRUE; m$in_mri_extract <- TRUE
m$session_id <- NULL; s$session_id <- NULL
j <- dplyr::left_join(r,m,by=c('participant_id','visit'),relationship='one-to-one')
j <- dplyr::left_join(j,s,by=c('participant_id','visit'),relationship='one-to-one')
stopifnot(nrow(j)==nrow(r)); check_keys(j,'Merged')
conflict <- !is.na(j[[GATE]]) & !is.na(j$account_new) & j[[GATE]]!=j$account_new
csv(j[conflict,c('participant_id','visit',GATE,'account_new')],'01_account_conflicts_RESTRICTED.csv')
if(any(conflict)) stop('Conflicting observed account responses; inspect 01_account_conflicts before proceeding')
# Primary gate remains the ORIGINAL response; new extract is validation, not silent backfill.
j$age <- j[[AGE]] # original units retained; scaling does not change a linear age adjustment
j$family <- clean_factor(j[[FAMILY]]); j$scanner <- clean_factor(j[[SCANNER]])
j$software <- clean_factor(j[[SOFTWARE]])
# One fixed randomly selected participant per family across ALL target waves.
fam <- unique(j[c('participant_id','family')]); fam <- fam[!is.na(fam$family),]
if(anyDuplicated(fam$participant_id)) stop('Participant has multiple family IDs; resolve before sampling')
fam <- fam[order(as.character(fam$family),fam$participant_id),]
set.seed(CFG$seed)
selected <- unlist(lapply(split(fam$participant_id,fam$family,drop=TRUE),function(ids) sample(sort(ids),1)),use.names=FALSE)
j$selected_family <- j$participant_id %in% selected
csv(fam,'01_family_membership_RESTRICTED.csv')
csv(data.frame(participant_id=selected),'01_selected_participants_RESTRICTED.csv')
# Check six-point SMA coding against source labels; no arbitrary recoding.
category_info <- list()
for(v in SMA) {
 labs <- metadata$raw[[v]]$labels
 # Documented six response categories are 1..6; 777 = Decline to answer.
 levels <- sort(unique(as.numeric(labs)[as.numeric(labs) %in% 1:6]))
 if(any(j[[v]] %in% 777)) j[[v]][j[[v]] %in% 777] <- NA_real_
 observed <- sort(unique(j[[v]][!is.na(j[[v]])]))
 if(length(levels)!=6 || any(!observed %in% levels)) stop('Verify six SMA category codes in metadata: ',v)
 category_info[[v]] <- levels
}
saveRDS(category_info,file.path(OUT,'01_sma_categories.rds'))
contrad <- !is.na(j[[GATE]]) & j[[GATE]]==0 & rowSums(!is.na(j[SMA]))>0
csv(j[contrad,c('participant_id','visit',GATE,SMA)],'01_no_account_with_SMA_RESTRICTED.csv')
csv(dplyr::bind_rows(lapply(WAVES,function(w) {
 z<-j[j$visit==w,]; data.frame(visit=w,N=nrow(z),MRI_record_matched=sum(z$in_mri_extract %in% TRUE),
 SMA_record_matched=sum(z$in_sma_extract %in% TRUE),new_gate_only=sum(is.na(z[[GATE]]) & !is.na(z$account_new)),
 original_gate_only=sum(!is.na(z[[GATE]]) & is.na(z$account_new)),
 no_account_with_SMA=sum(!is.na(z[[GATE]]) & z[[GATE]]==0 & rowSums(!is.na(z[SMA]))>0))
})),'01_join_audit.csv')
saveRDS(j,file.path(OUT,'01_merged_RESTRICTED.rds'))
csv(data.frame(source=c('Paper1','MRI_QC','SMA'),path=c(p_raw,p_mri,p_sma),
 md5=unname(tools::md5sum(c(p_raw,p_mri,p_sma)))),'01_input_provenance_RESTRICTED.csv')
save_environment('01'); message('01 complete: inspect 01_join_audit.csv and 01_visit_mapping.csv')
