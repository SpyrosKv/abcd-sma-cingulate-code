# Paper 1 only. Read-only audit: does not modify inputs or the primary analyses.
# Put beside 00; run in the same RStudio project. Select the large imagedataset.sav.
for(p in c('haven','tidyselect')) if(!requireNamespace(p,quietly=TRUE))
 stop('Install package ',p,' first.')
root <- getwd()
out_e <- file.path(root,'Paper1_revision_outputs','02e_motion_source_audit')
dir.create(out_e,recursive=TRUE,showWarnings=FALSE)
write_e <- function(x,n) write.csv(x,file.path(out_e,n),row.names=FALSE,na='')
message('Choose the FULL MRI source imagedataset.sav (the large file), NOT Novel tech, Paper1 or the QC extract.')
p_full <- file.choose()
header_e <- haven::read_sav(p_full,n_max=0,user_na=TRUE,.name_repair='check_unique')
attrs_text <- function(x,a) {
 v <- attr(x,a); if(is.null(v)) return('')
 paste(capture.output(dput(v)),collapse=' ')
}
dictionary_e <- do.call(rbind,lapply(names(header_e),function(v) {
 x <- header_e[[v]]; data.frame(variable=v,label=if(is.null(attr(x,'label'))) '' else attr(x,'label'),
  storage=typeof(x),classes=paste(class(x),collapse=';'),
  value_labels=attrs_text(x,'labels'),SPSS_missing_values=attrs_text(x,'na_values'),
  SPSS_missing_range=attrs_text(x,'na_range'),format=attrs_text(x,'format.spss'))
}))
write_e(dictionary_e,'02e_full_source_variable_dictionary.csv')
write_e(dictionary_e[grepl('dmri|dti|scan|series|run|repet|frame|motion|qc|education',
 paste(dictionary_e$variable,dictionary_e$label),ignore.case=TRUE),],
 '02e_relevant_variable_dictionary.csv')
ids_e <- c('participant_id','session_id')
motion_e <- c('mr_y_qc__raw__dmri__mot__ok_mean','mr_y_qc__raw__dmri__mot_mean')
qc_e <- c('mr_y_qc__incl__dmri_indicator','mr_y_qc__incl__smri__t2_indicator')
required_e <- c(ids_e,motion_e,qc_e)
if(length(setdiff(required_e,names(header_e)))) stop('Wrong/incomplete source; absent: ',
 paste(setdiff(required_e,names(header_e)),collapse=', '),'. Metadata exports are available.')
roi_e <- c('gcpd','gcpv','gsca','gscma','gscmp','scm')
fa_e <- unlist(lapply(c('lh','rh'),function(h)
 paste0('mr_y_dti__is__fa__wm__dst__',roi_e,'__',h,'_mean')))
support_e <- names(header_e)[grepl('^mr_y_qc__raw__dmri__|^mr_y_qc__post__aut__dmri__|^mr_y_qc__post__man__(dmri|fsurf)|^ab_g_dyn__design_mr__',names(header_e))]
selected_e <- unique(c(required_e,support_e,intersect(fa_e,names(header_e))))
write_e(data.frame(variable=selected_e),'02e_selected_columns.csv')
message('Reading ',length(selected_e),' selected columns. The file may still take time to scan.')
full_e <- haven::read_sav(p_full,col_select=tidyselect::all_of(selected_e),
 user_na=TRUE,.name_repair='check_unique')
# SPSS user-missing values are retained in metadata, converted to NA for comparisons.
numeric_e <- function(x) as.numeric(haven::zap_labels(haven::zap_missing(x)))
decode_e <- function(z,source) {
 x <- z$session_id; raw <- haven::zap_labels(haven::zap_missing(x))
 if(is.character(raw) && all(grepl('^ses-[0-9]{2}A$',raw[!is.na(raw)]))) visit <- raw else {
  lab <- attr(x,'labels')
  if(is.null(lab)||anyDuplicated(unname(lab))||anyDuplicated(names(lab)))
   stop(source,': absent/ambiguous visit labels. Do not guess visit codes.')
  visit <- names(lab)[match(raw,unname(lab))]
 }
 if(anyNA(visit)||any(!grepl('^ses-[0-9]{2}A$',visit))) stop(source,': invalid visit mapping.')
 z$participant_id <- as.character(haven::zap_labels(haven::zap_missing(z$participant_id)))
 z$visit <- visit
 if(anyNA(z$participant_id)||any(!nzchar(trimws(z$participant_id)))||
    anyDuplicated(z[c('participant_id','visit')])) stop(source,': missing/duplicate participant-visit key.')
 write_e(unique(data.frame(source=source,stored_code=as.character(raw),visit=visit)),
  paste0('02e_',source,'_visit_mapping.csv'))
 z
}
full_e <- decode_e(full_e,'full_source')
for(v in setdiff(selected_e,ids_e)) if(is.numeric(full_e[[v]])) full_e[[v]] <- numeric_e(full_e[[v]])
message('Choose the MRI_QC_original_codes.sav extract used in 01.')
p_extract <- file.choose()
extract_e <- haven::read_sav(p_extract,col_select=tidyselect::all_of(required_e),user_na=TRUE,
 .name_repair='check_unique')
extract_e <- decode_e(extract_e,'QC_extract')
# Match each source using its OWN visit labels; never equate stored numeric codes.
key_e <- function(z) paste(z$participant_id,z$visit,sep='\r')
pos_e <- match(key_e(extract_e),key_e(full_e))
write_e(data.frame(N_full=nrow(full_e),N_extract=nrow(extract_e),
 N_extract_matched=sum(!is.na(pos_e)),N_extract_unmatched=sum(is.na(pos_e)),
 N_full_not_in_extract=sum(!key_e(full_e)%in%key_e(extract_e))), '02e_key_match_summary.csv')
compare_e <- list()
for(v in c(motion_e,qc_e)) {
 a <- numeric_e(extract_e[[v]]); b <- numeric_e(full_e[[v]][pos_e])
 matched <- !is.na(pos_e); missing_different <- matched & xor(is.na(a),is.na(b))
 valid <- matched & is.finite(a) & is.finite(b)
 diff <- abs(a[valid]-b[valid])
 compare_e[[v]] <- data.frame(variable=v,N_matched=sum(matched),
  N_finite_pairs=sum(valid),N_missing_disagreement=sum(missing_different),
  N_numeric_disagreement=sum(diff>1e-10),max_abs_difference=if(length(diff)) max(diff) else NA_real_)
}
write_e(do.call(rbind,compare_e),'02e_source_vs_QC_extract.csv')
summary_e <- correlations_e <- candidates_e <- list()
for(w in sort(unique(full_e$visit))) {
 z <- full_e[full_e$visit==w,,drop=FALSE]
 pass <- z[[qc_e[1]]] %in% 1
 for(scope in c('all_source_records','dmri_QC_pass')) {
  zz <- if(scope=='dmri_QC_pass') z[pass,,drop=FALSE] else z
  for(v in motion_e) {
   x <- zz[[v]]; finite <- is.finite(x); q <- if(any(finite))
    quantile(x[finite],c(.5,.95,.99,1),names=FALSE) else rep(NA_real_,4)
   summary_e[[length(summary_e)+1]] <- data.frame(visit=w,scope=scope,variable=v,
    N=nrow(zz),missing=sum(is.na(x)),nonfinite=sum(!is.na(x)&!is.finite(x)),
    negative=sum(x<0,na.rm=TRUE),zero=sum(x==0,na.rm=TRUE),median=q[1],p95=q[2],p99=q[3],maximum=q[4])
  }
  good <- is.finite(zz[[motion_e[1]]]) & is.finite(zz[[motion_e[2]]])
  a <- zz[[motion_e[1]]][good]; b <- zz[[motion_e[2]]][good]
  correlations_e[[length(correlations_e)+1]] <- data.frame(visit=w,scope=scope,N_pairs=length(a),
   Pearson=if(length(a)>2 && sd(a)>0 && sd(b)>0) cor(a,b) else NA_real_,
   Spearman=if(length(a)>2 && sd(a)>0 && sd(b)>0) cor(a,b,method='spearman') else NA_real_,
   N_equal=sum(abs(a-b)<=1e-10),max_abs_difference=if(length(a)) max(abs(a-b)) else NA_real_)
 }
 x <- z[[motion_e[1]]]
 # Review zeros plus top 10 QC-passing records; this is NOT an exclusion cutoff.
 top <- which(pass & is.finite(x));top <- top[order(-x[top],z$participant_id[top])]
 review <- sort(unique(c(which(pass & x %in% 0),head(top,10))))
 if(length(review)) candidates_e[[length(candidates_e)+1]] <- z[review,,drop=FALSE]
}
write_e(do.call(rbind,summary_e),'02e_motion_summary_by_visit.csv')
write_e(do.call(rbind,correlations_e),'02e_motion_measure_agreement.csv')
saveRDS(full_e,file.path(out_e,'02e_selected_source_columns_RESTRICTED.rds'))
if(length(candidates_e)) write_e(do.call(rbind,candidates_e),'02e_zero_and_top_motion_with_supporting_fields_RESTRICTED.csv')
# Aggregate supporting values among zero records; no participant identifiers.
zero_audit_e <- list()
for(w in sort(unique(full_e$visit))) {
 z <- full_e[full_e$visit==w & full_e[[qc_e[1]]]%in%1 & full_e[[motion_e[1]]]%in%0,,drop=FALSE]
 for(v in setdiff(support_e,ids_e)) {
  x <- z[[v]]; if(!length(x)) next
  if(is.numeric(x)) zero_audit_e[[length(zero_audit_e)+1]] <- data.frame(visit=w,variable=v,N=nrow(z),
   missing=sum(is.na(x)),minimum=if(any(is.finite(x))) min(x[is.finite(x)]) else NA_real_,
   maximum=if(any(is.finite(x))) max(x[is.finite(x)]) else NA_real_)
 }
}
write_e(if(length(zero_audit_e)) do.call(rbind,zero_audit_e) else data.frame(),
 '02e_zero_motion_supporting_numeric_summary.csv')
writeLines(c(paste('Full-source filename:',basename(p_full)),
 paste('Full-source size bytes:',file.info(p_full)$size),paste('QC extract filename:',basename(p_extract)),
 paste('Full-source modified time:',file.info(p_full)$mtime),
 'File timestamps do not establish ABCD release or download provenance.',
 'Audit preserves source values and all original analysis outputs. No new exclusions or primary refits.',
 'Variable labels and SPSS missing codes are embedded metadata, not the complete official data dictionary.',
 'Source/extract equality verifies copying and joining, not whether zero/extreme motion is biologically or technically valid.',
 'Correlation or equal values between motion summaries does not establish correspondence to the FA-generating scan set.',
 'FA-generating scan correspondence requires explicit release documentation or series/processing provenance.',
 'No units are inferred from numeric magnitudes. No automatic declaration that motion is confirmed.',
 'Do not upload files marked RESTRICTED. Keep them in the authorised research environment.'),
 file.path(out_e,'02e_source_audit_limits.txt'))
capture.output(sessionInfo(),file=file.path(out_e,'02e_R_session_info.txt'))
message('02e complete: ',normalizePath(out_e),'. Send CSV/TXT except RESTRICTED. Source-scan correspondence is not automatically confirmed.')
