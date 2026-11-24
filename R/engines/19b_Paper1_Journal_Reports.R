# 19b: Generate journal reports after successful completion of 19a.
# Generate journal reports from completed revised analyses.
# No model fitting, imputation, network estimation or statistical testing is performed.
options(stringsAsFactors = FALSE)
packages <- c('xml2', 'base64enc', 'ggplot2', 'qgraph', 'pdftools')
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_packages)) stop('Install these packages first: ', paste(missing_packages, collapse=', '))
message('Choose 03_all_network_edges.csv from the FINAL revised results, not the original results.')
p_edges <- normalizePath(file.choose(), winslash='/', mustWork=TRUE)
if(basename(p_edges) != '03_all_network_edges.csv') stop('Select 03_all_network_edges.csv.')
message('Choose the ORIGINAL Supplementary File S2 HTML.')
p_s2 <- normalizePath(file.choose(), winslash='/', mustWork=TRUE)
message('Choose the ORIGINAL Supplementary File S3 HTML.')
p_s3 <- normalizePath(file.choose(), winslash='/', mustWork=TRUE)
message('Choose the NEW 19a_retention_model_summary.csv from SHAREABLE_RETENTION.')
p_retention <- normalizePath(file.choose(),winslash='/',mustWork=TRUE)
if(basename(p_retention)!='19a_retention_model_summary.csv') stop('Run 19a first and select its model summary.')
RETENTION <- dirname(p_retention)
if(!file.exists(file.path(RETENTION,'19a_complete.csv'))) stop('19a has not completed successfully.')
completion <- read.csv(file.path(RETENTION,'19a_complete.csv'))
paths <- file.path(RETENTION,completion$file)
if(!all(file.exists(paths)) || !identical(unname(tools::md5sum(paths)),completion$md5))
 stop('Retention outputs differ from the completed 19a export.')
provenance <- read.csv(file.path(RETENTION,'19a_input_provenance.csv'))
if(!identical(unname(tools::md5sum(p_edges)),provenance$md5[provenance$file=='03_all_network_edges.csv']))
 stop('Retention and report inputs come from different primary results.')

# A consolidated results package is preferred. The original completed output
# folders are also supported if they sit beside the selected results folder.
parent <- dirname(p_edges)
if(basename(parent) %in% c('01_Primary_and_sensitivity_results','SHAREABLE_SUMMARIES')) parent <- dirname(parent)
code_parent <- dirname(parent)
roots <- unique(c(parent,RETENTION,
 file.path(code_parent, 'Paper1_revision_outputs_original_methods_v2'),
 file.path(code_parent, 'Paper1_final_computational_checks'),
 file.path(code_parent, 'Paper1_original_required_checks'),
 file.path(code_parent, 'Paper1_VAS_remaining_checks')))
roots <- roots[dir.exists(roots)]
all_files <- unique(unlist(lapply(roots, list.files, recursive=TRUE, full.names=TRUE)))
all_files <- all_files[!grepl('/Journal_Supplementary_Reports/', gsub('\\\\','/',all_files), fixed=TRUE)]
find_file <- function(name) {
 aliases <- c('13_equivalence_TOST_two_bounds.csv'='13_equivalence_TOST_both_VAS_bounds.csv')
 accepted <- c(name, unname(aliases[name]))
 accepted <- accepted[!is.na(accepted)]
 hits <- if(startsWith(name,'19a_')) file.path(RETENTION,name) else all_files[basename(all_files) %in% accepted]
 hits <- hits[file.exists(hits)]
 if(!length(hits)) stop('Required result is absent: ',name,
  '. Use the consolidated final results package, or place all completed output folders together.')
 if(length(hits)>1L && length(unique(unname(tools::md5sum(hits))))>1L)
  stop('Different versions found for ',name,'. Use one consolidated final results package.')
 hits[1L]
}
read_result <- function(name) read.csv(find_file(name), check.names=FALSE, na.strings=c('NA','NaN',''))
required <- c('03_all_network_edges.csv','03_descriptive_centrality.csv',
 '03_node_dictionary.csv','03_residualisation_diagnostics.csv',
 '02_sequential_sample_flow.csv','02_structural_vs_item_missingness.csv',
 '02b_model_sample_and_convergence.csv','02b_bias_reduced_logistic.csv',
 '02_account_and_visit_transitions.csv','02_motion_audit.csv',
 '04_bootnet_centrality_CS.csv','04_bootnet_edge_intervals.csv',
 '05_original_NCT_omnibus.csv','05_original_NCT_edges.csv',
 '06_sensitivity_network_comparisons.csv','06_ordinal_matrix_diagnostics.csv',
 '07_original_PMM_eligibility_audit.csv','07_original_PMM_vs_primary_DESCRIPTIVE.csv',
 '12_Little_MCAR_eligible_AND.csv','12_primary_residual_reconstruction.csv',
 '12_wave_specific_edges.csv','12_wave_specific_model_diagnostics.csv',
 '12_ORIGINAL_posthoc_regression_sensitivity.csv',
 '13_primary_and_sensitivity_centrality.csv','13_sensitivity_NCT_omnibus.csv',
 '13_sensitivity_NCT_edges.csv','13_wave_specific_centrality_CS.csv',
 '13_wave_specific_edge_intervals.csv','13_unregularised_432_partial_correlations.csv',
 '13_equivalence_summary.csv','13_equivalence_TOST_two_bounds.csv',
 '14_wave_PMM_input_audit.csv','14_wave_PMM_vs_wave_CC_DESCRIPTIVE.csv',
 '15_parental_education_selection_SMD.csv','15_education_join_and_observed_counts.csv',
 '16_domain_strength_NCT_18_tests.csv','17_family_selection_audit.csv',
 '17_family_selection_network_comparisons.csv','17_family_selection_NCT_omnibus.csv',
 '18_SMA_distributions_with_kurtosis.csv','18_disjoint_exclusion_group_counts.csv',
 '18_included_vs_each_exclusion_group_SMD.csv','18_SMA_category_frequencies.csv',
 '17_family_selection_centrality_CS.csv','17_family_selection_edge_intervals.csv',
 '19a_retention_model_summary.csv','19a_retention_logistic_coefficients.csv',
 '19a_retention_wave_availability.csv','19a_retention_patterns.csv',
 '19a_retention_outcome_counts.csv','19a_retained_vs_not_retained_descriptives.csv',
 '19a_predictor_source_visits.csv','19a_retention_reference_categories.csv','19a_retention_definitions.csv')
D <- setNames(lapply(required, read_result), required)
get_result <- function(name) D[[name]]
primary_edges <- get_result('03_all_network_edges.csv')
wave_edges <- get_result('12_wave_specific_edges.csv')
nodes <- get_result('03_node_dictionary.csv')$node
keys <- unique(primary_edges$network)
if(length(nodes)!=18L || anyDuplicated(nodes) || length(keys)!=6L)
 stop('Expected 18 distinct nodes and six networks.')
validate_edges <- function(e) {
 if(nrow(e)!=918L || !setequal(unique(e$network),keys) || anyNA(e$weight) ||
    any(!is.finite(e$weight)) || any(table(e$network)!=153L) ||
    anyDuplicated(e[c('network','i','j')])) stop('Invalid network edge table.')
 if(any(e$i<1 | e$j>18 | e$i>=e$j) ||
    !all(e$node1==nodes[e$i]) || !all(e$node2==nodes[e$j])) stop('Node order mismatch.')
}
validate_edges(primary_edges); validate_edges(wave_edges)
models <- get_result('02b_model_sample_and_convergence.csv')
if(anyNA(models$converged) || !all(models$converged)) stop('Selection models did not all converge.')
diag_primary <- get_result('03_residualisation_diagnostics.csv')
if(any(grepl('visit|age.*:', diag_primary$formula))) stop('Unexpected primary adjustment formula.')
if(!all(get_result('12_primary_residual_reconstruction.csv')$retained_family) ||
 max(get_result('12_primary_residual_reconstruction.csv')$maximum_reconstruction_error)>1e-7)
 stop('Residual reconstruction check failed.')
if(!all(get_result('05_original_NCT_omnibus.csv')$permutations==1000L)) stop('Unexpected NCT permutation count.')
retention_models <- get_result('19a_retention_model_summary.csv')
if(nrow(retention_models)!=3L || !all(retention_models$converged)) stop('Three converged retention models are required.')
# Resolve every bootstrap figure before writing the reports.
boot_files <- setNames(lapply(keys, function(k) c(
 primary_edge=find_file(paste0('04_',k,'_edge_accuracy.pdf')),
 primary_case=find_file(paste0('04_',k,'_case_stability.pdf')),
 wave_edge=find_file(paste0('13_',k,'_wave_edge_accuracy.pdf')),
 wave_case=find_file(paste0('13_',k,'_wave_case_stability.pdf')))),keys)
check_template <- function(path, expected) {
 doc <- xml2::read_html(path)
 title <- xml2::xml_text(xml2::xml_find_first(doc,'//h1[contains(@class,"title")]'))
 if(!grepl(expected,title,fixed=TRUE)) stop('Wrong HTML template selected: ',expected)
 if(length(xml2::xml_find_all(doc,'//div[contains(@class,"toc-content")]'))!=1L)
  stop('The original HTML content container was not found.')
 doc
}
templates <- list(S2=check_template(p_s2,'Supplementary File S2'),
                  S3=check_template(p_s3,'Supplementary File S3'))
OUT <- file.path(parent, 'Journal_Supplementary_Reports')
dir.create(OUT, showWarnings=FALSE)
FIG <- file.path(OUT, 'S1_Figures'); dir.create(FIG, showWarnings=FALSE)

esc <- function(x) {
 x <- as.character(x); x[is.na(x)] <- '\u2014'
 x <- gsub('&','&amp;',x,fixed=TRUE); x <- gsub('<','&lt;',x,fixed=TRUE)
 x <- gsub('>','&gt;',x,fixed=TRUE); x <- gsub('"','&quot;',x,fixed=TRUE); x
}
pretty <- function(x) {
 if(is.numeric(x)) {
  y <- rep('\u2014',length(x)); ok <- is.finite(x)
  y[ok] <- format(x[ok],digits=6,trim=TRUE,scientific=FALSE); return(y)
 }
 as.character(x)
}
report_labels <- function(x) {
 replacements <- c(
 'conditional_family_removed_DIAGNOSTIC'='Conditional residuals with family contribution removed',
 'all_families_retained_DESCRIPTIVE'='All families retained (descriptive sensitivity)',
 'no_family_random_intercept'='Pooled adjustment without family random intercept',
 'wave_specific_no_family'='Wave-specific adjustment',
 'alternate_motion_common_sample'='Alternative motion covariate on common participants',
 'ordinal_vs_matched_Pearson'='Ordinal versus matched Pearson correlations',
 'ordinal_vs_primary'='Ordinal versus primary pooled inputs',
 'same_sample_different_SMA_adjustment_estimand'='Same participants; different SMA adjustment',
 'same_sample_same_age_adjustment_different_correlation'='Same participants and age adjustment; different correlation estimator',
 'different_sample_family_sensitivity'='Family sensitivity with a different participant sample',
 'matched_participants'='Matched participants',
 'wave_specific'='Wave-specific', 'no_family_pooled'='Pooled without family intercept',
 'retain_family'='Family contribution retained',
 'SMA_nonresponse'='SMA item nonresponse', 'FA_excluded'='Excluded from FA network',
 'T2w_excluded'='Excluded from T2-weighted network',
 'selected_family_account_holders'='Account holders in the selected family sample',
 'all_account_holders'='All account holders', 'ALL_SIX'='All six networks')
 for(key in names(replacements)) x <- gsub(key,replacements[[key]],x,fixed=TRUE)
 x
}
tab <- function(d, caption) {
 if(!is.data.frame(d) || !nrow(d)) stop('Empty report table: ',caption)
 names_display <- gsub('_',' ',names(d),fixed=TRUE)
 vals <- lapply(d, function(x) esc(report_labels(pretty(x))))
 rows <- vapply(seq_len(nrow(d)), function(i)
  paste0('<tr>',paste0('<td>',vapply(vals,`[`,character(1),i),'</td>',collapse=''),'</tr>'),character(1))
 paste0('<div style="overflow-x:auto"><table class="table table-condensed"><caption>',esc(caption),'</caption><thead><tr>',
  paste0('<th>',esc(names_display),'</th>',collapse=''),'</tr></thead><tbody>',
  paste(rows,collapse='\n'),'</tbody></table></div>')
}
p <- function(x) paste0('<p>',esc(x),'</p>')
section <- function(n, title, text) paste0('<div id="section-',n,'" class="section level1">',
 '<h1><span class="header-section-number">',n,'</span> ',esc(title),'</h1>',text,'</div>')
image_tag <- function(path,caption) paste0('<figure><img src="data:image/png;base64,',
 base64enc::base64encode(path),'" alt="',esc(caption),'" style="width:95%">',
 '<figcaption>',esc(caption),'</figcaption></figure>')
W_from_edges <- function(e,key) {
 z <- e[e$network==key,]; W <- matrix(0,18,18,dimnames=list(nodes,nodes))
 W[cbind(z$i,z$j)] <- z$weight; W <- W+t(W); W
}
# A single spring layout and edge scale across primary/sensitivity figures.
# This is a presentation choice; stored estimates are not changed.
Ws <- lapply(keys,function(k) W_from_edges(primary_edges,k))
set.seed(20261005)
layout <- qgraph::qgraph(Reduce(`+`,lapply(Ws,abs))/length(Ws),
 layout='spring',DoNotPlot=TRUE)$layout
maximum <- max(abs(c(primary_edges$weight,wave_edges$weight)))
if(!is.finite(maximum) || maximum<=0) stop('No valid nonzero network edges.')
central_primary <- get_result('03_descriptive_centrality.csv')
central_wave <- get_result('13_primary_and_sensitivity_centrality.csv')
central_wave <- central_wave[central_wave$variant=='wave_specific',]
if(nrow(central_wave)!=108L) stop('Expected 108 wave-specific centrality records.')
plot_central <- function(z,path) {
 z <- z[match(nodes,z$node),]
 long <- rbind(data.frame(node=z$node,value=z$strength,measure='Strength'),
  data.frame(node=z$node,value=z$expected_influence,measure='Expected influence'))
 long$node <- factor(long$node,levels=rev(nodes))
 plot <- ggplot2::ggplot(long,ggplot2::aes(x=value,y=node,group=measure))+
  ggplot2::geom_path(colour='#333333',linewidth=.35)+ggplot2::geom_point(size=1.4)+
  ggplot2::facet_wrap(~measure,scales='free_x')+ggplot2::theme_bw(base_size=10)+
  ggplot2::labs(x='Centrality (original scale)',y=NULL)
 ggplot2::ggsave(path,plot,width=8,height=6,dpi=300)
}
for(kind in c('primary','wave')) for(k in keys) {
 e <- if(kind=='primary') primary_edges else wave_edges
 ctab <- if(kind=='primary') central_primary else central_wave
 plot_central(ctab[ctab$network==k,],file.path(FIG,paste0(kind,'_',k,'_centrality.png')))
 gr <- function() qgraph::qgraph(W_from_edges(e,k),layout=layout,labels=nodes,
  groups=list(SMA=1:6,MRI=7:18),color=c('#E7AE55','#70A5CC'),
  posCol='#4575B4',negCol='#D73027',maximum=maximum,minimum=0,cut=0,
  vsize=5,label.cex=.7,title=gsub('_',' ',k),legend=FALSE)
 png(file.path(FIG,paste0(kind,'_',k,'_network.png')),width=2400,height=2400,res=300)
 gr(); dev.off()
 pdf(file.path(FIG,paste0(kind,'_',k,'_network.pdf')),width=8,height=8)
 gr(); dev.off()
 for(type in c('edge','case')) {
  src <- boot_files[[k]][paste(kind,type,sep='_')]
  destination <- file.path(FIG,paste0(kind,'_',k,'_',if(type=='edge') 'edge_accuracy' else 'case_stability','.png'))
  converted <- pdftools::pdf_convert(src,format='png',pages=1L,dpi=300,filenames=destination,verbose=FALSE)
  if(!file.exists(destination)) stop('Bootstrap image conversion failed.')
 }
}
# The process diagram is drawn separately; it is not an MRI result.
# Analytical workflow: family selection precedes eligibility and QC.
# Regenerates Figure 4 only; no statistical analyses are rerun.

local({

  # During the full report build, use its existing figure directory.
  # When run separately, select the previous Figure 4 PNG.
  if (exists("FIG", inherits = TRUE) &&
      dir.exists(get("FIG", inherits = TRUE))) {

    figure_dir <- get("FIG", inherits = TRUE)

  } else {

    message(
      "Select Figure_4_Analytical_Workflow.png ",
      "from the S1_Figures folder."
    )

    selected_figure <- normalizePath(
      file.choose(),
      winslash = "/",
      mustWork = TRUE
    )

    if (basename(selected_figure) !=
        "Figure_4_Analytical_Workflow.png") {
      stop("Please select Figure_4_Analytical_Workflow.png.")
    }

    figure_dir <- dirname(selected_figure)
  }

  workflow <- c(
    "ABCD observations at ses-02A, ses-04A and ses-06A",
    "Match participant ID and visit label across source files",
    "Select one participant per family, consistently across visits",
    "Restrict SMA eligibility to account holders at that visit",
    "Apply the FA and T2 inclusion flags separately",
    "Adjust outcomes; retain family contribution in pooled inputs",
    "Estimate separate complete-case EBICglasso networks",
    "Assess bootstrap uncertainty and paired network comparisons",
    "Assess missingness, ordinal, family-selection and PMM sensitivity"
  )

  draw_workflow <- function(path, format) {

    if (format == "png") {
      grDevices::png(
        path,
        width = 2400,
        height = 3300,
        res = 300
      )
    } else {
      grDevices::pdf(
        path,
        width = 8,
        height = 11
      )
    }

    on.exit(grDevices::dev.off(), add = TRUE)

    grid::grid.newpage()

    grid::grid.text(
      "Overview of the analytical workflow",
      y = .975,
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )

    for (i in seq_along(workflow)) {

      y <- .91 - (i - 1) * .098

      grid::grid.roundrect(
        x = .5,
        y = y,
        width = .9,
        height = .072,
        gp = grid::gpar(
          fill = c(
            "#F2E6C9", "#DFEAF3", "#DBEADB"
          )[(i - 1) %% 3 + 1],
          col = "#777777"
        )
      )

      grid::grid.text(
        paste(
          strwrap(
            paste0(i, ". ", workflow[i]),
            width = 62
          ),
          collapse = "\n"
        ),
        x = .5,
        y = y,
        gp = grid::gpar(fontsize = 11)
      )

      if (i < length(workflow)) {
        grid::grid.lines(
          x = c(.5, .5),
          y = c(y - .037, y - .058),
          arrow = grid::arrow(
            length = grid::unit(.08, "inches")
          ),
          gp = grid::gpar(col = "#777777")
        )
      }
    }
  }

  for (format in c("png", "pdf")) {

    output_path <- file.path(
      figure_dir,
      paste0("Figure_4_Analytical_Workflow.", format)
    )

    draw_workflow(output_path, format)

    message("Saved: ", output_path)
  }

})

# Heatmaps retain the orange/red palette and numeric labels of the original HTML.
heatmap <- function(mat,path,title) {
 if(!is.matrix(mat) || !nrow(mat) || !ncol(mat) ||
    length(rownames(mat))!=nrow(mat) || length(colnames(mat))!=ncol(mat) ||
    anyNA(rownames(mat)) || anyNA(colnames(mat)))
  stop('Heatmap requires a nonempty matrix with row and column labels: ',title)
 d <- expand.grid(row=rownames(mat),column=colnames(mat)); d$value <- as.vector(mat)
 d$label <- ifelse(abs(d$value)>0 & abs(d$value)<.001,format(d$value,scientific=TRUE,digits=2),sprintf('%.3f',d$value))
 d$row <- factor(d$row,levels=rev(rownames(mat)))
 g <- ggplot2::ggplot(d,ggplot2::aes(column,row,fill=abs(value)))+
  ggplot2::geom_tile(colour='white')+ggplot2::geom_text(ggplot2::aes(label=label),size=2.7)+
  ggplot2::scale_fill_gradient(low='#FFF5EB',high='#CB181D',name=if(grepl('p-values',title,fixed=TRUE)) 'p-value' else 'Absolute difference')+
  ggplot2::theme_minimal(base_size=10)+ggplot2::labs(x=NULL,y=NULL,title=title)+
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1),panel.grid=ggplot2::element_blank())
 ggplot2::ggsave(path,g,width=10,height=max(4,nrow(mat)*.25+2),dpi=200)
 path
}
# All summaries below are descriptive summaries of existing exports.
edge_summary <- function(e) {
 out <- lapply(split(e,list(e$network,e$domain),drop=TRUE),function(z)
  data.frame(network=z$network[1],domain=z$domain[1],possible_edges=nrow(z),
   retained_edges=sum(z$weight!=0),positive=sum(z$weight>0),negative=sum(z$weight<0),
   maximum_absolute_weight=max(abs(z$weight)),sum_absolute_weights=sum(abs(z$weight))))
 do.call(rbind,out)
}
variance_table <- function(d) {
 z <- d[c('modality','group','variable','family_variance','scanner_variance','residual_variance')]
 denominator <- rowSums(z[c('family_variance','scanner_variance','residual_variance')],na.rm=TRUE)
 z$family_variance_fraction <- z$family_variance/denominator
 z$scanner_variance_fraction <- z$scanner_variance/denominator; z
}
compact <- function(d,columns) d[intersect(columns,names(d))]
convergence <- compact(models,c('visit','target','N_eligible','N_model','events','parameters','converged','status'))
logistic <- get_result('02b_bias_reduced_logistic.csv')
logistic <- compact(logistic,c('visit','target','N','events','term','OR','lower','upper','p',
 'p_BH_all_nonintercept','p_Holm_all_nonintercept','sparse_event_model'))

build_report <- function(kind) {
 is_primary <- kind=='S2'
 e <- if(is_primary) primary_edges else wave_edges
 ctab <- if(is_primary) central_primary else central_wave
 dg <- if(is_primary) diag_primary else get_result('12_wave_specific_model_diagnostics.csv')
 nt <- if(is_primary) get_result('05_original_NCT_omnibus.csv') else {
  z <- get_result('13_sensitivity_NCT_omnibus.csv'); z[z$variant=='wave_specific',]
 }
 ne <- if(is_primary) get_result('05_original_NCT_edges.csv') else {
  z <- get_result('13_sensitivity_NCT_edges.csv'); z[z$variant=='wave_specific',]
 }
 cs <- get_result(if(is_primary) '04_bootnet_centrality_CS.csv' else '13_wave_specific_centrality_CS.csv')
 intervals <- get_result(if(is_primary) '04_bootnet_edge_intervals.csv' else '13_wave_specific_edge_intervals.csv')
 mi <- get_result(if(is_primary) '07_original_PMM_vs_primary_DESCRIPTIVE.csv' else '14_wave_PMM_vs_wave_CC_DESCRIPTIVE.csv')
 mi_audit <- get_result(if(is_primary) '07_original_PMM_eligibility_audit.csv' else '14_wave_PMM_input_audit.csv')
 label <- if(is_primary) 'primary' else 'wave'
 title <- if(is_primary) 'Supplementary File S2: Primary pooled-residualisation analysis' else
  'Supplementary File S3: Wave-specific residualisation sensitivity analysis'
 parts <- character()
 add <- function(n,t,s) parts <<- c(parts,section(n,t,s))
 add(1,'Overview',p(if(is_primary)
  'This report presents the revised primary complete-case networks and accompanying analyses. Covariate adjustment is estimated using pooled observations; separate between-person networks are estimated at each visit. The analysis does not estimate within-person developmental change.' else
  'This report presents the sensitivity analysis in which covariate adjustment is estimated separately at each visit, using the same participants as the primary networks. The resulting networks are between-person networks. The primary pooled analysis is reported in Supplementary File S2.'))
 decisions <- data.frame(Component=c('Population','Visits','Imaging quality','Family sampling','Adjustment',
  'Family contribution','Network estimation','Uncertainty','Wave comparisons','Missing data'),
  Approach=c('Social media account holders at the corresponding visit',
  'ses-02A, ses-04A and ses-06A','Separate recommended inclusion flags for FA and T2-weighted measurements',
  'One participant per family, selected once and retained across visits',
  if(is_primary) 'Pooled adjustment without visit fixed effects or age-by-visit interaction; see model formulas below' else
   'Separate adjustment within each visit; see model formulas below',
  if(is_primary) 'Fixed effects and scanner random effects removed; family contribution retained' else
   'Within a visit, one participant per family does not support estimation of a family random intercept',
  'Pearson correlations; EBICglasso tuning parameter 0.5; six SMA and twelve imaging nodes',
  '1,000 nonparametric and case-dropping bootstrap resamples',
  'Paired Network Comparison Tests on common participants; 1,000 permutations',
  'Complete-case primary networks; PMM sensitivity restricted to eligible account holders and modality-specific QC-pass records'))
 selection_definitions <- data.frame(Measure=c('SMA account eligibility','FA inclusion','T2-weighted inclusion','FA motion covariate'),
 Variable=c('nt_y_stq__socmed_001','mr_y_qc__incl__dmri_indicator',
 'mr_y_qc__incl__smri__t2_indicator','mr_y_qc__raw__dmri__mot__ok_mean'),
 Rule=c('Account ownership = 1 at the relevant visit; 0 is structurally ineligible; unknown status is excluded',
 'Inclusion flag = 1 at the relevant visit; 0 or missing does not pass',
 'Inclusion flag = 1 at the relevant visit; 0 or missing does not pass',
 'Continuous motion summary in millimetres; used in FA adjustment'))
 add(2,'Analytic decisions',paste0(tab(decisions,'Analytic design and estimation'),
 tab(selection_definitions,'Eligibility, imaging quality and motion variables'),
 p('Visits are joined by participant identifier and decoded visit label, rather than by source-specific numerical visit codes. The age covariate in network adjustment is measured in years. MRI variable names and node mappings are listed in Section 9; tissue definitions and acquisition details are provided in Supplementary File S1.')))
 if(is_primary) {
 add(3,'Sample flow, participation and retention',paste0(
  tab(get_result('02_sequential_sample_flow.csv'),'Sequential sample flow by visit and imaging modality'),
  tab(get_result('18_disjoint_exclusion_group_counts.csv'),'Mutually exclusive inclusion and exclusion groups'),
  tab(get_result('19a_retention_wave_availability.csv'),'Visit participation, eligibility and analytic-data availability'),
  tab(get_result('19a_retention_patterns.csv'),'Participation and account-status patterns; visit order is ses-02A, ses-04A, ses-06A'),
  tab(get_result('19a_retention_outcome_counts.csv'),'Sustained analytic-data availability across all three visits'),
  p('Retention here means sustained eligible analytic-data availability. It combines visit participation, account status, imaging quality and completeness; it is not a measure of study dropout alone. The person-level population consists of the locked family sample with an account at one or more target visits, retaining the original union-of-visits design while excluding never-eligible participants.'),
  tab(get_result('19a_retained_vs_not_retained_descriptives.csv'),'Demographic comparisons by sustained analytic-data availability'),
  p('Selection comparisons describe observed differences and do not establish representativeness. Parental education is compared only among records with available information at the same participant and visit; missing values are not carried forward.'),
  tab(get_result('15_education_join_and_observed_counts.csv'),'Parental education availability'),
  tab(get_result('15_parental_education_selection_SMD.csv'),'Parental education comparisons'),
  tab(get_result('18_included_vs_each_exclusion_group_SMD.csv'),'Included participants versus individual exclusion groups: standardised differences')))
 add(4,'Retention, selection and nonresponse models',paste0(
  p('The three retention models preserve the original outcomes: any eligible SMA-imaging data, complete FA data and complete T2-weighted data at all three visits. Predictors are age, sex, site, race, ethnicity and household income. As in the original script, ses-02A values are used when observed, otherwise the first observed target-visit value is used. The source-visit counts below identify this fallback; such values are not necessarily baseline measurements.'),
  tab(get_result('19a_predictor_source_visits.csv'),'Source visits for retention-model predictors'),
  tab(get_result('19a_retention_model_summary.csv'),'Retention-model samples and estimation methods'),
  tab(get_result('19a_retention_reference_categories.csv'),'Reference categories for retention-model predictors'),
  tab(get_result('19a_retention_logistic_coefficients.csv'),'Retention-model odds ratios and approximate Wald intervals'),
  p('Maximum-likelihood logistic regression is retained when estimable. Mean bias reduction is used only for separation or failed numerical checks, with the method identified in the table. Age is centered within each model; the retention-model age odds ratio is per month. Raw p-values and BH/Holm corrections across the non-intercept coefficients of the three retention models are reported as one family.'),
  p('The following visit-specific selection and SMA nonresponse models use mean bias-reduced logistic regression and approximate Wald intervals. Their age odds ratios are per year. Sparse-event flags identify limited information relative to the number of parameters. These are distinct from the three-visit retention models. BH/Holm corrections cover all non-intercept coefficients of the nine visit-specific models as a separate family.'),
  tab(convergence,'Visit-specific model sample sizes and convergence'),tab(logistic,'Visit-specific selection and nonresponse odds ratios')))
 add(5,'Missingness quantity and SMA distributions',paste0(
  tab(get_result('02_structural_vs_item_missingness.csv'),'Account status and completeness of the six SMA items'),
  p('Items skipped because the participant reported no account are structurally inapplicable and are excluded from SMA imputation. Missing account status is not assumed to indicate eligibility.'),
  tab(get_result('18_SMA_distributions_with_kurtosis.csv'),'SMA item distributions among eligible participants'),
  tab(get_result('18_SMA_category_frequencies.csv'),'Counts in each of the six response categories, with observed denominators'),
  tab(get_result('02_account_and_visit_transitions.csv'),'Account-status transitions across visits')))
 add(6,'Missingness assessment',paste0(
  tab(get_result('12_Little_MCAR_eligible_AND.csv'),'Little\u2019s MCAR tests on pooled-adjusted values in eligible AND samples'),
  p('These tests use the primary pooled-adjusted values within the indicated eligible AND samples. Rejection of MCAR does not establish MAR. PMM does not rule out MNAR or correct all selection bias. Imputation is restricted to eligible, QC-pass records.'),
  tab(mi_audit,'Pooled-adjustment PMM samples and imputation information'),
  tab(get_result('02_motion_audit.csv'),'Motion availability and descriptive checks')))
 } else {
 add(3,'Sample flow, participation and retention',paste0(
  p('Wave-specific complete-case networks use the same participants as the corresponding primary networks. Eligibility, imaging QC, family selection, participation patterns, retention outcomes and included-versus-excluded comparisons are shared with the primary analysis and reported in Supplementary File S2, Section 3.'),
  tab(get_result('02_sequential_sample_flow.csv'),'Shared complete-case sample flow by visit and modality')))
 add(4,'Retention, selection and nonresponse models',p('The retention, selection and nonresponse models describe the shared participant selection process and are reported in Supplementary File S2, Section 4. They are not re-estimated as wave-specific residualisation models.'))
 add(5,'Missingness quantity and SMA distributions',p('Account eligibility, structural missingness, item distributions and response-category frequencies are properties of the shared source sample and are reported in Supplementary File S2, Section 5.'))
 add(6,'Missingness assessment and wave-specific PMM',paste0(
  p('The MCAR tests in Supplementary File S2 use pooled-adjusted values. They are not presented here as independently repeated wave-specific tests. The PMM comparisons below use wave-specific adjusted values and their corresponding wave-specific complete-case networks, keeping residualisation consistent within each comparison.'),
  tab(mi_audit,'Wave-specific PMM samples and imputation information')))
 }
 add(7,'Clustering and variance components',paste0(
  p('Variance fractions below are calculated from the fitted model variance components. They describe the adjustment models and should not be interpreted as variance explained by the network.'),
  tab(variance_table(dg),'Model variance components and clustering fractions')))
 add(8,'Residualisation diagnostics',paste0(
  p(if(is_primary) 'The network input is the observed value minus the fixed-effect prediction and, for imaging variables, the scanner-specific random-intercept prediction. Family random-intercept predictions are retained. FA adjustment includes the selected dMRI motion measure. The pooled model contains no visit fixed effect.' else
  'Adjustment is estimated separately at each visit. Fixed-effect predictions and scanner random-intercept predictions are removed. With one child per family at a given visit, the family intercept cannot be identified separately from the observation-level residual.'),
  tab(dg,'Adjustment formulas and residualisation diagnostics'),
  if(is_primary) tab(get_result('12_primary_residual_reconstruction.csv'),'Reconstruction of the primary network inputs') else ''))
 add(9,if(is_primary) 'Primary network results' else 'Wave-specific sensitivity network results',paste0(
  tab(get_result('03_node_dictionary.csv'),'Node definitions'),tab(edge_summary(e),'Edges by domain and network'),
  paste(vapply(keys,function(k) image_tag(file.path(FIG,paste0(label,'_',k,'_network.png')),
   paste(gsub('_',' ',k),'network; blue edges are positive and red edges are negative.')),character(1)),collapse=''),
  tab(e,'All estimated edges; zero weights are retained in the table'),
  tab(ctab[c('network','node','strength','expected_influence')],'Strength and expected influence')))
 add(10,'Bridge centrality',paste0(
  p('SMA and MRI form the two predefined domains. Bridge strength is the sum of absolute cross-domain edge weights; bridge expected influence is their signed sum. These are descriptive network summaries, not clinical markers.'),
  tab(ctab[c('network','node','bridge_strength','bridge_expected_influence')],'Bridge centrality')))
 add(11,'Network comparisons across visits',paste0(
  p('Each comparison uses participants observed at both visits. Joint global strength combines SMA\u2013SMA, MRI\u2013MRI and SMA\u2013MRI edges; a difference does not by itself demonstrate altered symptom\u2013MRI coupling. Edge tests report BH and Holm corrections within each family of 153 edges.'),
  tab(nt,'Paired network structure and global strength comparisons'),
  tab(ne,'Individual-edge comparisons and multiplicity-adjusted p-values'),
  if(is_primary) paste0(p('Domain-specific permutation tests retain the full 18-node network and compare summed absolute edge weights within each domain. BH and Holm corrections are applied across the 18 domain comparisons.'),
   tab(get_result('16_domain_strength_NCT_18_tests.csv'),'Domain-specific network strength comparisons')) else ''))
 heatmaps <- character()
 for(metric in c('strength','expected_influence')) {
  columns <- unlist(lapply(c('FA','T2w'),function(mod) c(paste0(mod,': 02A\u219204A'),paste0(mod,': 04A\u219206A'),paste0(mod,': 02A\u219206A'))))
  mat <- matrix(NA_real_,18,6,dimnames=list(nodes,columns)); col <- 0L
  for(mod in c('FA','T2w')) for(pair in list(c('02A','04A'),c('04A','06A'),c('02A','06A'))) {
   col <- col+1L
   a <- ctab[ctab$network==paste0(mod,'_ses-',pair[1]),]; b <- ctab[ctab$network==paste0(mod,'_ses-',pair[2]),]
   mat[,col] <- b[[metric]][match(nodes,b$node)]-a[[metric]][match(nodes,a$node)]
  }
  path <- heatmap(mat,file.path(FIG,paste0(label,'_',metric,'_differences.png')),
   paste(gsub('_',' ',metric),'differences (later minus earlier visit)'))
  heatmaps <- c(heatmaps,image_tag(path,paste('Descriptive',gsub('_',' ',metric),'differences; values are signed.')))
 }
 comparisons <- unique(ne$comparison)
 identity <- paste(ne$node1,ne$node2,sep=' / ')
 top <- tapply(abs(ne$weight_B-ne$weight_A),identity,max)
 selected <- names(sort(top,decreasing=TRUE))[seq_len(min(20L,length(top)))]
 mat <- matrix(NA_real_,length(selected),length(comparisons),dimnames=list(selected,comparisons))
 for(j in seq_along(comparisons)) {
  z <- ne[ne$comparison==comparisons[j],]; id <- paste(z$node1,z$node2,sep=' / ')
  mat[,j] <- (z$weight_B-z$weight_A)[match(selected,id)]
 }
 path <- heatmap(mat,file.path(FIG,paste0(label,'_edge_differences.png')),'Largest absolute edge differences (later minus earlier visit)')
 heatmaps <- c(heatmaps,image_tag(path,'Twenty edges selected descriptively by the maximum absolute difference across comparisons. This display does not select statistically significant edges.'))
 if(!'comparison' %in% names(nt))
  nt$comparison <- paste0(nt$modality,'_',nt$visit_A,'_vs_',nt$visit_B)
 if(length(nt$comparison)!=nrow(nt) || anyNA(nt$comparison) || anyDuplicated(nt$comparison))
  stop('Invalid NCT comparison labels.')
 mat <- cbind('Structure p'=nt$p_structure,'Joint strength p'=nt$p_global)
 rownames(mat) <- nt$comparison
 path <- heatmap(mat,file.path(FIG,paste0(label,'_NCT_p_values.png')),'Raw NCT p-values')
 add(12,'Heatmaps',paste0(p('Orange/red intensity represents the absolute displayed value. Signed differences are printed in cells. These figures summarise the corresponding tables and add no statistical tests.'),
  paste(heatmaps,collapse=''),image_tag(path,'Unadjusted omnibus NCT p-values; see the full table for sample sizes.')))
 add(13,'Sensitivity analyses and cross-domain uncertainty',paste0(
  if(is_primary) paste0(
   p('Ordinal sensitivity uses ordered raw SMA items, adjusted MRI measures and age in a mixed correlation matrix: polychoric correlations for ordered-item pairs, polyserial correlations for ordered/continuous pairs, and Pearson correlations for continuous pairs. Age is partialled out before EBICglasso estimation. No automatic positive-definite matrix repair is applied. The matched Pearson comparison uses the same participants and age adjustment; comparison with the primary networks additionally changes SMA adjustment.'),
   tab(get_result('06_sensitivity_network_comparisons.csv'),'Sensitivity comparisons against primary pooled networks or the labelled matched Pearson reference'),
   tab(get_result('06_ordinal_matrix_diagnostics.csv'),'Ordinal correlation-matrix diagnostics')) else
   p('The ordinal, alternative-motion and family-intercept sensitivity comparisons use the primary pooled pipeline as their reference and are reported in Supplementary File S2, Section 13. They are not independent wave-specific sensitivity estimates.'),
  p('PMM comparisons are descriptive and report every imputation. AND and OR rules can yield different participant populations; their differences do not isolate imputation effects.'),
  tab(mi,'PMM versus the corresponding complete-case networks'),
  if(is_primary) paste0(tab(get_result('17_family_selection_audit.csv'),'Five additional family selections'),
   tab(get_result('17_family_selection_network_comparisons.csv'),'Family-selection network comparisons'),
   tab(get_result('17_family_selection_NCT_omnibus.csv'),'Paired NCT results across family selections'),
   tab(subset(get_result('13_sensitivity_NCT_omnibus.csv'),variant=='no_family_pooled'),
    'Paired NCT sensitivity without the pooled family random intercept'),
   tab(subset(get_result('13_sensitivity_NCT_edges.csv'),variant=='no_family_pooled'),
    'Edge-comparison sensitivity without the pooled family random intercept'),
   tab(get_result('12_ORIGINAL_posthoc_regression_sensitivity.csv'),'Selected-edge post hoc regression sensitivity'),
   p('Unregularised cross-domain partial correlations condition on the other 16 nodes. Approximate Fisher-z intervals are conditional on the saved preprocessing. These estimates and their tests differ from EBICglasso edge estimates.'),
   tab(get_result('13_unregularised_432_partial_correlations.csv'),'Unregularised cross-domain estimates and approximate intervals'),
   p('Exploratory equivalence analyses use bounds of plus or minus 0.05 and 0.03, chosen during revision. These bounds are not validated clinical thresholds. Holm correction is applied across 432 pairs separately at each bound; raw and within-network results are shown for transparency. Failure to establish equivalence does not establish a meaningful association.'),
   tab(get_result('13_equivalence_summary.csv'),'Exploratory equivalence summaries'),
   tab(get_result('13_equivalence_TOST_two_bounds.csv'),'Exploratory equivalence results for every pair')) else
   p('Family-selection, domain-specific strength and unregularised equivalence analyses of the primary pooled inputs are reported in Supplementary File S2.')))
 bootstrap_images <- paste(vapply(keys,function(k) paste0(
  image_tag(file.path(FIG,paste0(label,'_',k,'_case_stability.png')),paste(k,'case-dropping stability')),
  image_tag(file.path(FIG,paste0(label,'_',k,'_edge_accuracy.png')),paste(k,'bootstrap edge accuracy'))),character(1)),collapse='')
 add(14,'Bootstrap stability and accuracy',paste0(
  p('Nonparametric and case-dropping bootstrap procedures use 1,000 resamples. Edge intervals are pointwise and conditional on preprocessing. Regularised edge intervals are not used as significance tests. Centrality stability is reported for strength and expected influence.'),
  tab(cs,'Centrality correlation-stability coefficients'),tab(intervals,'Pointwise bootstrap edge intervals'),
  if(is_primary) paste0(tab(get_result('17_family_selection_centrality_CS.csv'),'Centrality stability across the five family selections'),
   tab(get_result('17_family_selection_edge_intervals.csv'),'Pointwise edge intervals across the five family selections')) else '',bootstrap_images))
 add(15,'Interpretation',p('The networks describe conditional associations among account-holding adolescents who met the indicated imaging and analytic criteria. Separate visit networks and paired comparisons do not estimate within-person change or causal effects. Small regularised cross-domain edges, selection differences and uncertainty limit clinical interpretation. Neither a zero EBICglasso weight nor an uncorrected null test proves practical absence of an association.'))
 definitions <- data.frame(Term=c('SMA\u2013SMA','MRI\u2013MRI','SMA\u2013MRI','Joint global strength','BH','Holm','AND','OR','PMM'),
  Definition=c('Edges connecting two SMA symptoms','Edges connecting two imaging measures',
  'Edges connecting a SMA symptom with an imaging measure','Sum of absolute weights of all 153 possible edges',
  'Benjamini\u2013Hochberg false discovery rate correction','Holm familywise error correction',
  'Observed information in both domains, within the eligible and QC-pass population',
  'Observed information in at least one domain, within the eligible and QC-pass population',
  'Predictive mean matching; sensitivity imputation on the adjusted scale'))
 add(16,'Reporting definitions',tab(definitions,'Definitions used in this report'))
 doc <- templates[[kind]]
 container <- xml2::xml_find_first(doc,'//div[contains(@class,"toc-content")]')
 xml2::xml_remove(xml2::xml_contents(container))
 content <- paste0('<div id="header"><h1 class="title toc-ignore">',esc(title),'</h1></div>',paste(parts,collapse='\n'))
 fragment <- xml2::read_html(paste0('<html><body>',content,'</body></html>'))
 for(node in xml2::xml_children(xml2::xml_find_first(fragment,'//body'))) xml2::xml_add_child(container,node)
 xml2::xml_set_text(xml2::xml_find_first(doc,'//title'),title)
 xml2::xml_remove(xml2::xml_find_all(doc,'//meta[@name="date"]'))
 # Retain original styles and scripts, but remove the obsolete code menu.
 text <- xml2::xml_text(container)
 if(grepl('pending|provisional|TODO|AI-generated|Using results folder',text,ignore.case=TRUE))
  stop('An editorial marker was found in the report.')
 if(length(xml2::xml_find_all(container,'.//div[contains(@class,"section level1")]'))!=16L)
  stop('The report must contain all 16 sections.')
 path <- file.path(OUT,paste0('Supplementary_File_',kind,'.building.html'))
 xml2::write_html(doc,path,options='format')
 path
}
reports <- c(build_report('S2'),build_report('S3'))
final_reports <- file.path(OUT,c('Supplementary_File_S2.html','Supplementary_File_S3.html'))
if(!all(file.copy(reports,final_reports,overwrite=TRUE))) stop('Could not save both final HTML reports.')
unlink(reports);reports <- final_reports
figure_index <- do.call(rbind,lapply(keys,function(k) data.frame(
 Figure=c('Supplementary Figure 1','Supplementary Figure 2','Supplementary Figure 3',
  if(k=='FA_ses-02A') 'Main manuscript network' else
   paste('Supplementary Figure',match(k,keys)+3L)),
 Panel=k,
 File=paste0('primary_',k,c('_centrality.png','_case_stability.png','_edge_accuracy.png','_network.png')))))
figure_index <- rbind(figure_index,data.frame(Figure='Supplementary Figure 4',
 Panel='Analytical workflow',File='Figure_4_Analytical_Workflow.png'))
write.csv(figure_index,file.path(OUT,'S1_figure_index.csv'),row.names=FALSE)
# Record sources privately; no participant IDs or local paths appear in the HTML.
write.csv(data.frame(file=required,md5=unname(tools::md5sum(vapply(required,find_file,character(1))))),
 file.path(OUT,'report_source_checksums.csv'),row.names=FALSE)
message('Reports and S1 statistical figures created: ',normalizePath(OUT,winslash='/'))
message('Send the two HTML reports for a content and appearance check before journal submission.')
