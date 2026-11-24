paper1_template <- function(cfg, kind, directory) {
  p <- cfg$templates[[kind]]
  if (!is.null(p) && nzchar(p)) return(paper1_path(p))
  p <- file.path(directory, paste0("template_", kind, ".html"))
  writeLines(c("<!doctype html><html><head><meta charset='utf-8'>",
    paste0("<title>Supplementary File ", kind, "</title>"),
    "<style>body{font-family:Arial,sans-serif;max-width:1200px;margin:36px auto;padding:0 20px;color:#222}table{border-collapse:collapse;display:block;overflow:auto;max-width:100%;margin:20px 0}th,td{border:1px solid #ddd;padding:6px;text-align:left;font-size:12px}th{background:#f1f3f5}caption{text-align:left;font-weight:bold;margin:8px 0}img{max-width:100%;height:auto}h1,h2,h3{line-height:1.3}p{line-height:1.6}</style></head><body>",
    paste0("<div class='toc-content'><h1 class='title'>Supplementary File ", kind, "</h1></div>"),
    "</body></html>"), p)
  p
}
paper1_reports <- function(cfg, revision = TRUE) {
  paper1_missing_packages(c("xml2", "base64enc", "ggplot2", "qgraph", "pdftools"))
  baseline <- paper1_baseline(cfg)
  target <- file.path(cfg$work_dir, if (revision) "R2_Reports" else "R1_Reports")
  snapshot <- file.path(target, "inputs"); dir.create(snapshot, recursive = TRUE, showWarnings = FALSE)
  engine <- readLines(file.path(.paper1_repo, "R", "engines", "19b_Paper1_Journal_Reports.R"), warn = FALSE)
  # Read the declared inputs from the report engine without executing the report.
  a <- grep("^required <- c\\(", engine); b <- grep("^D <- setNames", engine)
  if (length(a) != 1L || length(b) != 1L || b <= a) stop("Report input declaration changed.")
  ev <- new.env(parent = baseenv()); eval(parse(text = engine[a:(b - 1L)]), ev)
  required <- ev$required
  r2 <- file.path(cfg$work_dir, "R2_NCT_10000")
  if (revision) {
    marker <- file.path(r2, "RESTRICTED", "COMPLETE.rds")
    if (!file.exists(marker)) stop("Validated R2 comparisons are required.")
    completed <- readRDS(marker)
    if (!identical(unname(tools::md5sum(names(completed$checksums))), unname(completed$checksums)))
      stop("R2 results changed after validation.")
  }
  # Copy only named aggregate inputs/PDFs into a private report workspace.
  # Nothing is copied into the repository, and these tables still need disclosure review.
  for (name in required) {
    source <- if (revision && name == "05_original_NCT_omnibus.csv") {
      file.path(r2, "21_NCT10000_omnibus.csv")
    } else if (revision && name == "05_original_NCT_edges.csv") {
      file.path(r2, "21_NCT10000_edges.csv")
    } else {
      actual <- if (name == "13_equivalence_TOST_two_bounds.csv") "13_equivalence_TOST_both_VAS_bounds.csv" else name
      paper1_find(baseline$parent, actual)
    }
    if (revision && name == "05_original_NCT_omnibus.csv") {
      # Prior p-values remain in the private author audit, not the updated result table.
      z <- read.csv(source, check.names = FALSE)
      cols <- c("comparison", "modality", "visit_A", "visit_B", "paired_N", "permutations",
                "M", "p_structure", "global_A", "global_B", "p_global")
      if (!all(cols %in% names(z))) stop("Unexpected revised omnibus columns.")
      write.csv(z[cols], file.path(snapshot, name), row.names = FALSE)
    } else if (!file.copy(source, file.path(snapshot, name), overwrite = TRUE))
      stop("Could not stage report input: ", name)
  }
  keys <- unique(read.csv(file.path(snapshot, "03_all_network_edges.csv"))$network)
  for (k in keys) for (name in c(paste0("04_", k, "_edge_accuracy.pdf"),
      paste0("04_", k, "_case_stability.pdf"), paste0("13_", k, "_wave_edge_accuracy.pdf"),
      paste0("13_", k, "_wave_case_stability.pdf"))) {
    if (!file.copy(paper1_find(baseline$parent, name), file.path(snapshot, name), overwrite = TRUE))
      stop("Could not stage bootstrap figure.")
  }
  retention <- dirname(paper1_find(baseline$parent, "19a_retention_model_summary.csv"))
  s2 <- paper1_template(cfg, "S2", target); s3 <- paper1_template(cfg, "S3", target)
  if (revision) {
    line <- "if(!all(get_result('05_original_NCT_omnibus.csv')$permutations==1000L)) stop('Unexpected NCT permutation count.')"
    if (sum(engine == line) != 1L) stop("Report permutation-count check changed.")
    engine[engine == line] <- sub("==1000L", "==10000L", line, fixed = TRUE)
    line <- "  'Paired Network Comparison Tests on common participants; 1,000 permutations',"
    if (sum(engine == line) != 1L) stop("Report comparison description changed.")
    engine[engine == line] <- "  if(is_primary) 'Paired Network Comparison Tests on common participants; 10,000 permutations' else 'Paired Network Comparison Tests on common participants; 1,000 permutations',"
  }
  e <- paper1_scope(snapshot, c(file.path(snapshot, "03_all_network_edges.csv"), s2, s3,
                                file.path(retention, "19a_retention_model_summary.csv")))
  eval(parse(text = engine), e)
  generated <- file.path(snapshot, "Journal_Supplementary_Reports")
  # Apply the final October equivalence reporting (uncorrected pairwise + all-72 IUT).
  eq_sum <- paper1_find(baseline$parent, "VAS_equivalence_counts_percent_IUT_by_network.csv")
  eq_largest <- paper1_find(baseline$parent, "VAS_largest_partial_r_90CI_by_network.csv")
  eq <- paper1_scope(snapshot, c(file.path(generated, "Supplementary_File_S2.html"), eq_sum, eq_largest))
  sys.source(file.path(.paper1_repo, "R", "engines", "20_Paper1_Update_S2_and_Table33_VAS_Exact_v2.R"), eq)
  updated <- file.path(generated, "Supplementary_File_S2_VAS_IUT_UPDATED.html")
  if (!file.exists(updated)) stop("Final equivalence update did not produce S2.")
  docs <- list(S2 = xml2::read_html(updated), S3 = xml2::read_html(file.path(generated, "Supplementary_File_S3.html")))
  for (kind in names(docs)) {
    doc <- docs[[kind]]
    container <- xml2::xml_find_first(doc, '//div[contains(@class,"toc-content")]')
    note <- paste0("Omnibus tests address whether the joint network differs overall; domain-specific tests address differences in summed absolute weights within predefined edge domains. These distinct questions are reported in separate testing families. Omnibus p-values are unadjusted, remain exploratory, and do not control multiplicity across all omnibus claims. Results should be interpreted alongside the sensitivity analyses. ",
      if (kind == "S2") { if (revision) "The six primary paired comparisons use 10,000 permutations. " else "The six primary paired comparisons use 1,000 permutations. " } else "Wave-specific paired comparisons retain 1,000 permutations. ",
      "Family-selection and other R1 sensitivity comparisons retain 1,000 permutations. Their edge-level non-rejections must not be interpreted as evidence of absence given the available p-value resolution.")
    xml2::xml_add_child(container, "p", note)
    if (kind == "S2") xml2::xml_add_child(container, "p",
      "With one participant selected per family, the pooled family random intercept effectively represents participant-specific stable variation across assessments. Its prediction was not subtracted from the network inputs.")
    # Remove raw-code widgets inherited from user-supplied templates, not narrative tables.
    xml2::xml_remove(xml2::xml_find_all(doc, '//pre[contains(@class,"sourceCode")]'))
    body <- xml2::xml_text(xml2::xml_find_first(doc, "//body"))
    if (grepl("(?:sub-|NDAR_INV)[A-Z0-9]{6,}", body, perl = TRUE))
      stop("A participant-identifier pattern was detected in a report. Do not share it.")
    doc <- paper1_report_wording(doc, kind, cfg, revision)
    xml2::write_html(doc, file.path(target, paste0("Supplementary_File_", kind, ".html")), options = "format")
  }
  if (revision) {
    paper1_reviewer_summaries(cfg)
    summary_dir <- file.path(cfg$work_dir, "Reviewer_Revision_Summaries")
    for (name in c("22_family_selection_test_consistency.csv", "22_family_selection_network_similarity.csv")) {
      x <- read.csv(file.path(summary_dir, name), check.names = FALSE)
      # Add an explicit numbered-review follow-up table to S2 without inventing a S1 table number.
      doc <- xml2::read_html(file.path(target, "Supplementary_File_S2.html"))
      parent <- xml2::xml_find_first(doc, '//div[contains(@class,"toc-content")]')
      xml2::xml_add_child(parent, "h2", if (grepl("test_consistency", name)) "Family-selection inference sensitivity" else "Family-selection network similarity")
      if (grepl("test_consistency", name)) {
        focus <- x[x$modality == "FA" &
          ((x$visit_A == "ses-02A" & x$visit_B == "ses-04A" & x$test == "p_global") |
           (x$visit_A == "ses-02A" & x$visit_B == "ses-06A")), , drop = FALSE]
        statements <- vapply(seq_len(nrow(focus)), function(i) {
          z <- focus[i, , drop = FALSE]
          test_name <- if (z$test == "p_global") "global strength" else "network structure"
          sprintf("For FA %s versus %s, the unadjusted %s p-value was below .05 in %d of %d additional selections.",
            z$visit_A, z$visit_B, test_name, z$raw_p_below_05, z$repetitions)
        }, character(1))
        xml2::xml_add_child(parent, "p", paste(c(statements,
          "High similarity of edge weights does not imply identical test conclusions. These overlapping samples are sensitivity analyses, not independent replications."), collapse = " "))
      }
      tab <- xml2::xml_add_child(parent, "table")
      head <- xml2::xml_add_child(tab, "tr")
      for (nm in names(x)) xml2::xml_add_child(head, "th", nm)
      for (i in seq_len(nrow(x))) {
        row <- xml2::xml_add_child(tab, "tr")
        for (v in x[i, , drop = FALSE]) xml2::xml_add_child(row, "td", as.character(v))
      }
      xml2::write_html(doc, file.path(target, "Supplementary_File_S2.html"), options = "format")
    }
  }
  # Preserve S1 figures/CSV replacements locally. Word documents are not silently modified.
  message("Updated reports: ", target)
  message("S1 figures and replacement Table 33 are in inputs/Journal_Supplementary_Reports.")
  message("Update manuscript/S1 numerical references from these outputs and review small cells before journal/repository release.")
  invisible(target)
}

# Reporting-only changes; no statistical estimates or input files are modified.
paper1_report_wording <- function(doc, kind, cfg, revision) {
  if (kind == "S3") {
    cells <- xml2::xml_find_all(doc, "//td")
    labels <- trimws(xml2::xml_text(cells))
    replace <- cells[labels == "Family contribution retained"]
    if (length(replace)) xml2::xml_text(replace) <-
      "Fixed effects and, where fitted, scanner random-intercept predictions removed; no family random intercept fitted"
  }
  if (kind == "S2") {
    paragraphs <- xml2::xml_find_all(doc, "//p")
    previous <- grepl("Holm-across-432 columns", xml2::xml_text(paragraphs), fixed = TRUE)
    if (any(previous)) xml2::xml_text(paragraphs[previous]) <-
      "The detailed pairwise table also reports Holm-adjusted p-values. These address multiplicity-adjusted pairwise claims. The network-level all-72 equivalence claim uses the intersection-union criterion described above and is supported only when every pair passes its unadjusted TOST at alpha = .05."
    if (isTRUE(revision)) {
      e <- read.csv(file.path(cfg$work_dir, "R2_NCT_10000", "21_NCT10000_edges.csv"), check.names = FALSE)
      required <- c("comparison", "node1", "node2", "p_BH_within153", "p_Holm_within153")
      if (!all(required %in% names(e)) || nrow(e) != 918L) stop("Unexpected revised edge table for report wording.")
      bh <- e$p_BH_within153 < .05; holm <- e$p_Holm_within153 < .05
      if (anyNA(bh) || anyNA(holm)) stop("Missing corrected edge p-values.")
      sig <- e[bh | holm, , drop = FALSE]
      if (nrow(sig)) {
        descriptions <- vapply(seq_len(nrow(sig)), function(i) {
          z <- sig[i, , drop = FALSE]
          correction <- if (z$p_BH_within153 < .05 && z$p_Holm_within153 < .05) "both BH and Holm" else if (z$p_BH_within153 < .05) "BH" else "Holm"
          sprintf("%s: %s to %s (%s; BH-adjusted p = %.6g; Holm-adjusted p = %.6g)",
            z$comparison, z$node1, z$node2, correction, z$p_BH_within153, z$p_Holm_within153)
        }, character(1))
        summary <- paste0("Individual-edge differences surviving at least one correction across 153 edges within each comparison: ",
          paste(descriptions, collapse = "; "),
          ". These tests concern differences between estimates from paired-participant networks. They do not establish that either individual edge differs from zero, clinical importance, or within-person developmental change.")
      } else {
        summary <- "No individual-edge differences survived BH or Holm correction across 153 edges within each comparison. Non-rejection does not establish absence of a difference."
      }
      heading <- xml2::xml_find_first(doc, "//h1[contains(.,'Network comparisons across visits')]")
      if (inherits(heading, "xml_missing")) stop("Network-comparison report section not found.")
      xml2::xml_add_sibling(heading, "p", summary, .where = "after")
    }
  }
  doc
}
