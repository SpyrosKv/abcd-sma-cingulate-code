# 20_Paper1_Update_S2_and_Table33_VAS_Exact.R
# Updates ONLY the equivalence-reporting content requested by VAS:
#   - Supplementary Table 33 replacement
#   - Supplementary File S2 equivalence subsection
#
# It does NOT refit models, networks, TOSTs, NCTs, imputations, or bootstraps.
# It does NOT overwrite the original S2 HTML.

options(stringsAsFactors = FALSE)

message("Choose the CURRENT Supplementary File S2 HTML (latest pooled-primary S2).")
p_s2 <- normalizePath(file.choose(), winslash="/", mustWork=TRUE)

message("Choose: VAS_equivalence_counts_percent_IUT_by_network.csv")
p_sum <- normalizePath(file.choose(), winslash="/", mustWork=TRUE)
if (basename(p_sum) != "VAS_equivalence_counts_percent_IUT_by_network.csv") {
  stop("Wrong equivalence-summary file selected.")
}

message("Choose: VAS_largest_partial_r_90CI_by_network.csv")
p_largest <- normalizePath(file.choose(), winslash="/", mustWork=TRUE)
if (basename(p_largest) != "VAS_largest_partial_r_90CI_by_network.csv") {
  stop("Wrong largest-partial-r file selected.")
}

s <- read.csv(p_sum, check.names=FALSE, na.strings=c("NA","NaN",""))
l <- read.csv(p_largest, check.names=FALSE, na.strings=c("NA","NaN",""))

need_s <- c("network","bound_r","pairs","equivalent_n","equivalent_percent","all_72_pass_IUT")
need_l <- c("network","largest_abs_partial_r","signed_partial_r","lower90","upper90")

miss_s <- setdiff(need_s, names(s))
miss_l <- setdiff(need_l, names(l))
if (length(miss_s)) stop("Missing summary columns: ", paste(miss_s, collapse=", "))
if (length(miss_l)) stop("Missing largest-r columns: ", paste(miss_l, collapse=", "))

if (nrow(s) != 12L) stop("Expected 12 network-by-bound summary rows.")
if (length(unique(s$network)) != 6L) stop("Expected six networks in equivalence summary.")
if (!setequal(unique(s$bound_r), c(.05,.03))) stop("Expected bounds 0.05 and 0.03.")
if (any(s$pairs != 72L)) stop("Expected 72 pairs per network/bound.")
if (nrow(l) != 6L || length(unique(l$network)) != 6L) {
  stop("Expected six rows in largest-partial-r summary.")
}
if (!setequal(unique(s$network), unique(l$network))) {
  stop("Network labels do not match between the two VAS outputs.")
}

# Preserve the order produced by the VAS analysis.
net_order <- unique(s$network)
s$network <- factor(s$network, levels=net_order)
s <- s[order(s$bound_r, s$network, decreasing=FALSE), , drop=FALSE]
# Put 0.05 first, then 0.03, matching the VAS request/email order.
s <- rbind(
  s[s$bound_r == .05, , drop=FALSE],
  s[s$bound_r == .03, , drop=FALSE]
)
s$network <- as.character(s$network)

# Merge the requested largest |partial r| and 90% CI into Table 33.
m <- merge(s, l, by="network", all.x=TRUE, sort=FALSE)
m$network <- factor(m$network, levels=net_order)
m <- rbind(
  m[m$bound_r == .05, , drop=FALSE][order(m[m$bound_r == .05, "network"]), , drop=FALSE],
  m[m$bound_r == .03, , drop=FALSE][order(m[m$bound_r == .03, "network"]), , drop=FALSE]
)
m$network <- as.character(m$network)

# Human-readable labels only; no analysis is changed.
pretty_network <- function(z) {
  z <- sub("^FA_", "FA ", z)
  z <- sub("^T2w_", "T2w ", z)
  z
}

table33 <- data.frame(
  Network = pretty_network(m$network),
  Bound = sprintf("\u00b1%.02f", m$bound_r),
  Pairs = m$pairs,
  `Equivalent n` = m$equivalent_n,
  `Equivalent %` = sprintf("%.1f", m$equivalent_percent),
  `All 72 pass IUT` = ifelse(m$all_72_pass_IUT, "Yes", "No"),
  `Largest |partial r|` = sprintf("%.4f", m$largest_abs_partial_r),
  `Signed partial r` = sprintf("%.4f", m$signed_partial_r),
  `90% CI` = sprintf("[%.4f, %.4f]", m$lower90, m$upper90),
  check.names=FALSE,
  stringsAsFactors=FALSE
)

# Table requested separately by VAS: largest |partial r| and 90% CI per network.
largest_table <- data.frame(
  Network = pretty_network(l$network),
  `Largest |partial r|` = sprintf("%.4f", l$largest_abs_partial_r),
  `Signed partial r` = sprintf("%.4f", l$signed_partial_r),
  `90% CI` = sprintf("[%.4f, %.4f]", l$lower90, l$upper90),
  check.names=FALSE,
  stringsAsFactors=FALSE
)

outdir <- dirname(p_s2)

write.csv(
  table33,
  file.path(outdir, "Supplementary_Table_33_VAS_IUT_UPDATED.csv"),
  row.names=FALSE,
  na=""
)

write.csv(
  largest_table,
  file.path(outdir, "Supplementary_Table_33_Largest_r_90CI_VAS_UPDATED.csv"),
  row.names=FALSE,
  na=""
)

table33_note <- paste(
  "Note. Unregularised partial correlations condition on the other 16 nodes;",
  "approximate Fisher-z SE = 1/sqrt(N-19). Pairwise TOSTs used alpha = .05",
  "without multiplicity correction, with 90% confidence intervals and bounds of",
  "+/- .05 and +/- .03. Counts and percentages are reported per network and bound.",
  "For the network-level claim that all 72 SMA-MRI associations fall within a",
  "given bound, the intersection-union principle was used: the claim is supported",
  "only if all 72 pairwise tests pass at alpha = .05. The bounds were selected",
  "during revision and are not clinically validated thresholds."
)

writeLines(
  table33_note,
  file.path(outdir, "Supplementary_Table_33_VAS_IUT_UPDATED_NOTE.txt")
)

# ---------- Update only the equivalence subsection in S2 HTML ----------

html <- paste(readLines(p_s2, warn=FALSE, encoding="UTF-8"), collapse="\n")

# Safety check: ensure this is the pooled-primary S2, not S3.
if (!grepl("Supplementary File S2: Primary pooled-residualisation analysis",
           html, fixed=TRUE)) {
  stop("The selected HTML does not appear to be Supplementary File S2.")
}

# Identify the existing equivalence section by its table captions.
# This is more robust than relying on one exact paragraph string.
summary_caption <- "<caption>Exploratory equivalence summaries</caption>"
detail_caption  <- "<caption>Exploratory equivalence results for every pair</caption>"

summary_cap_pos <- regexpr(summary_caption, html, fixed=TRUE)[1]
detail_cap_pos  <- regexpr(detail_caption, html, fixed=TRUE)[1]

if (summary_cap_pos < 1L) {
  stop("This S2 does not contain the existing 'Exploratory equivalence summaries' table. ",
       "Please select the latest S2 that contains the current Holm-based equivalence section.")
}
if (detail_cap_pos < 1L) {
  stop("This S2 does not contain the detailed 'Exploratory equivalence results for every pair' table.")
}
if (detail_cap_pos <= summary_cap_pos) stop("Unexpected S2 equivalence-section order.")

# Find the beginning of the paragraph immediately preceding the old summary table.
before_summary <- substr(html, 1L, summary_cap_pos - 1L)
p_starts <- gregexpr("<p>", before_summary, fixed=TRUE)[[1]]
p_starts <- p_starts[p_starts > 0L]
if (!length(p_starts)) stop("Could not locate the paragraph before the equivalence summary.")
start_pos <- max(p_starts)

old_intro <- substr(html, start_pos, summary_cap_pos - 1L)
if (!grepl("equivalence", old_intro, ignore.case=TRUE)) {
  stop("The paragraph before the equivalence summary does not appear to be the equivalence paragraph.")
}

# Find the start of the detailed table container, so the full pairwise audit remains untouched.
before_detail <- substr(html, 1L, detail_cap_pos - 1L)
table_starts <- gregexpr('<div style="overflow-x:auto"><table', before_detail, fixed=TRUE)[[1]]
table_starts <- table_starts[table_starts > summary_cap_pos]
if (!length(table_starts)) {
  # The detailed table's <div> starts before its caption; use the last table container before the detail caption.
  table_starts <- gregexpr('<div style="overflow-x:auto"><table', before_detail, fixed=TRUE)[[1]]
  table_starts <- table_starts[table_starts > 0L]
}
if (!length(table_starts)) stop("Could not locate the detailed pairwise table container.")
detail_start_pos <- max(table_starts)

html_escape <- function(z) {
  z <- gsub("&", "&amp;", z, fixed=TRUE)
  z <- gsub("<", "&lt;", z, fixed=TRUE)
  z <- gsub(">", "&gt;", z, fixed=TRUE)
  z
}

make_html_table <- function(df, caption) {
  hdr <- paste0(
    "<thead><tr>",
    paste0("<th>", html_escape(names(df)), "</th>", collapse=""),
    "</tr></thead>"
  )

  body_rows <- apply(df, 1, function(r) {
    paste0(
      "<tr>",
      paste0("<td>", html_escape(as.character(r)), "</td>", collapse=""),
      "</tr>"
    )
  })

  paste0(
    '<div style="overflow-x:auto"><table class="table table-condensed">\n',
    "<caption>", html_escape(caption), "</caption>\n",
    hdr, "\n<tbody>\n",
    paste(body_rows, collapse="\n"),
    "\n</tbody>\n</table></div>"
  )
}

summary_html <- make_html_table(
  table33[, c("Network","Bound","Pairs","Equivalent n","Equivalent %",
              "All 72 pass IUT"), drop=FALSE],
  "Exploratory equivalence summary by network and bound"
)

largest_html <- make_html_table(
  largest_table,
  "Largest absolute unregularised SMA-MRI partial correlation and 90% CI by network"
)

new_paragraph <- paste0(
  "<p>Exploratory equivalence analyses used uncorrected two one-sided tests ",
  "(TOSTs) at alpha = .05 for each of the 432 unregularised SMA-MRI partial ",
  "correlations, using bounds of plus or minus 0.05 and 0.03. Pairwise ",
  "equivalence corresponds to the 90% Fisher-z confidence interval falling ",
  "within the stated bound. Counts and percentages are reported for each ",
  "network and bound. For the network-level claim that all 72 SMA-MRI ",
  "associations fall within a bound, the intersection-union principle was used: ",
  "the claim is supported only if every pair passes at alpha = .05, with no ",
  "multiplicity correction required for that global claim. The bounds were ",
  "selected during revision and are not validated clinical thresholds.</p>"
)

holm_audit_note <- paste0(
  "<p>The Holm-across-432 columns in the detailed pairwise audit below are ",
  "retained only as a record of the earlier approach. That correction was too ",
  "conservative for the revised network-level all-pairs equivalence question ",
  "and is not used for the revised equivalence conclusion.</p>"
)

new_block <- paste(
  new_paragraph,
  summary_html,
  largest_html,
  holm_audit_note,
  sep="\n"
)

before <- if (start_pos > 1L) substr(html, 1L, start_pos - 1L) else ""
after <- substr(html, detail_start_pos, nchar(html))
updated <- paste0(before, new_block, "\n", after)

# Final safety checks before writing.
if (grepl("Holm correction is applied across 432 pairs separately at each bound",
          updated, fixed=TRUE)) {
  stop("Old Holm main-equivalence wording remains in the updated S2.")
}
if (!grepl("intersection-union principle", updated, fixed=TRUE)) {
  stop("New IUT wording was not inserted.")
}
if (!grepl("Exploratory equivalence results for every pair", updated, fixed=TRUE)) {
  stop("Detailed pairwise audit was lost; refusing to write.")
}

out_html <- file.path(outdir, "Supplementary_File_S2_VAS_IUT_UPDATED.html")
writeLines(updated, out_html, useBytes=TRUE)

message("")
message("COMPLETE")
message("Created:")
message("  ", out_html)
message("  ", file.path(outdir, "Supplementary_Table_33_VAS_IUT_UPDATED.csv"))
message("  ", file.path(outdir, "Supplementary_Table_33_Largest_r_90CI_VAS_UPDATED.csv"))
message("  ", file.path(outdir, "Supplementary_Table_33_VAS_IUT_UPDATED_NOTE.txt"))
message("")
message("Nothing else in S2 was intentionally changed.")
