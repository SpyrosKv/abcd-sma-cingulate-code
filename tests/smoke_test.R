# Base-R tests only. No ABCD data and no model estimation.
source("run_analysis.R")
nodes <- paste0("node", seq_len(18))
p <- matrix(NA_real_, 18, 18); p[upper.tri(p)] <- seq_len(153) / 154
e <- paper1_extract_edge_p(list(einv.pvals = p), nodes)
stopifnot(nrow(e) == 153L, !anyDuplicated(e[c("i", "j")]))
tab <- data.frame(Var1 = nodes[e$i], Var2 = nodes[e$j], check.names = FALSE)
tab[["p-value"]] <- e$p_raw
e2 <- paper1_extract_edge_p(list(einv.pvals = tab), nodes)
stopifnot(isTRUE(all.equal(e, e2, check.attributes = FALSE)))
bad <- tab; bad$Var1[1] <- "unknown_node"
stopifnot(inherits(try(paper1_extract_edge_p(list(einv.pvals = bad), nodes), silent = TRUE), "try-error"))
bad <- rbind(tab[-1, ], tab[2, ])
stopifnot(inherits(try(paper1_extract_edge_p(list(einv.pvals = bad), nodes), silent = TRUE), "try-error"))
stopifnot(inherits(try(paper1_authorised(list(authorised_data_use = FALSE)), silent = TRUE), "try-error"))
for (f in list.files("R", "[.]R$", recursive = TRUE, full.names = TRUE)) parse(file = f)
message("Base-R smoke tests passed. No participant data were read.")
