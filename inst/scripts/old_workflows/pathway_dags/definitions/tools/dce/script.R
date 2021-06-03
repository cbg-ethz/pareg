library(tidyverse)
library(dce)

# read data
fname.cell <- commandArgs(TRUE)[1]
df.cell <- read_csv(fname.cell) %>% column_to_rownames("node")

fname.expr <- commandArgs(TRUE)[2]
df.expr <- read_csv(fname.expr) %>% column_to_rownames("node")

fname.info <- commandArgs(TRUE)[3]
df.info <- read_csv(fname.info)

# find pathways
pw.map <- list()
for (node in rownames(df.expr)) {
  pw <- strsplit(node, "_")[[1]][1]
  pw.map[[pw]] <- c(pw.map[[pw]], node)
}

# compute enrichment
purrr::map_dfr(pw.map, function (nodes) {
  graph <- igraph::igraph.to.graphNEL(
    igraph::graph_from_adjacency_matrix(
      df.cell[nodes, nodes] %>% as.matrix,
      weighted=TRUE
    )
  )

  X.wt <- df.expr[
    nodes,
    df.info %>%
      dplyr::filter(condition == "WT") %>%
      pull(sample)
  ]
  X.mt <- df.expr[
    nodes,
    df.info %>%
      dplyr::filter(condition == "MT") %>%
      pull(sample)
  ]

  pval <- dce::dce_nb(graph, t(X.wt), t(X.mt), p_method = "meanp")$pathway_pvalue

  ## this does throw up an error in the solver:
  ## Error in while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - Lm)/d1 +  :
  ## missing value where TRUE/FALSE needed
  ## pvals <- NULL
  ## dn <- t(X.wt)
  ## dt <- t(X.mt)
  ## pruns <- 100
  ## for (i in seq_len(pruns)) {
  ##     dnp <- dn
  ##     colnames(dnp) <- sample(colnames(dn), ncol(dn))
  ##     dtp <- dt
  ##     colnames(dtp) <- sample(colnames(dt), ncol(dt))
  ##     dnp <- dnp[, naturalsort::naturalorder(colnames(dnp))]
  ##     dtp <- dtp[, naturalsort::naturalorder(colnames(dtp))]
  ##     pvalp <- try(dce::dce_nb(graph, dnp, dtp)$pathway_pvalue)
  ##     if (length(grep("Error", pvalp)) > 0) {
  ##         save(dnp, dtp, graph, file = "temp.rda")
  ##     } else {
  ##         pvals <- c(pvals, pvalp)
  ##     }
  ## }
  ## pval <- sum(pvals <= pval)/pruns

  data.frame(pvalue=pval)
}, .id="pathway") %>%
  write_csv("result.csv")
