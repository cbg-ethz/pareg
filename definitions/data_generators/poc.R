library(tidyverse)


set.seed(42)

# helper functions
create <- function(num, lB, uB, pathway_name_template) {
  purrr::map(seq_len(num), function(i) {
    dce::create_random_DAG(
      10, prob=.8, lB=lB, uB=uB,
      node.labels = paste0(glue::glue(pathway_name_template), as.character(seq_len(10)))
    )
  })
}

resample <- function(graph_list, lB, uB) {
  purrr::map(graph_list, dce::resample_edge_weights, lB, uB)
}


# generate pathways which constitute cell
pathways.weak <- create(10, c(-0.1, 0), c(0, 0.1), "pw{i}weak_node")
pathways.medium <- create(10, c(-0.1, 0), c(0, 0.1), "pw{i}medium_node")
pathways.strong <- create(10, c(-0.1, 0), c(0, 0.1), "pw{i}strong_node")

pathways.weak.mt <- resample(pathways.weak, lB=c(-0.5, -0.1), uB=c(0.1, 0.5))
pathways.medium.mt <- resample(pathways.medium, lB=c(-1, -0.5), uB=c(0.5, 1))
pathways.strong.mt <- resample(pathways.strong, lB=c(-2, -1), uB=c(1, 2))

cell.wt <- dce::graph_union(c(pathways.weak, pathways.medium, pathways.strong))
cell.mt <- dce::graph_union(c(pathways.weak.mt, pathways.medium.mt, pathways.strong.mt))

# add background nodes
cell.bg <- as(matrix(0, nrow=10, ncol=10), "graphNEL")
nodes(cell.bg) <- paste0("bg_node", as.character(seq_len(10)))
edgemode(cell.bg) <- "directed"

cell.wt <- dce::graph_union(c(cell.wt, cell.bg))
cell.mt <- dce::graph_union(c(cell.mt, cell.bg))

# simulate data
X.wt <- dce::simulate_data(cell.wt, n=50)
X.mt <- dce::simulate_data(cell.mt, n=50)

# merge data
df.cts <- bind_cols(
  X.wt %>% t %>% as.data.frame %>% rename_all(list(~ paste0("WT_", .))),
  X.mt %>% t %>% as.data.frame %>% rename_all(list(~ paste0("MT_", .)))
)
rownames(df.cts) <- colnames(X.wt)

df.info <- data.frame(rowname=df.cts %>% colnames) %>%
  mutate(condition=rowname %>% str_extract("^..") %>% as.factor) %>%
  column_to_rownames("rowname")

stopifnot(all(rownames(df.info) == colnames(df.cts)))

# store data
fname.pathway <- commandArgs(TRUE)[1]
cell.wt %>%
  as("matrix") %>%
  as.data.frame %>%
  rownames_to_column("node") %>%
  write_csv(fname.pathway)

fname.expr <- commandArgs(TRUE)[2]
df.cts %>% rownames_to_column("node") %>% write_csv(fname.expr)

fname.info <- commandArgs(TRUE)[3]
df.info %>% rownames_to_column("sample") %>% write_csv(fname.info)
