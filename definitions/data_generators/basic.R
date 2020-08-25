library(tidyverse)
library(graph)

set.seed(42)

# generate pathways which constitute cell
pathway1 <- dce::create_random_DAG(
  10, prob=.8, lB=c(-1, 0), uB=c(0, 1),
  node_labels = paste0("pw1_node", as.character(seq_len(10)))
)
pathway2 <- dce::create_random_DAG(
  10, prob=.8, lB=c(-1, 0), uB=c(0, 1),
  node_labels = paste0("pw2_node", as.character(seq_len(10)))
)

pathway2.mt <- dce::resample_edge_weights(pathway2, tp=1, mineff=2, maxeff=3)

# add background nodes
cell.bg <- as(matrix(0, nrow=10, ncol=10), "graphNEL")
nodes(cell.bg) <- paste0("bg_node", as.character(seq_len(10)))
edgemode(cell.bg) <- "directed"

# merge into cell
cell.wt <- dce::graph_union(c(pathway1, pathway2, cell.bg))
cell.mt <- dce::graph_union(c(pathway1, pathway2.mt, cell.bg))

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
