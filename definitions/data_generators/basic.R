library(tidyverse)


# generate pathways which constitute cell
pathway1 <- dce::create_random_DAG(
    10, prob=.8, lB=c(0.5, 1), uB=c(0.5, 1),
    node.labels = paste0("pw1_node", as.character(seq_len(10)))
)
pathway2 <- dce::create_random_DAG(
    10, prob=.8, lB=c(0.5, 1), uB=c(0.5, 1),
    node.labels = paste0("pw2_node", as.character(seq_len(10)))
)

pathway2.mt <- dce::resample_edge_weights(pathway2, lB=c(1.5, 2), uB=c(1.5, 2))

cell.wt <- graph_union(pathway1, pathway2)
cell.mt <- graph_union(pathway1, pathway2.mt)

# simulate data
X.wt <- dce::simulate_data_better(cell.wt, n=10)
X.mt <- dce::simulate_data_better(cell.mt, n=10)

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
