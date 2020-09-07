library(tidyverse)


# parameters
node_num <- 100 # rows
sample_num <- 100 # columns

# create dummy pathway topology
cell_wt <- matrix(0, node_num, node_num)

# generate count data
count_mean <- 1000

df_cts <- matrix(
  rnbinom(node_num * sample_num, size = 100, mu = count_mean),
  node_num, sample_num
) %>%
  as.data.frame

rownames(df_cts) <- paste0("node_", 1:node_num)
colnames(df_cts) <- c(paste0("WT_", 1:50), paste0("MT_", 1:50))

# simulate pathway enrichment
condition_size <- round(sample_num / 2)
wt_range <- seq(condition_size + 1, sample_num)

df_cts[1:10, wt_range] <- rnbinom(
  10 * condition_size,
  size = 100, mu = count_mean * if_else(runif(100) > .75, 3, 1)
)
df_cts[11:20, wt_range] <- rnbinom(
  10 * condition_size,
  size = 100, mu = count_mean * if_else(runif(100) > .5, 3, 1)
)
df_cts[21:30, wt_range] <- rnbinom(
  10 * condition_size,
  size = 100, mu = count_mean * if_else(runif(100) > .25, 3, 1)
)

rownames(df_cts) <- c(
  paste0("pw0weak_", rownames(df_cts[1:10, ])),
  paste0("pw0medium_", rownames(df_cts[11:20, ])),
  paste0("pw0strong_", rownames(df_cts[21:30, ])),
  paste0("bg_", rownames(df_cts[31:node_num, ]))
)

rownames(cell_wt) <- rownames(df_cts)
colnames(cell_wt) <- rownames(df_cts)

# store sample information
df_info <- data.frame(rowname = df_cts %>% colnames) %>%
  mutate(condition = rowname %>% str_extract("^..") %>% as.factor) %>%
  column_to_rownames("rowname")

stopifnot(all(rownames(df_info) == colnames(df_cts)))

# summary
# df_cts %>%
#   as.matrix %>%
#   reshape2::melt(.) %>%
# ggplot(aes(x = Var2, y = Var1, fill = value)) +
#   geom_tile()

# store data
fname_pathway <- commandArgs(TRUE)[1]
cell_wt %>%
  as.data.frame %>%
  rownames_to_column("node") %>%
  write_csv(fname_pathway)

fname_expr <- commandArgs(TRUE)[2]
df_cts %>% rownames_to_column("node") %>% write_csv(fname_expr)

fname_info <- commandArgs(TRUE)[3]
df_info %>% rownames_to_column("sample") %>% write_csv(fname_info)
