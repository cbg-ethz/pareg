library(tidyverse)


# parameters
node_num <- 10 # rows, per pathway
pathway_num <- 3 # per category

sample_num <- 100 # columns

# helper functions
generate_count_matrix <- function(
  pathway_num, node_num, count_mean,
  perturbation_parameter, pathway_name_template
) {
  purrr::map_dfr(seq_len(pathway_num), function(i) {
    condition_size <- round(sample_num / 2)

    # compute wildtype counts
    wt_block <- matrix(
      rnbinom(
        node_num * condition_size,
        size = 100, mu = count_mean
      ),
      node_num, condition_size
    ) %>%
      as.data.frame
    colnames(wt_block) <- paste0("WT_", seq_len(condition_size))

    # compute mutant counts
    perturbation_vector <- if_else(
      runif(node_num * condition_size) > perturbation_parameter,
      3, 1
    )
    mt_block <- matrix(
      rnbinom(
        node_num * condition_size,
        size = 100, mu = count_mean * perturbation_vector
      ),
      node_num, condition_size
    ) %>%
      as.data.frame
    colnames(mt_block) <- paste0("MT_", seq_len(condition_size))

    # concat count blocks
    df_cts <- bind_cols(wt_block, mt_block)
    rownames(df_cts) <- paste0(
      glue::glue(pathway_name_template), seq_len(node_num)
    )

    df_cts
  })
}

# generate count data
count_mean <- 1000

df_cts <- bind_rows(
  generate_count_matrix(
    pathway_num, node_num, count_mean,
    .75, "pw{i}weak_node"
  ),
  generate_count_matrix(
    pathway_num, node_num, count_mean,
    .5, "pw{i}medium_node"
  ),
  generate_count_matrix(
    pathway_num, node_num, count_mean,
    .25, "pw{i}strong_node"
  ),

  generate_count_matrix(
    1, 500, count_mean,
    1, "bg{i}_node"
  )
)

# create dummy pathway topology
cell_wt <- matrix(0, dim(df_cts)[[1]], dim(df_cts)[[1]])
rownames(cell_wt) <- colnames(cell_wt) <- rownames(df_cts)

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
