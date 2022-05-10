library(msigdbr)
library(tidyverse)
library(ComplexHeatmap)

devtools::load_all()


# parameters
fname_terms <- snakemake@output$fname_terms
fname_sim <- snakemake@output$fname_sim

# load data
df_terms <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:KEGG"
) %>%
  select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol) %>%
  distinct(.keep_all = TRUE)

df_terms %>%
  head()

# compute similarities
term_list_list <- df_terms %>%
  select(term, gene) %>%
  pipe_split("term", "gene")

term_similarities <- 1 - proxy::dist(
  x = term_list_list,
  method = function(x, y) {
    1 - jaccard(x, y)
  },
  diag = TRUE, pairwise = TRUE
) %>%
  as.matrix()

# plot summary
png(
  file.path(dirname(fname_sim), "similarity_clustermap.png"),
  width = 20,
  height = 20,
  units = "in",
  res = 300
)
Heatmap(
  term_similarities,
  name = "similarity",
  col = circlize::colorRamp2(c(0, 1), c("white", "black")),
  show_row_names = FALSE,
  show_column_names = FALSE
)
dev.off()

# save results
df_terms %>%
  write_csv(fname_terms)

term_similarities %>%
  write.csv(fname_sim)
