library(tidyverse)


# parameters
fname_terms <- snakemake@input$fname_terms
fname_out <- snakemake@output$fname

# read data
df_terms <- read_csv(fname_terms, col_types = cols(gs_url = col_character()))

# define distance measure
jaccard <- function(x, y) {
  return(1 - length(intersect(x, y)) / length(union(x, y)))
}

# compute term overlaps
term_list_list <- df_terms %>%
  select(gs_name, gene_symbol) %>%
  { split(.$gene_symbol, .$gs_name) }

term_similarities <- 1 - proxy::dist(
  x = term_list_list,
  method = jaccard,
  diag = TRUE, pairwise = TRUE
) %>%
  as.matrix

# save result
term_similarities %>%
  write.csv(fname_out)
