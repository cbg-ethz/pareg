library(tidyverse)


# parameters
fname_terms <- snakemake@input$fname_terms
fname_out <- snakemake@output$fname

plotdir <- snakemake@output$plotdir
dir.create(plotdir, recursive = TRUE)

category <- snakemake@params$params$category
subcategory <- snakemake@params$params$subcategory

# read data
df_terms <- read_csv(
  fname_terms,
  col_types = cols(
    gs_url = col_character(),
    gs_exact_source = col_character(),
    gs_geoid = col_character(),
    gs_pmid = col_character()
  )
) %>%
  filter(gs_cat == category) %>%
  {
    if (subcategory != "nan") {
      print("Filtering subcategory")
      filter(., gs_subcat == subcategory)
    } else {
      print("Skipping subcategory filter")
      .
    }
  }

# define distance measure
jaccard <- function(x, y) {
  return(1 - length(intersect(x, y)) / length(union(x, y)))
}

# compute term overlaps
term_list_list <- df_terms %>%
  select(gs_name, gene_symbol) %>%
  {
    split(.$gene_symbol, .$gs_name)
  }

term_similarities <- 1 - proxy::dist(
  x = term_list_list,
  method = jaccard,
  diag = TRUE, pairwise = TRUE
) %>%
  as.matrix()

# save result
term_similarities %>%
  write.csv(fname_out)

# visualize result
df_sim <- data.frame(similarity = term_similarities[upper.tri(term_similarities)])

df_sim %>%
  arrange(desc(similarity)) %>%
  head()

df_sim %>%
  ggplot(aes(x = similarity)) +
  geom_histogram(bins = 100) +
  xlim(0, 1) +
  scale_y_sqrt() +
  theme_minimal()
ggsave(file.path(plotdir, "similarity_histogram.pdf"))
