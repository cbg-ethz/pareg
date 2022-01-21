library(tidyverse)


# parameters
fname_terms <- "results/term_database.csv"
category <- "C2"


fname_terms <- snakemake@input$fname_terms
fname_out <- snakemake@output$fname

plotdir <- snakemake@output$plotdir
dir.create(plotdir, recursive = TRUE)

# read data
df_terms <- read_csv(fname_terms)

# define distance measure
jaccard <- function(x, y) {
  return(1 - length(intersect(x, y)) / length(union(x, y)))
}

# compute term overlaps
term_list_list <- df_terms %>%
  select(term, gene) %>%
  {
    split(.$gene, .$term)
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

# simlarity histogram
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

# similarity clustermap
max_size <- 500
if (nrow(term_similarities) > max_size) {
  # subsample so creating clustermap doesn't take too long
  random_indices <- sample(nrow(term_similarities), max_size)
  term_similarities_subsample <- term_similarities[random_indices, random_indices]
} else {
  term_similarities_subsample <- term_similarities
}

png(
  file.path(plotdir, "similarity_clustermap.png"),
  width = 20,
  height = 20,
  units = "in",
  res = 300
)
ComplexHeatmap::Heatmap(
  term_similarities_subsample,
  name = "similarity",
  col = circlize::colorRamp2(c(0, 1), c("white", "black"))
)
dev.off()
