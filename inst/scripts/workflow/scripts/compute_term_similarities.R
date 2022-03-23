library(tidyverse)

devtools::load_all("../..")


# parameters
fname_terms <- snakemake@input$fname_terms
similarity_measure <- snakemake@params$params$similaritymeasure
fname_out <- snakemake@output$fname

plotdir <- snakemake@output$plotdir
dir.create(plotdir, recursive = TRUE)

# read data
df_terms <- read_csv(fname_terms)

# select similarity function
similarity_function <- match.fun(similarity_measure, descend = FALSE)
similarity_function

# compute term similarities
term_list_list <- df_terms %>%
  select(term, gene) %>%
  pipe_split("term", "gene")

term_similarities <- 1 - proxy::dist(
  x = term_list_list,
  method = function(x, y) {
    1 - similarity_function(x, y)
  },
  diag = TRUE, pairwise = TRUE
) %>%
  as.matrix()

# save result
term_similarities %>%
  write.csv(fname_out)

# similarity histogram
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
