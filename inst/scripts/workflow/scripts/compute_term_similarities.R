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
if (similarity_measure == "overlap") {
  # fix till https://github.com/snakemake/snakemake/pull/1299 lands
  similarity_function <- overlap_coefficient
} else if (similarity_measure == "semantic") {
  similarity_function <- function(term_list) {
    # convert to upper case because workflow assumes lower case
    term_list_orig <- str_replace_all(str_to_upper(term_list), "_", ":")

    # actually compute similarities
    mat <- simplifyEnrichment::GO_similarity(term_list_orig, measure = "Rel")

    # fix term names
    stopifnot(all(rownames(mat) == colnames(mat)))
    rownames(mat) <- str_replace_all(str_to_lower(rownames(mat)), ":", "_")
    colnames(mat) <- str_replace_all(str_to_lower(colnames(mat)), ":", "_")

    # include terms which have no associated semantic similarity
    # there must be a better way?!
    current_rownames <- rownames(mat)
    for (term in setdiff(term_list, current_rownames)) {
      print(paste("Adding missing GO term", term))

      mat <- rbind(cbind(mat, foo = 0), foo = 0)

      rownames(mat)[length(rownames(mat))] <- term
      colnames(mat)[length(colnames(mat))] <- term
    }

    # make sure matrix is ordered correctly
    stopifnot(setequal(rownames(mat), term_list))
    mat <- mat[term_list, term_list]

    mat
  }
} else {
  similarity_function <- match.fun(similarity_measure, descend = FALSE)
}

# compute term similarities
term_list_list <- df_terms %>%
  select(term, gene) %>%
  pipe_split("term", "gene")

if (similarity_measure == "semantic") {
  term_similarities <- similarity_function(names(term_list_list))
} else {
  term_similarities <- 1 - proxy::dist(
    x = term_list_list,
    method = function(x, y) {
      1 - similarity_function(x, y)
    },
    diag = TRUE, pairwise = TRUE
  ) %>%
    as.matrix()
}

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
  col = circlize::colorRamp2(c(0, 1), c("white", "black")),
  show_row_names = FALSE,
  show_column_names = FALSE
)
dev.off()
