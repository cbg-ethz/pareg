library(tidyverse)
library(enrichplot)

devtools::load_all()


# parameters
fname_enr <- snakemake@input$fname_enr
fname_obj <- snakemake@input$fname_obj

cancer_type <- snakemake@wildcards$cancer

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- read_csv(fname_enr)
fit <- read_rds(fname_obj)

df_enr %>%
  head()

# network plot
min_similarity <- 0.1
initial_term_count <- 50

enriched_terms <- df_enr %>%
  arrange(desc(abs(enrichment))) %>%
  head(initial_term_count) %>%
  pull(term)

mat_sub <- fit$term_network[enriched_terms, enriched_terms]
enriched_terms <- enriched_terms %>%
  map_dfr(function(x) {
    edge_weights <- mat_sub[x, !colnames(mat_sub) %in% c(x)]
    max_w <- max(edge_weights)

    data.frame(term = x, max_weight = max_w)
  }) %>%
  filter(max_weight >= min_similarity) %>%
  pull(term)
print(paste("#selected terms:", length(enriched_terms)))

plot(fit, term_subset = enriched_terms, min_similarity = min_similarity) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )

ggsave(
  file.path(outdir, glue("network_{cancer_type}.pdf")),
  width = 12,
  height = 12
)

# enrichplots
obj <- as_enrichplot_object(fit)

dotplot(obj, showCategory = 25) +
  scale_color_gradient2(
    low = "red",
    mid = "grey",
    high = "blue",
    midpoint = 0,
    na.value = "black",
    name = "Enrichment"
  ) +
  scale_size(range = c(3, 8), name = "Term size") +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )
ggsave(
  file.path(outdir, glue("dotplot_{cancer_type}.pdf")),
  width = 10,
  height = 11
)

treeplot(obj, showCategory = 30) +
  scale_colour_continuous(name = "Enrichment Score")
ggsave(
  file.path(outdir, glue("treeplot_{cancer_type}.pdf")),
  width = 15,
  height = 10
)

upplot <- upsetplot(obj, showCategory = 30)
upplot
ggsave(
  file.path(outdir, glue("upsetplot_{cancer_type}.pdf")),
  plot = upplot,
  width = 15,
  height = 10
)
