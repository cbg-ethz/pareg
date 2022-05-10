library(tidyverse)
library(enrichplot)

devtools::load_all()


# parameters
fname_enr <- snakemake@input$fname_enr
fname_obj <- snakemake@input$fname_obj

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- read_csv(fname_enr)
fit <- read_rds(fname_obj)

obj <- as_enrichplot_object(fit)

# plots
dotplot(obj, showCategory = 30) +
  scale_colour_continuous(name = "Enrichment Score")
ggsave(file.path(outdir, "dotplot.pdf"), width = 10, height = 15)

treeplot(obj, showCategory = 30) +
  scale_colour_continuous(name = "Enrichment Score")
ggsave(file.path(outdir, "treeplot.pdf"), width = 15, height = 10)

enriched_terms <- df_enr %>%
  arrange(desc(abs(enrichment))) %>%
  head(30) %>%
  pull(term)
plot(fit, term_subset = enriched_terms)
ggsave(file.path(outdir, "network.pdf"), width = 15, height = 15)

upplot <- upsetplot(obj, showCategory = 30)
upplot
ggsave(file.path(outdir, "upsetplot.pdf"), plot = upplot, width = 15, height = 10)
