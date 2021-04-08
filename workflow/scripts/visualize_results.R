library(tidyverse)
library(plotROC)


# parameters
fname_enr <- snakemake@input$fname_enr

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- read_csv(fname_enr)

df_enr %>%
  head

# comparison plot
df_enr %>%
  ggplot(aes(x = is_on_term, y = enrichment, fill = method)) +
    geom_boxplot() +
    ylab("Enrichment measure") +
    theme_minimal()
ggsave(file.path(outdir, "comparison.pdf"))

# ROC curves
df_enr %>%
  ggplot(aes(m = enrichment, d = is_on_term, color = method)) +
    geom_roc() +
    theme_minimal()
ggsave(file.path(outdir, "roc_curves.pdf"))
