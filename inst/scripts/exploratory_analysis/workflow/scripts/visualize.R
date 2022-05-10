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

# dotplot
dotplot(obj) +
  scale_colour_continuous(name = "Enrichment Score")
ggsave(file.path(outdir, "dotplot.pdf"))
