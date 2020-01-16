library(tidyverse)


# read data
df.cts <- read_csv(snakemake@input$expr_fname) %>% column_to_rownames("node")
df.info <- read_csv(snakemake@input$info_fname) %>% column_to_rownames("sample")

# do DEA
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData=df.cts,
  colData=df.info,
  design=~condition
)
dds <- DESeq2::DESeq(dds)
df.res <- DESeq2::results(dds)

# store output
df.res %>%
  as.data.frame %>%
  rownames_to_column("node") %>%
  arrange(padj) %>%
  write_csv(snakemake@output$dea_fname)
