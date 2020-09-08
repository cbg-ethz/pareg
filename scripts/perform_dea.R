library(tidyverse)


# read data
df_cts <- read_csv(snakemake@input$expr_fname) %>% column_to_rownames("node")
df_info <- read_csv(snakemake@input$info_fname) %>% column_to_rownames("sample")

out_dir <- snakemake@output$out_dir

stopifnot(all(rownames(df_info) == colnames(df_cts)))

# do DEA
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = df_cts,
  colData = df_info,
  design = ~ condition
)
dds <- DESeq2::DESeq(dds)
df_res <- DESeq2::results(dds)

# sanity checks
gene_list <- c(
  "bg1_node1",
  "pw1weak_node1", "pw1medium_node1", "pw1strong_node1"
)

purrr::map_dfr(gene_list, function(gene) {
  DESeq2::plotCounts(
    dds, gene = gene, intgroup = "condition", returnData = TRUE
  ) %>%
    mutate(gene = gene)
}) %>%
  mutate(
    condition = factor(condition, levels = c("WT", "MT")),
    gene = factor(gene, levels = gene_list)
  ) %>%
  ggplot(aes(x = condition, y = count)) +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  facet_wrap(~ gene) +
  theme_minimal()

ggsave(file.path(out_dir, "normalized_counts.pdf"))

# store output
df_res %>%
  as.data.frame %>%
  rownames_to_column("node") %>%
  arrange(padj) %>%
  write_csv(snakemake@output$dea_fname)
