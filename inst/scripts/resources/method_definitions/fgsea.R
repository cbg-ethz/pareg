# prepare environment
source(snakemake@params$setup_code_fname)

gene_list <- study$df %>%
  dplyr::select(gene, pvalue) %>%
  deframe()

term_list <- df_terms %>% {
  split(.$gene, .$term)
}

# run model
res <- fgsea::fgsea(
  pathways = term_list,
  stats = gene_list,
  scoreType = "std",
  nPermSimple = 10000
)

# extract results
df <- res %>%
  as.data.frame() %>%
  filter(!is.na(pval)) %>%
  mutate(enrichment = ifelse(
    pval > 0,
    -log10(pval),
    -log10(min(pval[pval > 0]))
  )) %>%
  rename(term = pathway) %>%
  select(term, enrichment)

# finalize
df <- df %>%
  mutate(
    method = "fgsea",
    term = str_to_lower(term)
  )

df %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
