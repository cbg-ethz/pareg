# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
df <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  truncate_response = TRUE
) %>%
  as.data.frame() %>%
  mutate(method = "pareg", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
