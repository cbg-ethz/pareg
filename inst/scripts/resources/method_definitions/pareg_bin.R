# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
df <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  family = netReg::poisson,
  response_column_name = "pvalue_notsig"
) %>%
  as.data.frame() %>%
  mutate(method = "pareg_bin", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
