# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
future::plan(future::multisession, worker = snakemake@threads)
fit <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  cv = TRUE,
  family = pareg::bernoulli,
  response_column_name = "pvalue_sig"
)
future::plan(future::sequential) # shut down workers

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_bernoulli", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# misc
pareg_post_processing(fit, dirname(snakemake@output$fname))

# save result
df %>%
  write_csv(snakemake@output$fname)
