# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
future::plan(future::multisession, worker = snakemake@threads)
fit <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  term_network = term_similarities_sub,
  cv = TRUE,
  family = pareg::bum
)
future::plan(future::sequential) # shut down workers

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_network_cv_bum", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
