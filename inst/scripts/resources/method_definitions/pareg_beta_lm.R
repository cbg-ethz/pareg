# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
fit <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  term_network = term_similarities_sub,
  cv = TRUE,
  family = pareg::beta_phi_lm,
  lasso_param_range = lasso_param_range,
  network_param_range = network_param_range,
  log_level = logger::TRACE
)

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_beta_lm", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# misc
pareg_post_processing(fit, dirname(snakemake@output$fname))

# save result
df %>%
  write_csv(snakemake@output$fname)
