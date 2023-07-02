# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
fit <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  cv = TRUE,
  cv_method = pareg_cv_method,
  family = pareg::beta_phi_var,
  lasso_param_range = lasso_param_range,
  network_param_range = c(0),
  tempdir = file.path(dirname(snakemake@output$fname), "cv_dump"),
  max_iteration = pareg_max_iteration
)

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_lasso", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# misc
pareg_post_processing(fit, dirname(snakemake@output$fname))

# save result
df %>%
  write_csv(snakemake@output$fname)
