# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
future::plan(future::multisession, worker = snakemake@threads)
fit <- pareg::pareg(
  study$df %>%
    select(gene, pvalue) %>%
    mutate(pvalue = -ifelse(
      pvalue > 0,
      log10(pvalue),
      log10(min(pvalue[pvalue > 0]))
    )),
  df_terms,
  term_network = term_similarities_sub,
  cv = TRUE,
  family = pareg::gaussian
)
future::plan(future::sequential) # shut down workers

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_ng", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# misc
pareg_post_processing(fit, dirname(snakemake@output$fname))

# save result
df %>%
  write_csv(snakemake@output$fname)
