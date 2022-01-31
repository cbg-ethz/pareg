# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
fit <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  network_param = 1, term_network = term_similarities_sub,
  family = netReg::bernoulli,
  response_column_name = "pvalue_sig"
)

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_network_ber", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# misc
do.call(rbind, fit$obj$loss_hist) %>%
  as_tibble() %>%
  unnest(everything()) %>%
  mutate(iteration = seq_along(fit$obj$loss_hist)) %>%
  pivot_longer(-iteration) %>%
  ggplot(aes(x = iteration, y = value, color = name)) +
    geom_line() +
    # geom_point() +
    theme_minimal()
ggsave(
  file.path(dirname(snakemake@output$fname), "loss.pdf"),
  width = 8,
  height = 6
)

# save result
df %>%
  write_csv(snakemake@output$fname)
