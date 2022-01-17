# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
fit <- pareg::pareg(
  study$df %>% select(gene, pvalue),
  df_terms,
  term_network = term_similarities_sub,
  truncate_response = TRUE,
  cv = TRUE,
  family = netReg::bum
)

df <- fit %>%
  as.data.frame() %>%
  mutate(method = "pareg_network_cv_bum", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# misc
data.frame(
  iteration = seq_along(fit$obj$loss_hist),
  loss = fit$obj$loss_hist
) %>%
  ggplot(aes(x = iteration, y = loss)) +
    geom_line() +
    geom_point() +
    theme_minimal()
ggsave(
  file.path(dirname(snakemake@output$fname), "loss.pdf"),
  width = 8,
  height = 6
)

# save result
df %>%
  write_csv(snakemake@output$fname)
