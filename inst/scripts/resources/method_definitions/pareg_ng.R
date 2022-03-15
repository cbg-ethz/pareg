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
saveRDS(fit, file.path(dirname(snakemake@output$fname), "fit.rds"))

df_loss <- do.call(rbind, fit$obj$fit$loss_hist) %>%
  as_tibble() %>%
  unnest(everything()) %>%
  mutate(iteration = seq_along(fit$obj$fit$loss_hist))
df_loss %>%
  write_csv(file.path(dirname(snakemake@output$fname), "loss.csv"))

df_loss %>%
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

ggplot(fit$obj$loss_grid, aes(x = lambda, y = psigx, fill = loss)) +
  geom_tile() +
  geom_text(aes(label = round(loss)), color = "white", size = 1) +
  theme_minimal()
ggsave(
  file.path(dirname(snakemake@output$fname), "loss_grid.pdf"),
  width = 8,
  height = 6
)

data.frame(
  pseudo_r_squared = fit$obj$fit$pseudo_r_squared
) %>%
  write_csv(file.path(dirname(snakemake@output$fname), "extra_stats.csv"))

# save result
df %>%
  write_csv(snakemake@output$fname)
