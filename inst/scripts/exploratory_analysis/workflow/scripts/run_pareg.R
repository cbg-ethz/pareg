library(ggfittext)
library(tidyverse)

devtools::load_all()


# parameters
fname_de <- snakemake@input$fname_de
fname_terms <- snakemake@input$fname_terms
fname_sim <- snakemake@input$fname_sim

fname_enr <- snakemake@output$fname_enr
fname_obj <- snakemake@output$fname_obj

outdir <- dirname(fname_enr)

# load data
df_de <- read_tsv(fname_de) %>%
  dplyr::select(`SYMBOL...2`, PVAL) %>%
  dplyr::rename(gene = `SYMBOL...2`, pvalue = PVAL)
df_de %>%
  head()

df_terms <- read_csv(fname_terms)
df_terms %>%
  head()

term_similarities <- read.csv(fname_sim, row.names = 1, check.names = FALSE)
term_names <- df_terms %>%
  distinct(term) %>%
  pull(term)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix()

# run pareg
fit <- pareg::pareg(
  df_de,
  df_terms,
  term_network = term_similarities_sub,
  cv = TRUE,
  cv_method = "grid_search_lsf",
  family = pareg::beta_phi_var,
  lasso_param_range = seq(0, 1, length.out = 2),
  network_param_range = seq(0, 200, length.out = 5),
  tempdir = file.path(outdir, "cv_dump"),
  max_iteration = 10000,
  log_level = TRACE
)

df <- fit %>%
  as.data.frame() %>%
  arrange(desc(abs(enrichment)))

df %>%
  head()

# post-processing
if (fit$cv) {
  netreg_obj <- fit$obj$fit
} else {
  netreg_obj <- fit$obj
}

df_loss <- do.call(rbind, netreg_obj$loss_hist) %>%
  as_tibble() %>%
  unnest(everything()) %>%
  mutate(iteration = seq_along(netreg_obj$loss_hist))
df_loss %>%
  write_csv(file.path(outdir, "loss.csv"))

df_loss %>%
  pivot_longer(-iteration) %>%
  ggplot(aes(x = iteration, y = value, color = name)) +
  geom_line() +
  theme_minimal()
ggsave(
  file.path(outdir, "loss.pdf"),
  width = 8,
  height = 6
)

if ("loss_grid" %in% names(fit$obj)) {
  fit$obj$loss_grid %>%
    write_csv(file.path(outdir, "loss_grid.csv"))

  ggplot(fit$obj$loss_grid, aes(x = lambda, y = psigx, fill = loss)) +
    geom_tile() +
    geom_fit_text(aes(label = round(loss, 2)), color = "white") +
    scale_fill_viridis_c(direction = -1) +
    theme_minimal()
  ggsave(
    file.path(outdir, "loss_grid.pdf"),
    width = 8,
    height = 6
  )
}

# save result
df %>%
  write_csv(fname_enr)

write_rds(fit, fname_obj)
