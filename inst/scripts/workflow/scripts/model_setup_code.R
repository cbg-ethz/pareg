library(tidyverse)
library(ggfittext)


# parameters
fname_study <- snakemake@input$fname_study
fname_terms <- snakemake@input$fname_terms
fname_term_sim <- snakemake@input$fname_term_sim

fname_out <- snakemake@output$fname

# read data
study <- readRDS(fname_study)

df_terms <- read_csv(fname_terms)

df_terms %>% dim()
df_terms %>% head()

term_similarities <- read.csv(fname_term_sim, row.names = 1, check.names = FALSE)

# prepare data
study_genes <- study$df %>%
  filter(pvalue <= 0.05) %>%
  pull(gene)

nonstudy_genes <- study$df %>%
  filter(pvalue > 0.05) %>%
  pull(gene)

all_genes <- df_terms %>%
  distinct(gene) %>%
  pull(gene)
# all_genes <- c(study_genes, nonstudy_genes)

term_names <- df_terms %>%
  distinct(term) %>%
  pull(term)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix()

# helper functions
pareg_post_processing <- function(fit, outdir) {
  # save whole fit object
  saveRDS(fit, file.path(outdir, "fit.rds"))

  # extract netreg object based on whether CV was done or not
  if (fit$cv) {
    netreg_obj <- fit$obj$fit
  } else {
    netreg_obj <- fit$obj
  }

  # loss trajectory
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
      # geom_point() +
      theme_minimal()
  ggsave(
    file.path(outdir, "loss.pdf"),
    width = 8,
    height = 6
  )

  # loss grid for CV calls
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

  # performance statistics
  df_stats <- data.frame(
    pseudo_r_squared = netreg_obj$pseudo_r_squared
  )

  if (!is.null(netreg_obj$mse)) {
    df_stats$mse <- netreg_obj$mse
  }

  df_stats %>%
    write_csv(file.path(outdir, "extra_stats.csv"))
}
