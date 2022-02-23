library(tidyverse)
devtools::load_all()


# parameters
fname_aucs <- snakemake@input$fname_aucs
fname_enr_list <- snakemake@input$fname_enr_list
fname_study_list <- snakemake@input$fname_study_list

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read loss data
df_loss_test <- pmap_dfr(
  list(fname_enr_list, fname_study_list),
  function(path_enr, path_study) {
    # parse parameters
    path_parts <- gtools::split_path(path_enr, depth_first = FALSE)
    param_str <- path_parts[[length(path_parts) - 4]]
    tmp <- list()
    for (param_pair in strsplit(param_str, "_")[[1]]) {
      parts <- strsplit(param_pair, "~")[[1]]
      tmp[parts[[1]]] <- parts[[2]]
    }

    replicate <- path_parts[[length(path_parts) - 2]]
    method <- path_parts[[length(path_parts) - 1]]

    if (str_detect(method, "pareg", negate = TRUE)) {
      return(NULL)
    }

    # read data
    fit <- readRDS(file.path(dirname(path_enr), "fit.rds"))
    mod <- keras::load_model_tf(fit$obj$model)
    study <- readRDS(path_study)

    # setup model inputs
    df_model <- pareg::create_model_df(study$df %>% select(gene, pvalue), fit$df_terms)
    covariates <- df_model %>%
      select(ends_with(".member")) %>%
      names()
    X_test <- df_model %>%
      select(all_of(covariates)) %>%
      as.matrix()
    Y_test <- df_model %>%
      select(fit$response_column_name) %>%
      as.matrix()

    # compute loss scores on test samples
    loss_test <- edgenet.loss(
      fit$params$lasso_param, fit$params$network_param, 0,
      cast_float(fit$term_network), NULL,
      fit$obj$family
    )(mod, cast_float(X_test), cast_float(transform_y(Y_test)))

    # finalize result
    data.frame(likelihood_test = loss_test$likelihood) %>%
      rowwise %>%
      mutate(
        method = method,
        replicate = replicate,
        !!!tmp
      )
})

df_loss <- fname_enr_list %>%
  map_dfr(function(path_enr) {
    # parse parameters
    path_parts <- gtools::split_path(path_enr, depth_first = FALSE)
    param_str <- path_parts[[length(path_parts) - 4]]
    tmp <- list()
    for (param_pair in strsplit(param_str, "_")[[1]]) {
      parts <- strsplit(param_pair, "~")[[1]]
      tmp[parts[[1]]] <- parts[[2]]
    }

    replicate <- path_parts[[length(path_parts) - 2]]
    method <- path_parts[[length(path_parts) - 1]]

    if (str_detect(method, "pareg", negate = TRUE)) {
      return(NULL)
    }

    # parse extra statistics
    path_stats <- file.path(dirname(path_enr), "extra_stats.csv")
    df_stats <- read_csv(path_stats)

    # parse loss
    path_loss <- file.path(dirname(path_enr), "loss.csv")

    df_lo <- read_csv(path_loss) %>%
      tail(1) # only use final loss value

    # assemble result
    df_lo %>%
      mutate(path = path_enr) %>%
      rowwise %>%
      mutate(
        method = str_split_fixed(path, "/", Inf)[[7]],
        replicate = replicate,
        !!!tmp,
        !!!df_stats
      ) %>%
      select(-path)
  })

# read ROC-AUCs
df_aucs <- read_csv(
  fname_aucs,
  col_types = cols(.default = "c")
) %>%
  filter(str_detect(method, "pareg")) %>%
  mutate(
    roc_auc = as.numeric(roc_auc),
    pr_auc = as.numeric(pr_auc)
  )

# combine data sources
df_res <- df_loss %>%
  inner_join(df_aucs) %>%
  inner_join(df_loss_test)

df_res %>%
  head()

# create summary plot
p <- cowplot::plot_grid(
  ggplot(df_res, aes(x = method, y = likelihood)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = likelihood_test)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = pseudo_r_squared)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ylim(0, 1) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = mse)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = roc_auc)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ylim(0, 1) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = pr_auc)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ylim(0, 1) +
    theme_minimal(),
  ncol = 1
)
p
cowplot::save_plot(file.path(outdir, "loss_vs_rocauc.pdf"), p, nrow = 3)

p <- cowplot::plot_grid(
  ggplot(df_res, aes(x = method, y = likelihood_test)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = roc_auc)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ylim(0, 1) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = pr_auc)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ylim(0, 1) +
    theme_minimal(),
  ncol = 1
)
p
cowplot::save_plot(file.path(outdir, "loss_vs_rocauc_subset.pdf"), p, nrow = 3)
