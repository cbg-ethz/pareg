library(tidyverse)
library(glue)

library(PRROC)
library(plotROC)
library(cowplot)


# parameters
fname_list <- snakemake@input$fname_list

fname_aucs <- snakemake@output$fname_aucs
outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- fname_list %>%
  map_dfr(function(path) {
    # TODO: make this better
    path_parts <- gtools::split_path(path, depth_first = FALSE)
    param_str <- path_parts[[length(path_parts) - 1]]
    tmp <- list()
    for (param_pair in strsplit(param_str, "_")[[1]]) {
      parts <- strsplit(param_pair, "~")[[1]]
      tmp[parts[[1]]] <- parts[[2]]
    }
    read_csv(path) %>%
      mutate(!!!tmp)
  })

df_enr %>%
  head

# plot ROC curves
df_enr %>%
  group_by(method) %>%
  group_walk(function(df_group, key) {
    print(key)
    parameter_columns <- setdiff(colnames(df_group), c("method", "term", "enrichment", "is_on_term", "replicate"))

    plot_list <- parameter_columns %>%
      map(function(param_name) {
        print(param_name)
        df_group %>%
          mutate(
            replicate = as_factor(replicate), # fix hue plotting
            is_on_term = recode(as.character(is_on_term), "FALSE" = 0, "TRUE" = 1) # fix warning: "D not labeled 0/1, assuming FALSE = 0 and TRUE = 1"
          ) %>%
          ggplot(aes_string(m = "enrichment", d = "is_on_term", color = param_name)) +
          geom_roc() +
          geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
          ggtitle(glue("Parameter: {param_name}")) +
          theme_minimal()
      })

    p <- plot_grid(plotlist = plot_list)
    save_plot(
      file.path(outdir, glue("rocs_{key}.pdf")),
      p,
      base_height = 10
    )
  })

# compute and plot AUCs
df_auc <- df_enr %>%
  group_by(method, replicate, termsource, alpha, beta, similarityfactor, ontermcount, siggenescaling) %>%
  group_modify(function(df_group, key) {
    auc <- roc.curve(
      scores.class0 = df_group$enrichment,
      weights.class0 = df_group$is_on_term,
      curve = FALSE
    )$auc

    data.frame(auc = auc)
  })

df_auc %>%
  write_csv(fname_aucs)

df_auc %>%
  head()

df_auc %>%
  ggplot(aes(x = method, y = auc, fill = method)) +
  geom_boxplot() +
  xlab("Method") +
  ylab("ROC-AUC") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_minimal()
ggsave(file.path(outdir, glue("roc_aucs.pdf")), width = 8, height = 6)
