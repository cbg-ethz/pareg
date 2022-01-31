library(tidyverse)


# parameters
fname_aucs <- snakemake@input$fname_aucs
fname_enr_list <- snakemake@input$fname_enr_list

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read loss data
df_loss <- fname_enr_list %>%
  map_dfr(function(path_enr) {
    path_loss <- file.path(dirname(path_enr), "loss.csv")

    path_parts <- gtools::split_path(path_loss, depth_first = FALSE)
    param_str <- path_parts[[length(path_parts) - 4]]
    tmp <- list()
    for (param_pair in strsplit(param_str, "_")[[1]]) {
      parts <- strsplit(param_pair, "~")[[1]]
      tmp[parts[[1]]] <- parts[[2]]
    }

    replicate <- path_parts[[length(path_parts) - 2]]
    method <- path_parts[[length(path_parts) - 1]]

    if (!"pareg" %in% method) {
      return(NULL)
    }

    read_csv(path_loss) %>%
      tail(1) %>% # only use final loss value
      mutate(path = path_loss) %>%
      rowwise %>%
      mutate(
        method = str_split_fixed(path, "/", Inf)[[7]],
        replicate = replicate,
        !!!tmp
      ) %>%
      select(-path)
  })

# read ROC-AUCs
df_aucs <- read_csv(
  fname_aucs,
  col_types = cols(.default = "c")
) %>%
  filter(str_detect(method, "pareg")) %>%
  mutate(auc = as.numeric(auc))

# combine data sources
df_res <- df_loss %>%
  inner_join(df_aucs)

df_res %>%
  head()

# create summary plot
p <- cowplot::plot_grid(
  ggplot(df_res, aes(x = method, y = total_loss)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal(),
  ggplot(df_res, aes(x = method, y = auc)) +
    geom_boxplot() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_minimal(),
  ncol = 1
)
p
cowplot::save_plot(file.path(outdir, "loss_vs_rocauc.pdf"), p, nrow = 2)
