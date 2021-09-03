library(tidyverse)
library(glue)

library(plotROC)
library(cowplot)


# parameters
fname_list <- snakemake@input$fname_list

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- fname_list %>%
  map_dfr(function(path) {
    # TODO: make this better
    param_str <- gtools::split_path(path, depth_first = FALSE)[[3]]
    tmp <- list()
    for(param_pair in strsplit(param_str, "_")[[1]]) {
      parts <- strsplit(param_pair, "~")[[1]]
      tmp[parts[[1]]] <- parts[[2]]
    }
    read_csv(path) %>%
      mutate(!!!tmp)
  })

df_enr %>%
  head

# visualize
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
