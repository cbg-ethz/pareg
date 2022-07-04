library(tidyverse)
library(glue)

library(PRROC)
library(plotROC)
library(cowplot)
library(scales)
library(latex2exp)


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

if ("pareg" %in% df_enr$method) {
  # put in first position
  df_enr <- df_enr %>%
    mutate(method = fct_relevel(method, "pareg"))
}
if ("random" %in% df_enr$method) {
  # put in last position
  df_enr <- df_enr %>%
    mutate(method = fct_relevel(method, "random", after = length(unique(df_enr$method))))
}

df_enr %>%
  head

# plot ROC/PR curves
df_enr %>%
  group_by(method) %>%
  group_walk(function(df_group, key) {
    print(key)
    parameter_columns <- setdiff(colnames(df_group), c("method", "term", "enrichment", "is_on_term", "replicate"))

    # ROC
    plot_list <- parameter_columns %>%
      map(function(param_name) {
        print("ROC")
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
          xlim(0, 1) +
          ylim(0, 1) +
          theme_minimal()
      })

    p <- plot_grid(plotlist = plot_list)
    save_plot(
      file.path(outdir, glue("rocs_{key}.pdf")),
      p,
      base_height = 10
    )

    # PR
    plot_list <- parameter_columns %>%
      map(function(param_name) {
        print("PR")
        print(param_name)
        df_group %>%
          mutate(
            replicate = as_factor(replicate), # fix hue plotting
          ) %>%
          group_by_at(vars(one_of(param_name))) %>%
          group_modify(function(group, key) {
            print(key)
            pr.curve(
              scores.class0 = group$enrichment,
              weights.class0 = group$is_on_term,
              curve = TRUE,
              minStepSize = 0.01
            )$curve %>%
              as.data.frame()
          }) %>%
          rename(recall = V1, precision = V2, threshold = V3) %>%
          ggplot(aes_string(x = "recall", y = "precision", color = param_name)) +
          geom_line() +
          ggtitle(glue("Parameter: {param_name}")) +
          xlim(0, 1) +
          ylim(0, 1) +
          theme_minimal()
      })

    p <- plot_grid(plotlist = plot_list)
    save_plot(
      file.path(outdir, glue("prs_{key}.pdf")),
      p,
      base_height = 10
    )
  })

# plot individual ROC/PR curves
roc_plot_dir <- file.path(outdir, "roc_plots")
dir.create(roc_plot_dir, showWarnings = FALSE, recursive = TRUE)

pr_plot_dir <- file.path(outdir, "pr_plots")
dir.create(pr_plot_dir, showWarnings = FALSE, recursive = TRUE)

df_enr %>%
  group_by(across(-all_of(c("term", "enrichment", "is_on_term", "replicate")))) %>%
  group_walk(function(df_group, key) {
    id_ <- apply(
      apply(key, 1, function(x) { n <- names(key); paste0(paste(n,x, sep = "~")) }),
      2,
      paste0, collapse = "_"
    )
    print(id_)

    df_group %>%
      mutate(
        replicate = as_factor(replicate), # fix hue plotting
        is_on_term = recode(as.character(is_on_term), "FALSE" = 0, "TRUE" = 1) # fix warning: "D not labeled 0/1, assuming FALSE = 0 and TRUE = 1"
      ) %>%
      ggplot(aes(m = enrichment, d = is_on_term, color = replicate)) +
      geom_roc() +
      geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
      ggtitle(id_) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(roc_plot_dir, glue("rocs_{id_}.pdf")), width = 8, height = 6)

    df_group %>%
      mutate(
        replicate = as_factor(replicate), # fix hue plotting
      ) %>%
      group_by(replicate) %>%
      group_modify(function(group, key) {
        print(key)
        pr.curve(
          scores.class0 = group$enrichment,
          weights.class0 = group$is_on_term,
          curve = TRUE,
          minStepSize = 0.01
        )$curve %>%
          as.data.frame()
      }) %>%
      rename(recall = V1, precision = V2, threshold = V3) %>%
      ggplot(aes(x = recall, y = precision, color = replicate)) +
      geom_line() +
      ggtitle("PR-curve") +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(pr_plot_dir, glue("prs_{id_}.pdf")), width = 8, height = 6)
  })

# plot individual ROC curves for subset of FPR values
roc_sub_plot_dir <- file.path(outdir, "roc_subset_plots")
dir.create(roc_sub_plot_dir, showWarnings = FALSE, recursive = TRUE)

df_enr %>%
  group_by(across(-all_of(c("term", "enrichment", "is_on_term", "replicate")))) %>%
  group_walk(function(df_group, key) {
    id_ <- apply(
      apply(key, 1, function(x) { n <- names(key); paste0(paste(n,x, sep = "~")) }),
      2,
      paste0, collapse = "_"
    )
    print(id_)

    df_group %>%
      mutate(
        replicate = as_factor(replicate), # fix hue plotting
        is_on_term = recode(as.character(is_on_term), "FALSE" = 0, "TRUE" = 1) # fix warning: "D not labeled 0/1, assuming FALSE = 0 and TRUE = 1"
      ) %>%
      ggplot(aes(m = enrichment, d = is_on_term, color = replicate)) +
      geom_roc() +
      geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
      ggtitle(id_) +
      scale_x_continuous(limits = c(0, 0.2), oob = squish) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(roc_sub_plot_dir, glue("rocs_{id_}.pdf")), width = 8, height = 6)
  })

# compute AUCs
df_auc <- df_enr %>%
  group_by(across(-all_of(c("term", "enrichment", "is_on_term")))) %>%
  group_modify(function(df_group, key) {
    roc_auc <- roc.curve(
      scores.class0 = df_group$enrichment,
      weights.class0 = df_group$is_on_term,
      curve = FALSE
    )$auc
    pr_auc <- pr.curve(
      scores.class0 = df_group$enrichment,
      weights.class0 = df_group$is_on_term,
      curve = FALSE
    )$auc.integral

    data.frame(roc_auc = roc_auc, pr_auc = pr_auc)
  })

df_auc %>%
  write_csv(fname_aucs)

df_auc %>%
  head()

# plot all AUCs
df_auc %>%
  ggplot(aes(x = method, y = roc_auc, fill = method)) +
  geom_boxplot() +
  geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
  xlab("Method") +
  ylab("ROC-AUC") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylim(0, 1) +
  theme_minimal()
ggsave(file.path(outdir, glue("roc_aucs.pdf")), width = 8, height = 6)

df_auc %>%
  ggplot(aes(x = method, y = pr_auc, fill = method)) +
  geom_boxplot() +
  geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
  xlab("Method") +
  ylab("PR-AUC") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylim(0, 1) +
  theme_minimal()
ggsave(file.path(outdir, glue("pr_aucs.pdf")), width = 8, height = 6)

# plot individual AUCs
auc_plot_dir <- file.path(outdir, "auc_plots")
dir.create(auc_plot_dir, showWarnings = FALSE, recursive = TRUE)

df_auc %>%
  group_by(across(-all_of(c("method", "replicate", "roc_auc", "pr_auc")))) %>%
  group_walk(function(df_group, key) {
    print(key)

    id_ <- apply(
      apply(key, 1, function(x) { n <- names(key); paste0(paste(n,x, sep = "~")) }),
      2,
      paste0, collapse = "_"
    )

    ggplot(df_group, aes(x = method, y = roc_auc, fill = method)) +
      geom_boxplot() +
      geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
      xlab("Method") +
      ylab("ROC-AUC") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(auc_plot_dir, glue("roc_aucs_{id_}.pdf")), width = 8, height = 6)

    ggplot(df_group, aes(x = method, y = pr_auc, fill = method)) +
      geom_boxplot() +
      geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
      xlab("Method") +
      ylab("PR-AUC") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(auc_plot_dir, glue("pr_aucs_{id_}.pdf")), width = 8, height = 6)
  })

# more final plots
parameter_columns <- setdiff(
  colnames(df_enr),
  c("method", "term", "enrichment", "is_on_term", "replicate")
)

parameter_columns %>%
  walk(function(param_name) {
    enr_grpd <- df_enr %>%
      group_by_at(vars(one_of(param_name, "method"))) %>%
      group_modify(function(group, key) {
        print(key)
        pr.curve(
          scores.class0 = group$enrichment,
          weights.class0 = group$is_on_term,
          curve = TRUE,
          minStepSize = 0.01
        )$curve %>%
          as.data.frame()
      }) %>%
      rename(recall = V1, precision = V2, threshold = V3)

    enr_grpd %>%
      ggplot(aes_string(
        x = "recall",
        y = "precision",
        color = param_name,
        linetype = "method"
      )) +
      geom_line() +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(outdir, glue("pr_aggregated_{param_name}.pdf")), width = 8, height = 6)

    enr_grpd %>%
      group_by_at(vars(one_of(param_name))) %>%
      group_walk(function(group, key) {
        print(key)

        group %>%
          ggplot(aes(
            x = recall,
            y = precision,
            color = method
          )) +
          geom_line() +
          xlim(0, 1) +
          ylim(0, 1) +
          theme_minimal()
        ggsave(file.path(outdir, glue("pr_aggregated_{param_name}_{key}.pdf")), width = 8, height = 6)
      })
  })

parameter_columns %>%
  walk(function(param_name) {
    df_auc %>%
      ggplot(aes_string(x = param_name, y = "pr_auc", fill = "method")) +
      geom_boxplot(outlier.colour = NA) +
      geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
      xlab(TeX(case_when(
        param_name == "beta" ~ "Noise level \\eta",  # \eta is not a typo
        param_name == "similarityfactor" ~ "Similarity factor \\rho",
        TRUE ~ param_name
      ))) +
      ylab("PR-AUC") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      ylim(0, 1) +
      theme_minimal()
    ggsave(file.path(outdir, glue("pr_aucs_aggregated_{param_name}.pdf")), width = 8, height = 6)
  })

# comparison for given precision/recall values
recall_value <- 0.8
precision_value <- 0.8

parameter_columns %>%
  walk(function(param_name) {
    # gather data
    enr_grpd <- df_enr %>%
      group_by_at(vars(one_of(param_name, "method", "replicate"))) %>%
      group_modify(function(group, key) {
        pr.curve(
          scores.class0 = group$enrichment,
          weights.class0 = group$is_on_term,
          minStepSize = 0.01,
          curve = TRUE
        )$curve %>%
          as.data.frame()
      }) %>%
      rename(recall = V1, precision = V2, threshold = V3)

    # compute threshold values
    selection_threshold <- 0.01
    recall_selection <- enr_grpd %>%
      filter(abs(recall - recall_value) < selection_threshold) %>%
      slice_min(abs(recall - recall_value), with_ties = FALSE)
    precision_selection <- enr_grpd %>%
      filter(abs(precision - precision_value) < selection_threshold) %>%
      slice_min(abs(precision - precision_value), with_ties = FALSE)

    # `recall_selection` will always have all entries, but `precision_selection`
    # might lack some. Let's add them
    precision_selection_missing <- recall_selection %>%
      distinct(method, replicate) %>%
      pmap_dfr(function (...) {
        row <- list(...)

        df_match <- plyr::match_df(
          precision_selection,
          data.frame(
            method = row$method,
            replicate = row$replicate
          ),
          on = c("method", "replicate")
        )

        if (dim(df_match)[[1]] == 0) {
          # entry in `recall_selection` is missing from `precision_selection`
          new_row <- data.frame(
            method = row$method,
            replicate = row$replicate,
            recall = 0,
            precision = precision_value,
            threshold = NA
          )
          new_row[[param_name]] <- row[[param_name]]

          return(new_row)
        }

        NULL
      })

    precision_selection <- bind_rows(
      precision_selection,
      precision_selection_missing
    ) %>%
      arrange(method, replicate)

    # recall plot
    if (dim(recall_selection)[[1]]) {
      print(glue("Recall summary for {param_name}"))
      recall_selection %>%
        group_by_at(vars(one_of(param_name, "method"))) %>%
        summarize(
          precision_mean = mean(precision),
          precision_median = median(precision),
          precision_sd = sd(precision)
        ) %>%
        print()

      recall_selection %>%
        ggplot(aes_string(x = param_name, y = "precision", fill = "method")) +
        geom_boxplot(outlier.colour = NA) +
        geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
        xlab(TeX(case_when(
          param_name == "beta" ~ "Noise level \\eta",  # \eta is not a typo
          param_name == "similarityfactor" ~ "Similarity factor \\rho",
          TRUE ~ param_name
        ))) +
        ylab("Precision") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        ylim(0, 1) +
        theme_minimal()
      ggsave(file.path(outdir, glue("pr_recall_threshold_{param_name}_{recall_value}.pdf")), width = 8, height = 6)
    }

    # precision plot
    if (dim(precision_selection)[[1]]) {
      print(glue("Precision summary for {param_name}"))
      precision_selection %>%
        group_by_at(vars(one_of(param_name, "method"))) %>%
        summarize(
          recall_mean = mean(recall),
          recall_median = median(recall),
          recall_sd = sd(recall)
        ) %>%
        print()

      precision_selection %>%
        ggplot(aes_string(x = param_name, y = "recall", fill = "method")) +
        geom_boxplot(outlier.colour = NA) +
        geom_point(color = "black", size = 0.4, alpha = 0.9, position = position_jitterdodge()) +
        xlab(TeX(case_when(
          param_name == "beta" ~ "Noise level \\eta",  # \eta is not a typo
          param_name == "similarityfactor" ~ "Similarity factor \\rho",
          TRUE ~ param_name
        ))) +
        ylab("Recall") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        ylim(0, 1) +
        theme_minimal()
      ggsave(file.path(outdir, glue("pr_precision_threshold_{param_name}_{precision_value}.pdf")), width = 8, height = 6)
    }
  })
