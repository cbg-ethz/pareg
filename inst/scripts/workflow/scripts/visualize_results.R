library(tidyverse)
library(plotROC)
library(PRROC)


# parameters
fname_enr <- snakemake@input$fname_enr

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- read_csv(fname_enr)

df_enr %>%
  head()

# comparison plot
df_enr %>%
  ggplot(aes(x = is_on_term, y = enrichment, fill = method)) +
  geom_boxplot() +
  ylab("Enrichment measure") +
  facet_wrap(~method, scales = "free") +
  theme_minimal()
ggsave(file.path(outdir, "comparison.pdf"))

# ROC curves
df_enr %>%
  mutate(
    is_on_term = recode(as.character(is_on_term), "FALSE" = 0, "TRUE" = 1)
  ) %>%
  ggplot(aes(m = enrichment, d = is_on_term, color = method)) +
  geom_roc() +
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
  theme_minimal()
ggsave(file.path(outdir, "roc_curves.pdf"))

# other approach
cowplot::plot_grid(
  df_enr %>%
    group_by(method) %>%
    group_modify(function(group, key) {
      roc.curve(
        scores.class0 = group$enrichment,
        weights.class0 = group$is_on_term,
        curve = TRUE
      )$curve %>%
        as.data.frame()
    }) %>%
    rename(fpr = V1, tpr = V2, threshold = V3) %>%
    ggplot(aes(x = fpr, y = tpr, color = method)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
    ggtitle("ROC-curve") +
    theme_minimal(),
  df_enr %>%
    group_by(method) %>%
    group_modify(function(group, key) {
      pr.curve(
        scores.class0 = group$enrichment,
        weights.class0 = group$is_on_term,
        curve = TRUE
      )$curve %>%
        as.data.frame()
    }) %>%
    rename(recall = V1, precision = V2, threshold = V3) %>%
    ggplot(aes(x = recall, y = precision, color = method)) +
    geom_line() +
    ggtitle("PR-curve") +
    theme_minimal()
)
cowplot::save_plot(
  file.path(outdir, "performance_curves.pdf"), last_plot(),
  ncol = 2, nrow = 1
)
