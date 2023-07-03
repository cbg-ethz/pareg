library(tidyverse)
library(ggsignif)
library(plotROC)
library(PRROC)


# parameters
fname_enr <- snakemake@input$fname_enr
fname_benchmark <- snakemake@input$fname_benchmark

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# read data
df_enr <- read_csv(fname_enr)
df_enr %>%
  head()

df_runtime <- read_csv(fname_benchmark)
df_runtime %>%
  head()

# comparison plot
df_enr %>%
  ggplot(aes(x = is_on_term, y = enrichment)) +
  geom_boxplot() +
  geom_jitter(shape = ".") +
  geom_signif(comparisons = list(c("FALSE", "TRUE"))) +
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
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) +
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
    xlim(0, 1) +
    ylim(0, 1) +
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
    xlim(0, 1) +
    ylim(0, 1) +
    theme_minimal()
)
cowplot::save_plot(
  file.path(outdir, "performance_curves.pdf"), last_plot(),
  ncol = 2, nrow = 1
)


# pairwise scatterplots
pairwise_dir <- file.path(outdir, "pairwise")
dir.create(pairwise_dir, showWarnings = FALSE, recursive = TRUE)

df_enr %>%
  group_by(replicate) %>%
  group_walk(function(df_group, key) {
    p <- df_group %>%
      pivot_wider(
        id_cols = c("term", "is_on_term"),
        names_from = c("method"),
        values_from = c("enrichment")
      ) %>%
    GGally::ggpairs(
      aes(color = is_on_term),
      columns = df_enr %>% pull(method) %>% unique
    )
    p
    ggsave(file.path(pairwise_dir, glue::glue("pairwise_enrichment_{key}.pdf")), p)
  })

# benchmark plot
ggplot(df_runtime, aes(x = method, y = duration_seconds)) +
  geom_boxplot() +
  scale_y_log10(labels=function(x) {format(as.POSIXct(x), format = "%H:%M:%OS", tz = "UTC")}) +
  labs(x = "Method", y = "Runtime [hh:mm:ss]") +
  theme_minimal()
ggsave(file.path(outdir, "runtime.pdf"))
