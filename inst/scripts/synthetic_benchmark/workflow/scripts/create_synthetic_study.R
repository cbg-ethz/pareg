library(tidyverse)

devtools::load_all("../..")


# parameters
fname_terms <- snakemake@input$fname_terms
fname_sim <- snakemake@input$fname_sim

fname_rds <- snakemake@output$fname_rds
plotdir <- snakemake@output$plotdir
dir.create(plotdir, recursive = TRUE)

alpha <- snakemake@params$params$alpha # false positive rate
beta <- snakemake@params$params$beta # false negative rate
similarity_factor <- snakemake@params$params$similarityfactor
on_term_count <- snakemake@params$params$ontermcount
sig_gene_scaling <- snakemake@params$params$siggenescaling
model <- snakemake@params$params$model


# setup
set.seed(snakemake@wildcards$replicate)


# read data
df_terms <- read_csv(fname_terms)
sim_mat <- read.csv(fname_sim, row.names = 1, check.names = FALSE)

# select activated terms
# on_terms <- df_terms %>%
#   distinct(term) %>%
#   pull(term) %>%
#   sample(on_term_count, replace = TRUE)
on_terms <- pareg::similarity_sample(
  sim_mat,
  on_term_count,
  similarity_factor = similarity_factor
)

on_terms <- unique(on_terms)

# model DE p-values
if (model == "mgsa") {
  # extract activated genes from activated terms
  study_genes_orig <- df_terms %>%
    filter(term %in% on_terms) %>%
    distinct(gene) %>%
    pull(gene)

  # remove false negatives from activated genes
  fn_genes <- sample(seq_along(study_genes_orig), size = length(study_genes_orig) * beta)
  if (length(fn_genes) > 0) {
    study_genes <- study_genes_orig[-fn_genes]
  } else {
    study_genes <- study_genes_orig
  }

  # add false positives to activated genes
  other_genes_orig <- df_terms %>%
    filter(!(term %in% on_terms)) %>%
    distinct(gene) %>%
    pull(gene) %>%
    setdiff(study_genes_orig)
  stopifnot(length(intersect(study_genes_orig, other_genes_orig)) == 0)

  fp_genes <- sample(seq_along(other_genes_orig), size = length(study_genes_orig) * alpha)
  study_genes <- c(study_genes, other_genes_orig[fp_genes])

  # find inactive genes
  nonstudy_genes <- df_terms %>%
    distinct(gene) %>%
    filter(!(gene %in% study_genes)) %>%
    pull(gene)

  # compute how many pathways each gene is a member of
  if (sig_gene_scaling == "membercount") {
    member_count <- purrr::map_dfc(study_genes, function(current_gene) {
      df_terms %>%
        group_by(term) %>%
        summarise("{current_gene}" := current_gene %in% gene) %>%
        select({{ current_gene }})
    }) %>%
      summarize(across(everything(), sum)) %>%
      t %>%
      as.data.frame %>%
      pull(V1)

    # +2 to have reasonable parameter with member_count=1
    rbeta_shape1_param <- 2 ^ -(2 + member_count)
  } else {
    rbeta_shape1_param <- as.numeric(sig_gene_scaling)
  }

  # compute "DE" p-values for active and inactive genes
  study_pvalues <- rbeta(length(study_genes), rbeta_shape1_param, 1) # peak at 0
  nonstudy_pvalues <- rbeta(length(nonstudy_genes), 1, 1) # uniform

  # aggregate results
  df <- data.frame(
    gene = c(study_genes, nonstudy_genes),
    pvalue = c(study_pvalues, nonstudy_pvalues),
    in_study = c(
      rep(TRUE, length(study_genes)),
      rep(FALSE, length(nonstudy_genes))
    ),
    in_orig_study = c(
      study_genes %in% study_genes_orig,
      nonstudy_genes %in% study_genes_orig
    )
  )

  # sanity checks
  stopifnot(all(sort(unique(c(study_genes, nonstudy_genes))) == sort(unique(df_terms$gene))))
} else if (model == "linear") {
  # generate ground truth coefficient vector
  terms <- df_terms %>%
    distinct(term) %>%
    pull(term)

  term_vector <- rep(0, times = length(terms))
  names(term_vector) <- terms

  term_vector[on_terms] <- -1

  # generate term-membership matrix
  mat_X <- df_terms %>%
    group_by(.data$term) %>%
    mutate(member = TRUE) %>%
    ungroup() %>%
    mutate_at(vars(.data$gene), as.character) %>%
    pivot_wider(
      names_from = .data$term,
      values_from = .data$member,
      values_fill = FALSE
    ) %>%
    select(
      # use `any_of` to handle case where NA does not exist
      # (happens when a gene appears in no term)
      -any_of("NA")
    ) %>%
    mutate_at(
      vars(.data$gene),
      factor # gene is character if select statement is executed
    )

  genes <- mat_X$gene
  mat_X <- mat_X %>%
    select(-gene)

  mat_X <- mat_X %>%
    mutate(across(everything(), as.numeric))

  # visualize raw matrix
  tmp <- as.matrix(mat_X)
  rownames(tmp) <- genes

  ht_raw <- ComplexHeatmap::Heatmap(
    tmp,
    col = c("white", "black"),
    column_title = "Raw matrix",
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  # introduce noise
  for (term in on_terms) {
    neg_matches <- which(mat_X[[term]] == 0)
    pos_matches <- which(mat_X[[term]] == 1)

    # false positives
    fp_count <- as.integer(alpha * length(neg_matches))
    fp_ind <- sample(neg_matches, fp_count, replace = FALSE)
    mat_X[fp_ind, term] <- 1

    # false negatives
    fn_count <- as.integer(beta * length(pos_matches))
    fn_ind <- sample(pos_matches, fn_count, replace = FALSE)
    mat_X[fn_ind, term] <- 0
  }

  # visualize noisy matrix
  tmp <- as.matrix(mat_X)
  rownames(tmp) <- genes

  ht_noisy <- ComplexHeatmap::Heatmap(
    tmp,
    col = c("white", "black"),
    column_title = "Noisy matrix",
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )

  png(
    file.path(plotdir, "member_matrix.png"),
    width = 40,
    height = 20,
    units = "in",
    res = 300
  )
  ComplexHeatmap::draw(
    ht_raw + ht_noisy,
    row_title = "Genes",
    column_title = "Terms"
  )
  dev.off()

  # simulate p-values
  term_vector <- term_vector[colnames(mat_X)]

  eta <- as.matrix(mat_X) %*% term_vector

  linkinv <- pareg::beta_phi_var()$linkinv
  mu <- linkinv(eta)$numpy()
  phi <- 1
  p_values <- rbeta(length(mu), mu * phi, (1 - mu) * phi)
  p_values <- transform_y(p_values)
  names(p_values) <- genes

  # aggregate results
  study_genes <- df_terms %>%
    filter(term %in% on_terms) %>%
    pull(gene)

  df <- data.frame(
    gene = names(p_values),
    pvalue = unname(p_values),
    in_study = names(p_values) %in% study_genes,
    in_orig_study = names(p_values) %in% study_genes # degenerate case
  )
} else {
  stop("Invalid model")
}

# summarize results
table(df$pvalue <= 0.05, df$in_study, dnn = c("sig. p-value", "in study"))
table(
  df$in_study,
  df$in_orig_study,
  dnn = c("gene in generated study", "gene in original study")
)

data.frame(
  min_pvalue = min(df$pvalue),
  median_pvalue = median(df$pvalue),
  max_pvalue = max(df$pvalue),
  summary = paste(as.vector(summary(df$pvalue)), collapse = "__")
) %>%
  write_csv(file.path(plotdir, "pvalue_stats.csv"))

# overview plots
ggplot(df, aes(x = pvalue)) +
  geom_histogram() +
  facet_wrap(~ in_study) +
  theme_minimal()
ggsave(
  file.path(plotdir, "study_pvalue_histogram.pdf"),
  width = 16
)

ggplot(df, aes(x = pvalue)) +
  geom_histogram() +
  facet_wrap(~ in_orig_study) +
  theme_minimal()
ggsave(
  file.path(plotdir, "orig_study_pvalue_histogram.pdf"),
  width = 16
)

# save result
saveRDS(
  list(on_terms = on_terms, df = df),
  file = fname_rds
)
