library(tidyverse)

devtools::load_all("../..")


# parameters
fname_terms <- snakemake@input$fname_terms
fname_sim <- snakemake@input$fname_sim
fname_rds <- snakemake@output$fname_rds

category <- snakemake@params$params$category
subcategory <- snakemake@params$params$subcategory
alpha <- snakemake@params$params$alpha # false positive rate
beta <- snakemake@params$params$beta # false negative rate
similarity_factor <- snakemake@params$params$similarityfactor
on_term_count <- snakemake@params$params$ontermcount
sig_gene_scaling <- snakemake@params$params$siggenescaling

# read data
df_terms <- read_csv(
  fname_terms,
  col_types = cols(
    gs_url = col_character(),
    gs_exact_source = col_character(),
    gs_geoid = col_character(),
    gs_pmid = col_character()
  )
) %>%
  filter(gs_cat == category) %>%
  {
    if (subcategory != "None") {
      print("Filtering subcategory")
      filter(., gs_subcat == subcategory)
    } else {
      print("Skipping subcategory filter")
      .
    }
  }

sim_mat <- read.csv(fname_sim, row.names = 1)

# select activated terms
# on_terms <- df_terms %>%
#   distinct(gs_name) %>%
#   pull(gs_name) %>%
#   sample(on_term_count, replace = TRUE)
on_terms <- pareg::similarity_sample(
  sim_mat,
  on_term_count,
  similarity_factor = similarity_factor
)

on_terms <- unique(on_terms)

# extract activated genes from activated terms
study_genes_orig <- df_terms %>%
  filter(gs_name %in% on_terms) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

# remove false negatives from activated genes
fn_genes <- sample(seq_along(study_genes_orig), size = length(study_genes_orig) * beta)
if (length(fn_genes) > 0) {
  study_genes <- study_genes_orig[-fn_genes]
} else {
  study_genes <- study_genes_orig
}

# add false positives to activated genes
other_genes_orig <- df_terms %>%
  filter(!(gs_name %in% on_terms)) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol) %>%
  setdiff(study_genes_orig)
stopifnot(length(intersect(study_genes_orig, other_genes_orig)) == 0)

fp_genes <- sample(seq_along(other_genes_orig), size = length(other_genes_orig) * alpha)
study_genes <- c(study_genes, other_genes_orig[fp_genes])

# find inactive genes
nonstudy_genes <- df_terms %>%
  distinct(gene_symbol) %>%
  filter(!(gene_symbol %in% study_genes)) %>%
  pull(gene_symbol)

# compute how many pathways each gene is a member of
if (sig_gene_scaling == "membercount") {
  member_count <- purrr::map_dfc(study_genes, function(gene) {
    df_terms %>%
      group_by(gs_name) %>%
      summarise("{gene}" := gene %in% gene_symbol) %>%
      select({{ gene }})
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

df <- data.frame(
  gene = c(study_genes, nonstudy_genes),
  pvalue = c(study_pvalues, nonstudy_pvalues),
  in_study = c(
    rep(TRUE, length(study_genes)),
    rep(FALSE, length(nonstudy_genes))
  )
)
table(df$pvalue <= 0.05, df$in_study, dnn = c("sig. p-value", "in study"))

# save result
saveRDS(
  list(on_terms = on_terms, df = df),
  file = fname_rds
)
