library(tidyverse)


# parameters
fname_terms <- snakemake@input$fname_terms
fname_rds <- snakemake@output$fname_rds

alpha <- snakemake@params$params$alpha # false positive rate
beta <- snakemake@params$params$beta # false negative rate

on_term_count <- snakemake@params$on_term_count

# read data
df_terms <- read_csv(fname_terms, col_types = cols(gs_url = col_character()))

# create synthetic studies
on_terms <- df_terms %>%
  distinct(gs_name) %>%
  pull(gs_name) %>%
  sample(on_term_count)

study_genes_orig <- df_terms %>%
  filter(gs_name %in% on_terms) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

fn_genes <- sample(seq_along(study_genes_orig), size = length(study_genes_orig) * beta)
study_genes <- study_genes_orig[-fn_genes]

other_genes_orig <- df_terms %>%
  filter(!(gs_name %in% on_terms)) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol) %>%
  setdiff(study_genes_orig)
stopifnot(length(intersect(study_genes_orig, other_genes_orig)) == 0)

fp_genes <- sample(seq_along(other_genes_orig), size = length(other_genes_orig) * alpha)
study_genes <- c(study_genes, other_genes_orig[fp_genes])

nonstudy_genes <- df_terms %>%
  distinct(gene_symbol) %>%
  filter(!(gene_symbol %in% study_genes)) %>%
  pull(gene_symbol)

# save result
saveRDS(
  list(on_terms = on_terms, study_genes = study_genes, nonstudy_genes = nonstudy_genes),
  file = fname_rds
)
