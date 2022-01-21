library(tidyverse)
library(magrittr) # for %<>%
library(msigdbr)


# parameters
fname_out <- snakemake@output$fname

term_filter_params <- snakemake@params$term_filter
termsource <- snakemake@params$params$termsource

parts <- strsplit(termsource, ";")[[1]]

if (parts[[1]] == "msigdb") {
  category <- parts[[2]]
  subcategory <- parts[[3]]

  # overview (http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
  msigdbr_species()
  msigdbr_collections() %>%
    arrange(desc(num_genesets))
  msigdbr_collections() %>%
    group_by(gs_cat) %>%
    summarize(num_genesets = sum(num_genesets)) %>%
    arrange(desc(num_genesets))

  # retrieve terms
  df_terms <- msigdbr(species = "Homo sapiens")

  # subset terms
  df_terms <- df_terms %>%
    filter(gs_cat == category) %>%
    {
      if (subcategory != "None") {
        print("Filtering subcategory")
        filter(., gs_subcat == subcategory)
      } else {
        print("Skipping subcategory filter")
        .
      }
    } %>%
    select(gs_name, gene_symbol) %>%
    rename(term = gs_name, gene = gene_symbol) %>%
    distinct(.keep_all = TRUE)
} else {
  stop(paste("Unknown term database source:", parts[[1]]))
}

# select terms of reasonable size
df_term_tally <- df_terms %>%
  group_by(term) %>%
  tally() %>%
  arrange(desc(n))

if (!is.null(term_filter_params$min_size)) {
  print("Filtering pathway by min_size")
  df_term_tally %<>%
    filter(term_filter_params$min_size < n)
}
if (!is.null(term_filter_params$max_size)) {
  print("Filtering pathway by max_size")
  df_term_tally %<>%
    filter(n < term_filter_params$max_size)
}
if (!is.null(term_filter_params$sample_num)) {
  print("Sampling pathways")
  df_term_tally %<>%
    group_by(gs_cat) %>%
    sample_n(min(term_filter_params$sample_num, n())) %>%
    ungroup()
}

term_selection <- df_term_tally %>%
  pull(term)
term_selection %>%
  head()

df_sub <- df_terms %>%
  filter(term %in% term_selection)

# data overview
df_sub %>%
  head()

print("Distinct term count:")
df_sub %>%
  distinct(term) %>%
  dim()

print("Genes per term:")
df_sub %>%
  group_by(term) %>%
  tally() %>%
  arrange(desc(n))

# save result
df_sub %>%
  write_csv(fname_out)
