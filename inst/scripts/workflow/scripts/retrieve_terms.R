library(tidyverse)
library(magrittr) # for %<>%
library(msigdbr)


# parameters
fname_out <- snakemake@output$fname

term_filter_params <- snakemake@params$term_filter

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

df_terms %>%
  head()

df_terms %>%
  group_by(gs_name) %>%
  tally() %>%
  arrange(desc(n))

# select terms of reasonable size
gs_selection <- df_terms %>%
  group_by(gs_cat, gs_name) %>%
  tally() %>%
  arrange(desc(n))

if (!is.null(term_filter_params$min_size)) {
  print("Filtering pathway by min_size")
  gs_selection %<>%
    filter(term_filter_params$min_size < n)
}
if (!is.null(term_filter_params$max_size)) {
  print("Filtering pathway by max_size")
  gs_selection %<>%
    filter(n < term_filter_params$max_size)
}
if (!is.null(term_filter_params$sample_num)) {
  print("Sampling pathways")
  gs_selection %<>%
    group_by(gs_cat) %>%
    sample_n(min(term_filter_params$sample_num, n())) %>%
    ungroup()
}

gs_selection %<>%
  pull(gs_name)

gs_selection %>%
  head()

df_sel <- df_terms %>%
  filter(gs_name %in% gs_selection)

df_sel %>%
  head()

df_sel %>%
  distinct(gs_name) %>%
  dim()

df_sel %>%
  distinct(gs_cat)

df_sel %>%
  distinct(gs_subcat)

# save result
df_sel %>%
  write_csv(fname_out)
