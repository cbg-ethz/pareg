library(tidyverse)
library(magrittr) # for %<>%
library(msigdbr)


# parameters
fname_out <- snakemake@output$fname

term_filter_params <- snakemake@params$term_filter

# overview
msigdbr_species()
msigdbr_collections()

# retrieve terms
df_terms <- msigdbr(category = "C2")

df_terms %>%
  head

df_terms %>%
  group_by(gs_name) %>%
  tally %>%
  arrange(desc(n))

# select terms of reasonable size
gs_selection <- df_terms %>%
  group_by(gs_name) %>%
  tally

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
  print("Sampling pathway")
  gs_selection %<>%
    sample_n(min(term_filter_params$sample_num, n()))
}

gs_selection %<>%
  pull(gs_name)

gs_selection %>%
  head

df_sel <- df_terms %>%
  filter(gs_name %in% gs_selection)

df_sel %>%
  head

df_sel %>%
  distinct(gs_name) %>%
  dim

# save result
df_sel %>%
  write_csv(fname_out)
