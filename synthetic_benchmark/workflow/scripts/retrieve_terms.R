library(tidyverse)
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
  tally %>%
  filter(term_filter_params$min_size < n & n < term_filter_params$max_size) %>%
  sample_n(min(term_filter_params$sample_num, n())) %>%
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
