# prepare environment
source(snakemake@params$setup_code_fname)

# run model
term_list <- df_terms %>%
  pull(term) %>%
  unique()

df <- data.frame(
  method = "random",
  term = term_list,
  enrichment = -log10(runif(length(term_list), min = 0, max = 1))
)

df %>%
  arrange(desc(enrichment)) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
