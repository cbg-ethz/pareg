# prepare environment
source(snakemake@params$setup_code_fname)

# run model
df <- pareg::pareg(
  study$df %>% select(-in_study),
  df_terms,
  term_network = term_similarities_sub,
  truncate_response = TRUE,
  cv = TRUE
) %>%
  as.data.frame() %>%
  mutate(method = "pareg_network_cv", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
