# prepare environment
source(snakemake@params$setup_code_fname)

devtools::load_all("../..")

# run model
df <- pareg::pareg(
  study$df %>% select(-in_study),
  df_terms,
  network_param = 0.9, term_network = term_similarities_sub,
  truncate_response = TRUE
) %>%
  as.data.frame() %>%
  mutate(method = "pareg_network", enrichment = abs(enrichment))

df %>%
  arrange(desc(abs(enrichment))) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
