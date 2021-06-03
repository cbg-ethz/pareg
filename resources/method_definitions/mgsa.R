# prepare environment
source(snakemake@params$setup_code_fname)

# run model
term_list <- df_terms %>%
  { split(.$gene, .$term) }

fit <- mgsa::mgsa(study_genes, term_list)

df <- fit@setsResults %>%
  rownames_to_column("term") %>%
  mutate(method = "MGSA") %>%
  rename(enrichment = estimate) %>%
  select(-inPopulation, -inStudySet, -std.error)

df %>%
  head

# save result
df %>%
  write_csv(snakemake@output$fname)
