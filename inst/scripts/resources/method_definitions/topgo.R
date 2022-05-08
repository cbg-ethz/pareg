# parameters
termsource <- snakemake@params$params$termsource

parts <- strsplit(termsource, "@")[[1]]
db <- parts[[1]]
category <- parts[[2]]

# topGO only works with Gene Ontology
if (db != "msigdb" || category != "C5") {
  print("Skipping topGO method due to invalid database or category")
  file.create(snakemake@output$fname) # create empty csv which is skipped by aggregation
  quit(save = "no")
}

subcategory <- substring(parts[[3]], 4) # remove "GO:" prefix

# prepare environment
source(snakemake@params$setup_code_fname)

library(topGO)

# run model
gene_list <- study$df %>%
  dplyr::select(gene, pvalue) %>%
  deframe()

tmp <- df_terms %>%
  mutate(term = term %>% str_to_upper %>% str_replace("_", ":"))
terms_gene2GO <- split(tmp$term, tmp$gene)
terms_gene2GO %>%
  head(1)

gene_selector_func <- function(input_list) {
  return(input_list < 0.05)
}
topgo_obj <- new(
  "topGOdata",
  description = "Benchmark run",
  ontology = subcategory,
  allGenes = gene_list,
  geneSel = gene_selector_func,
  annot = annFUN.gene2GO,
  gene2GO = terms_gene2GO
)

res_elim <- runTest(topgo_obj, algorithm = "elim", statistic = "ks")

df <- GenTable(topgo_obj, enrichment = res_elim) %>%
  as_tibble %>%
  dplyr::select(GO.ID, enrichment) %>%
  dplyr::rename(term = GO.ID) %>%
  mutate(
    term = term %>% str_to_lower %>% str_replace(":", "_"),
    enrichment = as.numeric(if_else(
      enrichment == "< 1e-30",
      "1e-30",
      enrichment
    ))
  ) %>%
  mutate(method = "topgo", enrichment = -log10(enrichment))

# there might be some NA terms because topGO uses its own ontology
df %>%
  filter(is.na(term))
df <- df %>%
  drop_na(term)

df %>%
  arrange(desc(enrichment)) %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
