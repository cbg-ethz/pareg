# parameters
category <- snakemake@params$params$category
subcategory <- snakemake@params$params$subcategory

# topGO only works with Gene Ontology
if (category != "C5") {
  print("Skipping topGO method due to invalid category")
  file.create(snakemake@output$fname) # create empty csv which is skipped by aggregation
  quit(save = "no")
}

subcategory <- substring(subcategory, 4) # remove "GO:" prefix

# prepare environment
source(snakemake@params$setup_code_fname)

library(topGO)

# run model
gene_list <- study$df %>%
  dplyr::select(gene, pvalue) %>%
  deframe()

name_map <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = subcategory) %>%
  dplyr::select(gs_name, gs_exact_source) %>%
  distinct(gs_exact_source, .keep_all = TRUE)

df_terms_goid <- name_map %>%
  right_join(df_terms, by = c("gs_name" = "term"))
df_terms_goid %>%
  head

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
  gene2GO = split(df_terms_goid$gs_exact_source, df_terms_goid$gene)
)

res_elim <- runTest(topgo_obj, algorithm = "elim", statistic = "ks")

df <- GenTable(topgo_obj, enrichment = res_elim) %>%
  as_tibble %>%
  dplyr::select(GO.ID, enrichment) %>%
  left_join(name_map, by = c("GO.ID" = "gs_exact_source")) %>%
  dplyr::rename(term = gs_name) %>%
  dplyr::select(term, enrichment) %>%
  mutate(method = "topGO", enrichment = -log10(as.numeric(enrichment)))

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
