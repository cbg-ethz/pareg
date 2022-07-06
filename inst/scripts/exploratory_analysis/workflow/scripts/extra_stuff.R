library(tidyverse)
library(topGO)
library(org.Hs.eg.db)
library(GO.db)

devtools::load_all()


# parameters
fname_de <- snakemake@input$fname_de
fname_terms <- snakemake@input$fname_terms
fname_sim <- snakemake@input$fname_sim

outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# load data
df_de <- read_tsv(fname_de) %>%
  dplyr::select(`SYMBOL...2`, PVAL) %>%
  dplyr::rename(gene = `SYMBOL...2`, pvalue = PVAL)
df_de %>%
  head()

df_terms <- read_csv(fname_terms)
df_terms %>%
  head()

term_similarities <- read.csv(fname_sim, row.names = 1, check.names = FALSE)
term_names <- df_terms %>%
  distinct(term) %>%
  pull(term)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix()

# run pareg without network regularization
fit_nonet <- pareg::pareg(
  df_de,
  df_terms,
  cv = FALSE,
  family = pareg::beta_phi_var,
  lasso_param = 1,
  network_param = 0,
  max_iteration = 10000,
  log_level = TRACE
)

df_enr <- fit_nonet %>%
  as.data.frame() %>%
  arrange(desc(abs(enrichment)))
head(df_enr)

# visualize pareg result
min_similarity <- 0.1
initial_term_count <- 50

enriched_terms <- df_enr %>%
  arrange(desc(abs(enrichment))) %>%
  head(initial_term_count) %>%
  pull(term)

plot(fit_nonet, term_subset = enriched_terms, min_similarity = min_similarity) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )

ggsave(
  file.path(outdir, glue("overview_pareg_nonet.pdf")),
  width = 15,
  height = 15
)

# run singular FET enrichment
gene_list <- df_de %>%
  dplyr::select(gene, pvalue) %>%
  deframe()
head(gene_list)

tmp <- AnnotationDbi::select(org.Hs.eg.db, df_de$gene, c("GO"), "SYMBOL")
terms_gene2GO <- split(tmp$GO, tmp$SYMBOL)
head(terms_gene2GO, 1)

gene_selector_func <- function(input_list) {
  return(input_list < 0.05)
}

topgo_obj <- new(
  "topGOdata",
  description = "Benchmark run",
  ontology = "BP",
  allGenes = gene_list,
  geneSel = gene_selector_func,
  annot = annFUN.gene2GO,
  gene2GO = terms_gene2GO
)

res_fet <- runTest(topgo_obj, algorithm = "classic", statistic = "fisher")

df_fet <- score(res_fet) %>%
  enframe() %>%
  dplyr::rename(id = name, pvalue = value) %>%
  arrange(pvalue)

df_fet$term <- AnnotationDbi::select(GO.db, df_fet$id, c("TERM"), "GOID") %>%
  mutate(
    TERM = TERM %>% str_to_upper() %>% str_replace_all(" ", "_")
  ) %>%
  pull(TERM)

head(df_fet)

# visualize FET results


ggsave(
  file.path(outdir, glue("overview_fet.pdf")),
  width = 10,
  height = 12
)
