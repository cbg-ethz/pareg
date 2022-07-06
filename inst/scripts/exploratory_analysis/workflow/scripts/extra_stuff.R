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

df_fet <- df_fet %>%
  filter(term %in% df_terms$term)

head(df_fet)

# compute term similarities
term_list_list <- AnnotationDbi::select(
    org.Hs.eg.db,
    df_de$gene,
    c("GO"),
    "SYMBOL"
  ) %>%
  inner_join(
    df_fet,
    by = c("GO" = "id")
  ) %>%
  pipe_split("term", "SYMBOL")
length(term_list_list)

term_similarities <- 1 - proxy::dist(
  x = term_list_list,
  method = function(x, y) {
    1 - jaccard(x, y)
  },
  diag = TRUE, pairwise = TRUE
) %>%
  as.matrix()

df_fet <- df_fet %>%
  filter(df_fet$term %in% colnames(term_similarities))

# visualize FET results
term_sizes <- df_terms %>%
  group_by(.data$term) %>%
  summarize(size = n())

min_similarity <- 0.1
initial_term_count <- 50

enriched_terms <- df_fet %>%
  arrange(pvalue) %>%
  head(initial_term_count) %>%
  pull(term)

term_network_tmp <- term_similarities[enriched_terms, enriched_terms]
term_network_tmp[term_network_tmp < min_similarity] <- 0

term_graph <- as_tbl_graph(graph_from_adjacency_matrix(
  term_network_tmp,
  weighted = TRUE
)) %>%
  activate("nodes") %>%
  mutate(
    pvalue = data.frame(term = .data$name) %>%
      left_join(df_fet, by = "term") %>%
      pull(.data$pvalue),
    enrichment = -log10(pvalue),
    term_size = data.frame(term = .data$name) %>%
      left_join(term_sizes, by = "term") %>%
      pull(.data$size),
  ) %>%
  activate("edges") %>%
  mutate(
    term_similarity = .data$weight
  )

term_graph %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(alpha = .data$term_similarity)) +
  geom_node_point(
    aes(size = .data$term_size, color = .data$enrichment)
  ) +
  scale_size(range = c(2, 10), name = "Term size") +
  scale_color_gradient2(
    low = "red",
    mid = "grey",
    high = "blue",
    midpoint = 0,
    na.value = "black",
    name = "Enrichment"
  ) +
  scale_edge_alpha(name = "Term similarity") +
  geom_text_repel(
    aes(label = .data$name, x = .data$x, y = .data$y),
    color = "black",
    bg.color = "white"
  ) +
  coord_fixed() +
  theme(
    panel.background = element_rect(fill = "white")
  )

ggsave(
  file.path(outdir, glue("overview_fet.pdf")),
  width = 15,
  height = 15
)
