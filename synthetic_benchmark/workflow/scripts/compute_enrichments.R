library(tidyverse)
library(pareg)


# parameters
fname_study <- snakemake@input$fname_study
fname_terms <- snakemake@input$fname_terms
fname_term_sim <- snakemake@input$fname_term_sim

fname_out <- snakemake@output$fname

# read data
study <- readRDS(fname_study)
study_genes <- study$study_genes
nonstudy_genes <- study$nonstudy_genes

df_terms <- read_csv(fname_terms, col_types = cols(gs_url = col_character()))

term_similarities <- read.csv(fname_term_sim, row.names = 1)

all_genes <- df_terms %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

#all_genes <- c(study_genes, nonstudy_genes)

# single-term methods
df_enr <- df_terms %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  group_map(function(sub, key) {
    term <- key %>% pull(gs_name)
    term_genes <- sub %>% pull(gene_symbol)

    nonterm_genes <- setdiff(all_genes, term_genes)

    # fisher's exact test
    cont_table <- matrix(c(
      length(intersect(term_genes, study_genes)), length(intersect(term_genes, nonstudy_genes)),
      length(intersect(nonterm_genes, study_genes)), length(intersect(nonterm_genes, nonstudy_genes))
    ), 2, 2)

    fet_p_value <- fisher.test(cont_table)$p.value
    # print(paste(term, p_value))

    # store results
    bind_rows(
      data.frame(method = "FET", term = term, enrichment = -log10(fet_p_value))
    )
  }) %>%
  { do.call(rbind.data.frame, .) }

# pareg no regularization
df_pareg <- pareg::pareg(
  data.frame(
    gene = all_genes,
    pvalue = ifelse(
      all_genes %in% study_genes,
      rbeta(length(all_genes), 0.1, 1),
      rbeta(length(all_genes), 1, 1)
    )
  ),
  df_terms %>%
    select(gs_name, gene_symbol) %>%
    rename(name = gs_name, gene = gene_symbol),
  truncate_response = TRUE
) %>%
  mutate(method = "pareg", enrichment = abs(enrichment)) %>%
  rename(term = name)
df_pareg %>% arrange(desc(abs(enrichment))) %>% head

# pareg with network
term_names <- df_terms %>%
  distinct(gs_name) %>%
  pull(gs_name)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix

df_pareg_network <- pareg::pareg(
  data.frame(
    gene = all_genes,
    pvalue = ifelse(
      all_genes %in% study_genes,
      rbeta(length(all_genes), 0.1, 1),
      rbeta(length(all_genes), 1, 1)
    )
  ),
  df_terms %>%
    select(gs_name, gene_symbol) %>%
    rename(name = gs_name, gene = gene_symbol),
  network_param = 0.9, term_network = term_similarities_sub,
  truncate_response = TRUE
) %>%
  mutate(method = "pareg_network", enrichment = abs(enrichment)) %>%
  rename(term = name)
df_pareg_network %>% arrange(desc(abs(enrichment))) %>% head

# MGSA
term_list <- df_terms %>%
  { split(.$gene_symbol, .$gs_name) }

fit <- mgsa::mgsa(study_genes, term_list)

df_mgsa <- fit@setsResults %>%
  rownames_to_column("term") %>%
  mutate(method = "MGSA") %>%
  rename(enrichment = estimate) %>%
  select(-inPopulation, -inStudySet, -std.error)
df_mgsa %>% head

# combine results
df_enr_all <- bind_rows(
  df_enr,
  df_pareg,
  df_pareg_network,
  df_mgsa
)

df_enr_all %>%
  head

# save result
df_enr_all %>%
  write_csv(fname_out)
