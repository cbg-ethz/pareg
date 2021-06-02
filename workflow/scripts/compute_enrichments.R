library(tidyverse)
library(pareg)


# parameters
fname_study <- snakemake@input$fname_study
fname_terms <- snakemake@input$fname_terms
fname_term_sim <- snakemake@input$fname_term_sim

fname_out <- snakemake@output$fname

category <- snakemake@params$params$category

# read data
study <- readRDS(fname_study)

df_terms <- read_csv(
  fname_terms,
  col_types = cols(
    gs_url = col_character(),
    gs_exact_source = col_character(),
    gs_geoid = col_character(),
    gs_pmid = col_character()
  )
) %>%
  filter(gs_cat == category) %>%
  select(gs_name, gene_symbol) %>%
  rename(name = gs_name, gene = gene_symbol) %>%
  distinct(.keep_all = TRUE)

df_terms %>% dim
df_terms %>% head

term_similarities <- read.csv(fname_term_sim, row.names = 1)

# prepare data
study_genes <- study$df %>%
  filter(pvalue <= 0.05) %>%
  pull(gene)

nonstudy_genes <- study$df %>%
  filter(pvalue > 0.05) %>%
  pull(gene)

all_genes <- df_terms %>%
  distinct(gene) %>%
  pull(gene)
# all_genes <- c(study_genes, nonstudy_genes)

# single-term methods
df_enr <- df_terms %>%
  group_by(name) %>%
  group_map(function(sub, key) {
    term <- key %>% pull(name)
    term_genes <- sub %>% pull(gene)

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
  df_terms,
  truncate_response = TRUE
) %>%
  mutate(method = "pareg", enrichment = abs(enrichment)) %>%
  rename(term = name)
df_pareg %>% arrange(desc(abs(enrichment))) %>% head

# pareg with network
term_names <- df_terms %>%
  distinct(name) %>%
  pull(name)
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
  df_terms,
  network_param = 0.9, term_network = term_similarities_sub,
  truncate_response = TRUE
) %>%
  mutate(method = "pareg_network", enrichment = abs(enrichment)) %>%
  rename(term = name)
df_pareg_network %>% arrange(desc(abs(enrichment))) %>% head

# pareg with network and cv
term_names <- df_terms %>%
  distinct(name) %>%
  pull(name)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix

df_pareg_network_cv <- pareg::pareg(
  data.frame(
    gene = all_genes,
    pvalue = ifelse(
      all_genes %in% study_genes,
      rbeta(length(all_genes), 0.1, 1),
      rbeta(length(all_genes), 1, 1)
    )
  ),
  df_terms,
  term_network = term_similarities_sub,
  truncate_response = TRUE,
  cv = TRUE
) %>%
  mutate(method = "pareg_network_cv", enrichment = abs(enrichment)) %>%
  rename(term = name)
df_pareg_network_cv %>% arrange(desc(abs(enrichment))) %>% head

# pareg with network and cv and BUM
term_names <- df_terms %>%
  distinct(name) %>%
  pull(name)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix

df_pareg_network_cv_bum <- pareg::pareg(
  data.frame(
    gene = all_genes,
    pvalue = ifelse(
      all_genes %in% study_genes,
      rbeta(length(all_genes), 0.1, 1),
      rbeta(length(all_genes), 1, 1)
    )
  ),
  df_terms,
  term_network = term_similarities_sub,
  truncate_response = TRUE,
  cv = TRUE,
  family = netReg::bum
) %>%
  mutate(method = "pareg_network_cv_bum", enrichment = abs(enrichment)) %>%
  rename(term = name)
df_pareg_network_cv_bum %>% arrange(desc(abs(enrichment))) %>% head

# MGSA
term_list <- df_terms %>%
  { split(.$gene, .$name) }

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
  df_pareg_network_cv,
  df_pareg_network_cv_bum,
  df_mgsa
)

df_enr_all %>%
  head

# save result
df_enr_all %>%
  write_csv(fname_out)
