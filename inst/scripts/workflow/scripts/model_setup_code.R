library(tidyverse)


# parameters
fname_study <- snakemake@input$fname_study
fname_terms <- snakemake@input$fname_terms
fname_term_sim <- snakemake@input$fname_term_sim

fname_out <- snakemake@output$fname

# read data
study <- readRDS(fname_study)

df_terms <- read_csv(fname_terms)

df_terms %>% dim()
df_terms %>% head()

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

term_names <- df_terms %>%
  distinct(term) %>%
  pull(term)
term_similarities_sub <- term_similarities[term_names, term_names] %>%
  as.matrix()
