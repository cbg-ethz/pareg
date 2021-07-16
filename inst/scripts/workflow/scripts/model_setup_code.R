library(tidyverse)


# parameters
fname_study <- snakemake@input$fname_study
fname_terms <- snakemake@input$fname_terms
fname_term_sim <- snakemake@input$fname_term_sim

fname_out <- snakemake@output$fname

msc_category <- snakemake@params$params$category
msc_subcategory <- snakemake@params$params$subcategory

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
  filter(gs_cat == msc_category) %>%
  {
    if (msc_subcategory != "nan") {
      print("Filtering subcategory")
      filter(., gs_subcat == msc_subcategory)
    } else {
      print("Skipping subcategory filter")
      .
    }
  } %>%
  select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol) %>%
  distinct(.keep_all = TRUE)

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
