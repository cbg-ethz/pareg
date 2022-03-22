devtools::load_all()

library(tidyverse)
library(proxy)


# parameters
similarity_function_list <- list(jaccard = jaccard, overlap_coefficient = overlap_coefficient)
pathway_db_list <- c("C2@CP:KEGG", "C2@CP:REACTOME")

# retrieve pathways
df_terms <- msigdbr::msigdbr(species = "Homo sapiens") %>%
  mutate(
    pathway_db = str_c(gs_cat, gs_subcat, sep = "@")
  )

# compute pathways pathway similarities
df_grp <- df_terms %>%
  filter(pathway_db %in% pathway_db_list) %>%
  select(pathway_db, gs_name, ensembl_gene) %>%
  rename(term = gs_name, gene = ensembl_gene) %>%
  group_by(pathway_db)

pathway_similarities <- df_grp %>%
  group_split() %>%
  map(function(df_sub) {
    term_list_list <- df_sub %>%
      select(term, gene) %>%
      pipe_split("term", "gene")

    map(similarity_function_list, function(func) {
      # double "1-" because func is similarity, bt we assume distance
      1 - dist(
        x = term_list_list,
        method = function(x, y) {
          1 - func(x, y)
        },
        diag = TRUE, pairwise = TRUE
      ) %>%
        as.matrix()
    }) %>%
      set_names(names(similarity_function_list))
  }) %>%
  set_names(df_grp %>% group_keys() %>% pull(pathway_db))

# store data
usethis::use_data(pathway_similarities, overwrite = TRUE)
