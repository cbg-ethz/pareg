devtools::load_all()

library(tidyverse)
library(proxy)


# parameters
similarity_function_list <- list(
  jaccard = jaccard,
  overlap_coefficient = overlap_coefficient,
  semantic = function(term_list) {
    simplifyEnrichment::GO_similarity(term_list, measure = "Rel")
  }
)
pathway_db_list <- c("C2@CP:KEGG", "C5@GO:BP")
max_term_count <- 1000

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

    # subset to limit resource consumption
    term_list_list <- sample(term_list_list, min(length(term_list_list), max_term_count))

    imap(similarity_function_list, function(func, name) {
      if (name == "semantic") {
        # semantic similarity returns whole matrix right away
        go_list <- names(term_list_list)

        # check if term list is compatible (i.e., from Gene Ontology)
        res <- try(simplifyEnrichment::guess_ont(go_list), silent = TRUE)
        if (inherits(res, "try-error")) {
          # could not determine GO category, assume no GO terms
          mat <- NULL
        } else {
          mat <- func(go_list)
        }
      } else {
        # double "1-" because func is similarity, but we assume distance
        mat <- 1 - dist(
          x = term_list_list,
          method = function(x, y) {
            1 - func(x, y)
          },
          diag = TRUE, pairwise = TRUE
        ) %>%
          as.matrix()
      }

      if (!is.null(mat)) {
        mat[lower.tri(mat, diag = TRUE)] <- 0
        mat <- as(mat, "sparseMatrix")
      }

      mat
    }) %>%
      set_names(names(similarity_function_list))
  }) %>%
  set_names(df_grp %>% group_keys() %>% pull(pathway_db))

# store data
usethis::use_data(pathway_similarities, overwrite = TRUE)
