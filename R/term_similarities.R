jaccard <- function(x, y) {
  return(length(intersect(x, y)) / length(union(x, y)))
}

compute_term_similarities <- function(
  df_terms,
  similarity_function = jaccard,
  max_similarity = 1
) {
  term_list_list <- df_terms %>%
    select(term, gene) %>%
    { split(.$gene, .$term) }

  term_similarities <- proxy::simil(
    x = term_list_list,
    method = similarity_function,
    pairwise = TRUE
  ) %>%
    as.matrix()
  diag(term_similarities) <- max_similarity

  return(term_similarities)
}
