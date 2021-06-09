#' @title Jaccard similarity.
#'
#' @description Compute Jaccard similarity between two sets.
#'
#' @param x First set.
#' @param y Second set.
#'
#' @return Jaccard similarity between set x and y.
#'
#' @examples
#' jaccard(c(1, 2, 3), c(2, 3, 4))
jaccard <- function(x, y) {
  return(length(intersect(x, y)) / length(union(x, y)))
}


#' @title Term similarity computation.
#'
#' @description Generate similarity matrix for input terms.
#'
#' @export
#'
#' @param df_terms Dataframe storing pathway database.
#' @param similarity_function Function to compute similarity between two sets.
#' @param max_similarity Value to fill diagonal with.
#'
#' @return Symmetric matrix of similarity scores.
#'
#' @examples
#' df_terms <- data.frame(
#'   term = c("A", "A", "B", "B", "B", "C", "C", "C"),
#'   gene = c("a", "b", "a", "b", "c", "a", "c", "d")
#' )
#' compute_term_similarities(df_terms)
#' @importFrom rlang .data
#' @importFrom dplyr select
#' @importFrom proxy simil
compute_term_similarities <- function(
  df_terms,
  similarity_function = jaccard,
  max_similarity = 1
) {
  term_list_list <- df_terms %>%
    select(term, gene) %>% {
        split(.$gene, .$term)
    }

  term_similarities <- proxy::simil(
    x = term_list_list,
    method = similarity_function,
    pairwise = TRUE
  ) %>%
    as.matrix()
  diag(term_similarities) <- max_similarity

  return(term_similarities)
}
