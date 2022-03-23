#' Collection of pathway similarity matrices.
#'
#' Contains matrices for various pathway databases and similarity measures.
#' Note that the matrices are sparse and upper triangular.
#' To use with \code{pareg}, transform them using \code{pareg::as_dense_sim}.
#'
#' @format A list of lists of matrices.
#' * Pathway database 1
#'   * Similarity measure 1
#'   * Similarity measure 2
#'   * ...
#' * Pathway database 2
#' * ...
"pathway_similarities"
