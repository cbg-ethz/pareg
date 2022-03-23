#' Collection of pathway similarity matrices.
#'
#' Contains matrices for various pathway databases and similarity measures.
#' Note that the matrices are sparse and upper triangular.
#' To use with \code{pareg}, transform them using \code{pareg::as_dense_sim}.
#'
#' @format A list of lists of matrices.
#' \describe{
#'   \itemize{
#'     \item Pathway database 1
#'     \itemize{
#'       \item Similarity measure 1
#'       \item Similarity measure 2
#'       \item ...
#'     }
#'     \item Pathway database 2
#'     \item ...
#'   }
#' }
"pathway_similarities"
