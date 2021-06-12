library(tidyverse)
library(netReg)


#' @title Pathway enrichment using a regularized regression approach.
#'
#' @description Run model to compute pathway enrichments.
#'  Can model inter-pathway relations, cross-validation and much more.
#'
#' @export
#'
#' @param df_genes Dataframe storing gene names and DE p-values.
#' @param df_terms Dataframe storing pathway database.
#' @param lasso_param Lasso regularization parameter.
#' @param network_param Network regularization parameter.
#' @param term_network Term similarity network as adjacency matrix.
#' @param truncate_response Make sure response is in (0, 1).
#' @param cv Estimate best regularization parameters using cross-validation.
#' @param family Distribution family of response.
#'
#' @return An object of class \code{pareg}.
#'
#' @examples
#' df_genes <- data.frame(
#'   gene = paste("g", 1:20, sep = ""),
#'   pvalue = c(
#'     rbeta(10, .1, 1),
#'     rbeta(10, 1, 1)
#'   )
#' )
#' df_terms <- rbind(
#'   data.frame(
#'     term = "foo",
#'     gene = paste("g", 1:10, sep = "")
#'   ),
#'   data.frame(
#'     term = "bar",
#'     gene = paste("g", 11:20, sep = "")
#'   )
#' )
#' pareg(df_genes, df_terms)
#' @import tidyverse
#' @importFrom netReg edgenet cv.edgenet beta
#' @importFrom glue glue_collapse
pareg <- function(
  df_genes, df_terms,
  lasso_param = NA_real_, network_param = NA_real_,
  term_network = NULL, truncate_response = FALSE,
  cv = FALSE,
  family = netReg::beta
) {
  # generate design matrix
  df_model <- create_model_df(df_genes, df_terms)

  # setup data
  covariates <- df_model %>%
    select(ends_with(".member")) %>%
    names()

  X <- df_model %>%
    select(all_of(covariates)) %>%
    as.matrix()
  Y <- df_model %>%
    select("pvalue") %>%
    as.matrix()

  if (!is.null(term_network)) {
    ordered_terms <- sapply(
      strsplit(covariates, ".", fixed = TRUE),
      function(x) {
        glue_collapse(x[1:(length(x)-1)], sep = ".")
      }
    )

    term_diff <- setdiff(ordered_terms, rownames(term_network))
    if (length(term_diff) > 0) {
      stop(paste("The following covariates do not appear in term network:", glue_collapse(term_diff, sep = ", ")))
    }

    term_network <- term_network[ordered_terms, ordered_terms]
  }

  # truncate response if requested
  if (truncate_response) {
    eps <- .Machine$double.eps * 1e9 # see ?mgcv::betar
    Y[Y > 1 - eps] <- 1 - eps
    Y[Y < eps] <- eps
  }

  # fit model
  if (cv) {
    fit_func <- netReg::cv.edgenet
  } else {
    fit_func <- netReg::edgenet

    if (is.na(lasso_param)) {
      lasso_param <- 0
    }
    if (is.na(network_param)) {
      network_param <- 0
    }
  }

  fit <- fit_func(
    X, Y,
    G.X = term_network,
    lambda = lasso_param, psigx = network_param, psigy = 0,
    family = family
  )

  # return structured object
  return(structure(list(
    obj = fit,
    term_network = term_network,
    df_terms = df_terms,
    covariates = covariates,
    params = list(
      lasso_param = lasso_param,
      network_param = network_param,
      cv = cv,
      truncate_response = truncate_response
    )
  ), class = "pareg"))
}
