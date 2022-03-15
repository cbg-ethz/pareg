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
#' @param cv Estimate best regularization parameters using cross-validation.
#' @param family Distribution family of response.
#' @param response_column_name Which column of model dataframe
#' to use as response.
#' @param max_iterations How many iterations to maximally run optimizer for.
#' @param ... Further arguments to pass to `(cv.)edgenet`.
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
#' pareg(df_genes, df_terms, max_iterations = 10)
#' @importFrom dplyr select ends_with all_of
#' @importFrom glue glue_collapse
#' @importFrom magrittr %>%
pareg <- function(
  df_genes,
  df_terms,
  lasso_param = NA_real_,
  network_param = NA_real_,
  term_network = NULL,
  cv = FALSE,
  family = beta,
  response_column_name = "pvalue",
  max_iterations = 1e5,
  ...
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
    select(response_column_name) %>%
    as.matrix()

  if (!is.null(term_network)) {
    ordered_terms <- vapply(
      strsplit(covariates, ".", fixed = TRUE),
      function(x) {
        glue_collapse(x[seq_len(length(x) - 1)], sep = ".")
      },
      FUN.VALUE = character(1)
    )

    term_diff <- setdiff(ordered_terms, rownames(term_network))
    if (length(term_diff) > 0) {
      msg <- paste(
        "The following covariates do not appear in term network:",
        glue_collapse(term_diff, sep = ", ")
      )
      stop(msg)
    }

    term_network <- term_network[ordered_terms, ordered_terms]
  }

  if (family()$family %in% c("beta", "beta_phi_lm", "beta_phi_var")) {
    # transform response from [0, 1] to (0, 1) if needed
    eps <- .Machine$double.eps * 1e9
    if (min(Y) < eps || max(Y) > 1 - eps) {
      Y <- transform_y(Y)
    }
  }

  # fit model
  if (cv) {
    fit_func <- cv.edgenet
  } else {
    fit_func <- edgenet

    if (is.na(lasso_param)) {
      lasso_param <- 0
    }
    if (is.na(network_param)) {
      network_param <- 0
    }
  }

  fit <- fit_func(
    X,
    Y,
    G.X = term_network,
    lambda = lasso_param,
    psigx = network_param,
    psigy = 0,
    family = family,
    maxit = max_iterations,
    ...
  )

  # update parameters to cross-validation estimates if needed
  if (cv) {
    lasso_param <- fit$lambda
    network_param <- fit$psigx
  }

  # return structured object
  return(structure(
    list(
      obj = fit,
      df_genes = df_genes,
      df_terms = df_terms,
      term_network = term_network,
      covariates = covariates,
      X = X,
      Y = Y,
      response_column_name = response_column_name,
      params = list(
        lasso_param = lasso_param,
        network_param = network_param,
        cv = cv
      )
    ),
    class = "pareg"
  ))
}
