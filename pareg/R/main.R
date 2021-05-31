library(tidyverse)
library(netReg)


pareg <- function(
  df_genes, df_terms,
  lasso_param = NA_real_, network_param = NA_real_,
  term_network = NULL, truncate_response = FALSE,
  cv = FALSE
) {
  # generate design matrix
  df_model <- create_model_df(df_genes, df_terms)

  # setup data
  covariates <- df_model %>%
    select(ends_with(".member")) %>%
    names

  X <- df_model %>% select(all_of(covariates)) %>% as.matrix
  Y <- df_model %>% select("pvalue") %>% as.matrix

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
    family = netReg::beta
  )

  # extract coefficients
  df_enrich <- as.data.frame(coef(fit)) %>%
    mutate(rowname = c("intercept", covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "name", "(.*).member") %>%
    rename(enrichment = `y[1]`)

  return(df_enrich)
}
