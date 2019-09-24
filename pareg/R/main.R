library(tidyverse)
library(netReg)


pareg <- function (
  df.genes, df.terms,
  lasso.param=0, network.param=0,
  term.network=NULL, truncate.response=FALSE
) {
  # generate design matrix
  df.model <- create_model_df(df.genes, df.terms)

  # setup data
  covariates <- df.model %>%
    select(ends_with(".member")) %>%
    names

  X <- df.model %>% select(covariates) %>% as.matrix
  Y <- df.model %>% select("pvalue") %>% as.matrix

  # truncate response if requested
  if (truncate.response) {
    eps <- .Machine$double.eps * 1e9 # see ?mgcv::betar
    Y[Y > 1 - eps] <- 1 - eps
    Y[Y < eps] <- eps
  }

  # fit model
  fit <- netReg::edgenet(
    X, Y,
    G.X=term.network,
    lambda=lasso.param, psigx=network.param, psigy=0,
    family=mgcv::betar())

  # extract coefficients
  df.enrich <- as.data.frame(coef(fit)) %>%
    mutate(rowname=c("intercept", covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "name", "(.*).member") %>%
    rename(enrichment=`y[1]`)

  return(df.enrich)
}
