library(tidyverse)
library(netReg)


pareg <- function (df.genes, df.terms, ...) {
  # generate design matrix
  df.model <- create_model_df(df.genes, df.terms)
  # print(df.model)

  # setup data
  covariates <- df.model %>%
    select(ends_with(".member")) %>%
    names

  X <- df.model %>% select(covariates) %>% as.matrix
  Y <- df.model[,"pvalue"] %>% as.matrix

  # fit model
  fit <- edgenet(X, Y, lambda=0, psigx=0, psigy=0, family=mgcv::betar())

  df.enrich <- as.data.frame(coef(fit)) %>%
    mutate(rowname=c("intercept", covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "name", "(.*).member") %>%
    rename(enrichment=`y[1]`)

  return(df.enrich)
}
