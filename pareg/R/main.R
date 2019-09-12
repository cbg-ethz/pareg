library(tidyverse)
library(netReg)


pareg <- function (df.genes, df.terms, term.network=NULL) {
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
  fit <- netReg::edgenet(X, Y, G.X=term.network, family=mgcv::betar())

  df.enrich <- as.data.frame(coef(fit)) %>%
    mutate(rowname=c("intercept", covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "name", "(.*).member") %>%
    rename(enrichment=`y[1]`)

  return(df.enrich)
}
