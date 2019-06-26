library(dplyr)
library(tidyr)
library(tibble)
library(rstanarm)

pareg <- function (df.genes, df.terms, ...) {
  # generate design matrix
  df.model <- create_model_df(df.genes, df.terms)
  # print(df.model)

  # create formula
  covariates <- df.model %>%
    select(ends_with(".member")) %>%
    names
  form <- as.formula(paste("pvalue ~", paste(covariates, collapse="+")))

  # fit model
  bfit <- rstanarm::stan_betareg(
    formula=form, data=df.model,
    ...
  )
  # print(summary(bfit))

  # extract enrichments
  df.enrich <- as.data.frame(coef(bfit)) %>%
    rownames_to_column %>%
    filter(grepl(".memberTRUE$", rowname)) %>%
    extract(rowname, "name", "(.*).memberTRUE") %>%
    rename(enrichment=`coef(bfit)`)

  return(df.enrich)
}
