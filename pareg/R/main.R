library(tidyr)
library(dplyr)

# library(betareg)


pareg <- function (df.genes, df.terms) {
  df.model <- df.genes %>%
    mutate(
      member=gene %in% df.terms$gene
    )
  # print(df.model)

  # bfit <- betareg::betareg(
  #   pvalue ~ member, df.model,
  #   link="logit", type="ML",
  #   control=betareg::betareg.control(method="L-BFGS-B", fsmaxit=0))

  bfit <- rstanarm::stan_betareg(pvalue ~ member, df.model)
  # print(summary(bfit))

  enrich <- coef(bfit)["memberTRUE"]

  return(data.frame(
    name=df.terms$name[1],
    enrichment=enrich
  ))
}
