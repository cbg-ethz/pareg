library(pareg)

# setup data
df_genes <- read.csv("genes.csv")
df_terms <- read.csv("terms.csv")

# run model
df_enr <- pareg::pareg(df_genes, df_terms, truncate_response = TRUE)

write.csv(df_enr, file = "enrichment_result.csv", row.names = FALSE)
