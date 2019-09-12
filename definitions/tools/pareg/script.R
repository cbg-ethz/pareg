library(pareg)

# setup data
df.genes = read.csv("genes.csv")
df.terms = read.csv("terms.csv")

# run model
df.enr <- pareg::pareg(df.genes, df.terms, truncate.response=TRUE)

write.csv(df.enr, file="enrichment_result.csv", row.names=FALSE)
