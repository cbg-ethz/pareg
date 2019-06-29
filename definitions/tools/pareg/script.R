library(pareg)

# setup data
df.genes = read.csv("genes.csv")
df.terms = read.csv("terms.csv")

# fix for "invalid dependent variable, all observations must be in (0, 1)"
df.genes[df.genes$pvalue == 1, "pvalue"] <- 1 - 1e-7

# run model
df.enr <- pareg::pareg(df.genes, df.terms)

write.csv(df.enr, file="enrichment_result.csv", row.names=FALSE)
