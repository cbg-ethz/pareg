library(pareg)

# setup data
df.genes = read.csv("genes.csv")
df.terms = read.csv("terms.csv")

# run model
df.enr <- pareg::pareg(df.genes, df.terms)

write.csv(df.enr, file="result.csv", row.names=FALSE)
