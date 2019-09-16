library(pareg)

# setup data
df.genes <- read.csv("genes.csv")
df.terms <- read.csv("terms.csv")

term.network <- as.matrix(read.csv("term_network.csv", row.names=1))

# run model
df.enr <- pareg::pareg(
    df.genes, df.terms,
    term.network=term.network, truncate.response=TRUE)

write.csv(df.enr, file="enrichment_result.csv", row.names=FALSE)
