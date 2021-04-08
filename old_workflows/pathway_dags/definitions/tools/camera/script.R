library(tidyverse)


# read data
fname.expr <- commandArgs(TRUE)[1]
df.cts <- read_csv(fname.expr) %>% column_to_rownames("node")

fname.info <- commandArgs(TRUE)[2]
df.info <- read_csv(fname.info) %>% column_to_rownames("sample")

# prepare data
design <- model.matrix(~ condition, data=df.info)
df.v <- limma::voom(df.cts, design)

# find pathways
index <- list()
for (node in rownames(df.cts)) {
  pw <- strsplit(node, "_")[[1]][1]
  index[[pw]] <- c(index[[pw]], node)
}

# compute enrichment
df.res <- limma::camera(df.v, index=index, design=design)

# store result
df.res %>% rownames_to_column("pathway") %>% write_csv("result.csv")
