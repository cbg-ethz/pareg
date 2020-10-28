library(tidyverse)


# parameters
fname <- snakemake@output$fname

# msigdb statistics
msigdbr::msigdbr_collections()

# retrieve raw data
df_all <- msigdbr::msigdbr(species = "Homo sapiens")

# save raw data
df_all %>%
  write_csv(fname)
