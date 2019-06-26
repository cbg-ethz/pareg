library(tidyverse)

library(SummarizedExperiment)

library(GSEABenchmarkeR)
library(EnrichmentBrowser)

## differential expression
# load data
geo2kegg <- GSEABenchmarkeR::loadEData("geo2kegg")
geo2kegg <- GSEABenchmarkeR::maPreproc(geo2kegg) # probes -> genes (ENTREZ Gene IDs [EnrichmentBrowser::idMap])
names(geo2kegg)

#tcga <- loadEData("tcga", nr.datasets=2)

# run differential expression analysis
geo2kegg <- GSEABenchmarkeR::runDE(geo2kegg, de.method="limma", padj.method="flexible")
sapply(geo2kegg, function (x) dim(rowData(x)))

df.dea <- purrr::imap(geo2kegg, function (entry, dataset.name) {
  rowData(entry) %>%
    as.data.frame %>%
    rownames_to_column("gene") %>%
    rename(ADJ.PVAL="p_value") %>%
    dplyr::select(gene, p_value) %>%
    mutate(dataset=dataset.name)
}) %>%
  map_df(bind_rows)

## genesets
kegg.gs <- getGenesets(org="hsa", db="kegg")
df.gs <- kegg.gs %>%
  enframe %>%
  unnest %>%
  rename(name="term", value="gene") %>%
  extract(term, "term")

## disease relevance rankings
# geneset ranking per disease
data.dir <- system.file("extdata", package="GSEABenchmarkeR")
mala.kegg.file <- file.path(data.dir, "malacards", "KEGG.rds")
mala.kegg <- readRDS(mala.kegg.file)
sapply(mala.kegg, nrow)

# associate disease with dataset
d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
d2d.map <- readDataId2diseaseCodeMap(d2d.file)
head(d2d.map)

## save to disk
root.dir <- snakemake@output[["out_dir"]]
store_csv <- function (data, group) {
  # setup
  dataset.name <- group[[1]]
  dname <- file.path(root.dir, dataset.name)

  unlink(dname, recursive=TRUE)
  dir.create(dname, recursive=TRUE)

  # save DE and GS data
  data %>% write_csv(file.path(dname, "input.csv"))
  df.gs %>% write_csv(file.path(dname, "terms.csv"))

  # save validation data (expected disease-pathway associations)
  disease <- d2d.map[[dataset.name]]
  mala.kegg[[disease]] %>%
    as.data.frame %>%
    rownames_to_column("term") %>%
    rename(REL.SCORE="relevance.score") %>%
    dplyr::select(term, relevance.score) %>%
    write_csv(file.path(dname, "expected_terms.csv"))
}

df.dea %>%
  group_by(dataset) %>%
  group_walk(store_csv)
