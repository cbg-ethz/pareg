library(GSEABenchmarkeR)


# parameters
outdir <- snakemake@output$outdir

# download data
tcga <- loadEData(
  "tcga",
  nr.datasets = 2,
  mode = "ehub",
  data.dir = file.path(outdir, "data")
)

# run DEA
tcga <- runDE(
  tcga,
  de.method = "limma",
  padj.method = "flexible"
)

# save results
writeDE(tcga, outdir)
