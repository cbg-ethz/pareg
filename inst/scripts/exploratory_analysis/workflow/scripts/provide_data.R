library(BiocParallel)
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
multicoreParam <- MulticoreParam(workers = snakemake@threads)

tcga <- runDE(
  tcga,
  de.method = "limma",
  padj.method = "flexible",
  parallel = multicoreParam
)

# save results
writeDE(tcga, outdir)
