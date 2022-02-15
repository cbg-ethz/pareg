# pareg

[![lintr](https://github.com/cbg-ethz/pareg/actions/workflows/lintr.yaml/badge.svg)](https://github.com/cbg-ethz/pareg/actions/workflows/lintr.yaml)
[![Mega-Linter](https://github.com/cbg-ethz/pareg/actions/workflows/mega-linter.yaml/badge.svg)](https://github.com/cbg-ethz/pareg/actions/workflows/mega-linter.yaml)
[![check-bioc](https://github.com/cbg-ethz/pareg/actions/workflows/check-bioc.yaml/badge.svg)](https://github.com/cbg-ethz/pareg/actions/workflows/check-bioc.yaml)
[![pkgdown](https://github.com/cbg-ethz/pareg/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/cbg-ethz/pareg/actions/workflows/pkgdown.yaml)

Pathway enrichment computations using a regularized regression approach to incorporate inter-pathway relations in the statistical model.


## Installation

Install the latest stable version from Bioconductor:
```r
BiocManager::install("pareg")
```

Install the latest development version from GitHub:
```r
remotes::install_github("cbg-ethz/pareg")
```


## Project structure

* `.`: R package
* `inst/scripts/`: Snakemake workflows and other utilities


## Dev notes

* Update `NAMESPACE` and write man pages: `Rscript -e "devtools::document()"`
* Format code: `Rscript -e "styler::style_dir('.', transformers = biocthis::bioc_style())"`
