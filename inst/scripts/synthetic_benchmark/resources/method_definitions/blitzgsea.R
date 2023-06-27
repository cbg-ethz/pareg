# prepare environment
source(snakemake@params$setup_code_fname)

work_dir <- dirname(snakemake@output$fname)

gene_rank_fname <- file.path(work_dir, "genes.csv")
study$df %>%
  select(gene, pvalue) %>%
  write_csv(gene_rank_fname, col_names = FALSE)

term_fname <- file.path(work_dir, "terms.csv")
df_terms %>%
  write_csv(term_fname, col_names = FALSE)

# run blitzGSEA
fname_out <- file.path(work_dir, "raw_result.csv")
system2("python", paste(
    file.path(dirname(snakemake@input$script), "blitzgsea_code.py"),
    gene_rank_fname,
    term_fname,
    fname_out
))

# format result
df <- read_csv(fname_out)

df %>%
    head()

df %>%
    rename(term = Term, enrichment = nes) %>%
    select(term, enrichment) %>%
    mutate(
        method = "blitzgsea",
        enrichment = abs(enrichment)
    ) %>%
    write_csv(snakemake@output$fname)
