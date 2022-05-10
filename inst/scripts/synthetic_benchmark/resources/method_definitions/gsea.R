# prepare environment
source(snakemake@params$setup_code_fname)

work_dir <- dirname(snakemake@output$fname)

gene_rank_fname <- file.path(work_dir, "genes.rnk")
study$df %>%
  select(gene, pvalue) %>%
  write_tsv(gene_rank_fname, col_names = FALSE)

term_fname <- file.path(work_dir, "terms.gmx")
df_terms %>%
  pivot_wider(names_from = term, values_from = gene, values_fn = list) %>%
  pivot_longer(everything()) %>%
  mutate(value = map(value, `length<-`, max(lengths(value)))) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  unnest(everything()) %>%
  write_tsv(term_fname, na = "")

results_dir <- file.path(work_dir, "results")
gsea_exec <- file.path(dirname(snakemake@input$script), "gsea-3.0.jar")

# run model
system2("java", paste(
  "-cp", gsea_exec,
  "-Xmx512m", "xtools.gsea.GseaPreranked",
  "-gmx", term_fname,
  "-norm", "", "-nperm", 1000, "-scoring_scheme", "classic",
  "-rpt_label", "my_analysis",
  "-create_svgs", "false", "-make_sets", "true", "-plot_top_x", 20,
  "-rnd_seed", 42,
  "-set_max", 9999999, "-set_min", 0,
  "-zip_report", "false",
  "-out", results_dir,
  "-gui", "false",
  "-rnk", gene_rank_fname
))

# extract results
analysis_dirs <- Sys.glob(file.path(results_dir, "/*"))
stopifnot(length(analysis_dirs) == 1)

df <- Sys.glob(file.path(analysis_dirs[[1]], "/gsea_report_for_na_*.xls")) %>%
  map_dfr(function(path) {
    read_tsv(path)
  }) %>%
  rename(term = NAME, pvalue = "NOM p-val") %>%
  mutate(enrichment = ifelse(
    pvalue > 0,
    -log10(pvalue),
    -log10(min(pvalue[pvalue > 0]))
  )) %>%
  select(term, enrichment)

# finalize
df <- df %>%
  mutate(
    method = "gsea",
    term = str_to_lower(term)
  )

df %>%
  head()

# save result
df %>%
  write_csv(snakemake@output$fname)
