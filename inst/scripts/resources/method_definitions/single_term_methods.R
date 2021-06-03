# prepare environment
source(snakemake@params$setup_code_fname)

# run models
df <- df_terms %>%
  group_by(term) %>%
  group_map(function(sub, key) {
    term <- key %>% pull(term)
    term_genes <- sub %>% pull(gene)

    nonterm_genes <- setdiff(all_genes, term_genes)

    # fisher's exact test
    cont_table <- matrix(c(
      length(intersect(term_genes, study_genes)), length(intersect(term_genes, nonstudy_genes)),
      length(intersect(nonterm_genes, study_genes)), length(intersect(nonterm_genes, nonstudy_genes))
    ), 2, 2)

    fet_p_value <- fisher.test(cont_table)$p.value
    # print(paste(term, p_value))

    # store results
    bind_rows(
      data.frame(method = "FET", term = term, enrichment = -log10(fet_p_value))
    )
  }) %>%
  { do.call(rbind.data.frame, .) }

df %>%
  head

# save result
df %>%
  write_csv(snakemake@output$fname)
