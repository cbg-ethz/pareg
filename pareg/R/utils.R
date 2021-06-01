#' Create design matrix
#'
#' Store term membership for each gene
create_model_df <- function(df_genes, df_terms) {
  df_terms %>%
    group_by(name) %>%
    mutate(member = gene %in% df_genes$gene) %>%
    ungroup %>%
    mutate_at(vars(gene), as.character) %>%
    right_join(
      df_genes %>% mutate_at(vars(gene), as.character), by = "gene"
    ) %>%
    pivot_wider(names_from = name, values_from = member, values_fill = FALSE) %>%
    select(
      -one_of("<NA>") # use `one_of` to handle case where <NA> does not exist
    ) %>%
    rename_at(vars(-gene, -pvalue), ~ paste0(., ".member")) %>%
    mutate_at(
      vars(gene), factor # gene is character if select statement is executed
    )
}


#' @export
as.data.frame.pareg <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(row.names) || optional) {
    stop("row.names and optional arguments not supported")
  }

  as.data.frame(coef(x$fit)) %>%  # nolint
    mutate(rowname = c("intercept", x$covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "name", "(.*).member") %>%
    rename(enrichment = `y[1]`)
}
