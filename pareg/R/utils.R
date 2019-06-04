#' Create design matrix
#'
#' Store term membership for each gene
create_model_df <- function (df.genes, df.terms) {
  # TODO: fix factor/character warning
  df.model <- df.terms %>%
    group_by(name) %>%
    mutate(member=gene %in% df.genes$gene) %>%
    ungroup %>%
    right_join(df.genes, by="gene") %>%
    spread(key=name, value=member, fill=FALSE) %>%
    select(-one_of("<NA>")) %>% # use `one_of` to handle case where <NA> does not exist
    rename_at(vars(-gene, -pvalue), ~ paste0(., ".member")) %>%
    mutate_at(vars(gene), factor) # gene is character if select statement is executed

  return(df.model)
}
