#' @title Create design matrix.
#'
#' @description Store term membership for each gene.
#'
#' @export
#'
#' @param df_genes Dataframe storing gene names and DE p-values.
#' @param df_terms Dataframe storing pathway database.
#'
#' @return Dataframe.
#'
#' @examples
#' df_genes <- data.frame(
#'   gene = c("g1", "g2"),
#'   pvalue = c(0.1, 0.2)
#' )
#' df_terms <- data.frame(
#'   term = c("A", "A", "B", "B", "C"),
#'   gene = c("g1", "g2", "g1", "g2", "g2")
#' )
#' create_model_df(df_genes, df_terms)
#'
#' @import tidyverse
#' @importFrom rlang .data
create_model_df <- function(df_genes, df_terms) {
  df_terms %>%
    group_by(term) %>%
    mutate(member = gene %in% df_genes$gene) %>%
    ungroup %>%
    mutate_at(vars(gene), as.character) %>%
    right_join(
      df_genes %>% mutate_at(vars(gene), as.character), by = "gene"
    ) %>%
    pivot_wider(names_from = term, values_from = member, values_fill = FALSE) %>%
    select(
      -one_of("NA") # use `one_of` to handle case where NA does not exist (happens when a gene appears in no term)
    ) %>%
    rename_at(vars(-gene, -pvalue), ~ paste0(., ".member")) %>%
    mutate_at(
      vars(gene), factor # gene is character if select statement is executed
    )
}


#' @title as.data.frame for an object of class \code{pareg}.
#'
#' @description Retrieve dataframe with enrichment information.
#'
#' @export
#'
#' @param x An object of class \code{pareg}.
#' @param row.names Optional character vector of rownames.
#' @param optional Allow optional arguments.
#' @param ... Additional arguments.
#'
#' @return Dataframe containing enrichment score and name for each pathway.
#'
#' @examples
#' df_genes <- data.frame(
#'   gene = paste("g", 1:20, sep = ""),
#'   pvalue = c(
#'     rbeta(10, .1, 1),
#'     rbeta(10, 1, 1)
#'   )
#' )
#' df_terms <- rbind(
#'   data.frame(
#'     term = "foo",
#'     gene = paste("g", 1:10, sep = "")
#'   ),
#'   data.frame(
#'     term = "bar",
#'     gene = paste("g", 11:20, sep = "")
#'   )
#' )
#' fit <- pareg(df_genes, df_terms)
#' as.data.frame(fit)
#'
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter rename arrange desc
#' @importFrom tidyr extract
#' @importFrom stats coef
as.data.frame.pareg <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(row.names) || optional) {
    stop("row.names and optional arguments not supported")
  }

  as.data.frame(coef(x$fit)) %>%  # nolint
    rownames_to_column %>%
    mutate(rowname = c("intercept", x$covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "term", "(.*).member") %>%
    rename(enrichment = `y[1]`) %>%
    arrange(desc(enrichment))
}
