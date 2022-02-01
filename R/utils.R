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
#' @import tidyverse
#' @importFrom rlang .data
create_model_df <- function(df_genes, df_terms, pvalue_threshold = 0.05) {
  df_terms %>%
    group_by(term) %>%
    mutate(member = gene %in% df_genes$gene) %>%
    ungroup() %>%
    mutate_at(vars(gene), as.character) %>%
    right_join(
      df_genes %>% mutate_at(vars(gene), as.character),
      by = "gene"
    ) %>%
    pivot_wider(names_from = term, values_from = member, values_fill = FALSE) %>%
    select(
      -one_of("NA") # use `one_of` to handle case where NA does not exist (happens when a gene appears in no term)
    ) %>%
    rename_at(vars(-gene, -pvalue), ~ paste0(., ".member")) %>%
    mutate_at(
      vars(gene), factor # gene is character if select statement is executed
    ) %>%
    mutate(
      pvalue_sig = pvalue <= pvalue_threshold,
      pvalue_notsig = pvalue > pvalue_threshold
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
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter rename arrange desc
#' @importFrom tidyr extract
#' @importFrom stats coef
as.data.frame.pareg <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(row.names) || optional) {
    stop("row.names and optional arguments not supported")
  }

  as.data.frame(coef(x$obj)) %>% # nolint
    rownames_to_column() %>%
    mutate(rowname = c("intercept", x$covariates)) %>%
    filter(grepl(".member$", rowname)) %>%
    extract(rowname, "term", "(.*).member") %>%
    rename(enrichment = `y[1]`) %>%
    arrange(desc(enrichment))
}


#' @title Sample elements based on similarity structure.
#'
#' @description Choose similar object more often,
#' depending on `similarity_factor`.
#'
#' @export
#' @param sim_mat Similarity matrix with samples as row/col names.
#' @param size How many samples to draw.
#' @param similarity_factor Uniform sampling for 0. Weights mixture of uniform
#' and similarity vector for each draw.
#'
#' @return Vector of samples.
#'
#' @examples
#' similarity_sample(matrix(runif(100), nrow = 10, ncol = 10), 3)
#' @import progress
similarity_sample <- function(sim_mat, size, similarity_factor = 1) {
  # sanity check
  stopifnot(all(rownames(sim_mat) == colnames(sim_mat)))

  # handle matrices without row/col names
  if (is.null(rownames(sim_mat))) {
    sample_list <- seq_len(dim(sim_mat)[[1]])
    rownames(sim_mat) <- colnames(sim_mat) <- sample_list
  } else {
    # trying to extract numeric rownames out will result in characters!?
    sample_list <- rownames(sim_mat)
  }
  sample_num <- length(sample_list)

  # prepare sampling
  selected_samples <- c()
  current_sample <- sample(sample_list, 1)

  pb <- progress::progress_bar$new(total = size)
  pb$tick(0)

  # draw samples
  for (i in seq_len(size)) {
    pb$tick()

    prob <- (1 - similarity_factor) * rep(1, sample_num) / sample_num + similarity_factor * sim_mat[current_sample, ]
    current_sample <- sample(sample_list, 1, prob = prob)
    selected_samples <- c(selected_samples, current_sample)
  }

  selected_samples
}


#' @title Similarity matrix generation.
#'
#' @description Generate block-structured similarity matrices
#' corresponding to cluster structures.
#'
#' @export
#' @param cluster_sizes List of cluster sizes.
#'
#' @return Similarity matrix with samples as row-/colnames.
#'
#' @examples
#' generate_similarity_matrix(c(1, 2, 3))
#' @importFrom Matrix bdiag
generate_similarity_matrix <- function(cluster_sizes) {
  as.matrix(bdiag(
    lapply(cluster_sizes, function(x) matrix(runif(x * x), nrow = x, ncol = x))
  ))
}


#' @title Transform vector from [0, 1] to (0, 1).
#'
#' @description Make (response) vector conform to Beta assumptions
#' as described in section 2 of the betareg vignette
#' https://cran.r-project.org/web/packages/betareg/vignettes/betareg.pdf.
#'
#' @export
#' @param y Numeric vector in [0, 1]^N
#'
#' @return Numeric vector in (0, 1)^N
#'
#' @examples
#' transform_y(c(0, 0.5, 1))
transform_y <- function(y) {
  n_obs <- sum(!is.na(y))
  (y * (n_obs - 1) + 0.5) / n_obs
}
