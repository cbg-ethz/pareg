#' @title Create design matrix.
#'
#' @description Store term membership for each gene.
#'
#' @export
#'
#' @param df_genes Dataframe storing gene names and DE p-values.
#' @param df_terms Dataframe storing pathway database.
#' @param pvalue_threshold P-value threshold to create binary columns
#' `pvalue_sig` and `pvalue_notsig`.
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
#' @importFrom dplyr group_by mutate ungroup mutate_at
#' @importFrom dplyr vars right_join select any_of rename_at
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom rlang .data
create_model_df <- function(df_genes, df_terms, pvalue_threshold = 0.05) {
  df_terms %>%
    group_by(.data$term) %>%
    mutate(member = .data$gene %in% df_genes$gene) %>%
    ungroup() %>%
    mutate_at(vars(.data$gene), as.character) %>%
    right_join(
      df_genes %>% mutate_at(vars(.data$gene), as.character),
      by = "gene"
    ) %>%
    pivot_wider(
      names_from = .data$term,
      values_from = .data$member,
      values_fill = FALSE
    ) %>%
    select(
      # use `any_of` to handle case where NA does not exist
      # (happens when a gene appears in no term)
      -any_of("NA")
    ) %>%
    rename_at(vars(-.data$gene, -.data$pvalue), ~ paste0(., ".member")) %>%
    mutate_at(
      vars(.data$gene),
      factor # gene is character if select statement is executed
    ) %>%
    mutate(
      pvalue_sig = .data$pvalue <= pvalue_threshold,
      pvalue_notsig = .data$pvalue > pvalue_threshold
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
#' fit <- pareg(df_genes, df_terms, max_iterations = 10)
#' as.data.frame(fit)
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter rename arrange desc
#' @importFrom tidyr extract
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom stats coef
as.data.frame.pareg <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(row.names) || optional) {
    stop("row.names and optional arguments not supported")
  }

  as.data.frame(coef(x$obj)) %>% # nolint
    rownames_to_column() %>%
    mutate(rowname = c("intercept", x$covariates)) %>%
    filter(grepl(".member$", .data$rowname)) %>%
    extract(.data$rowname, "term", "(.*).member") %>%
    rename(enrichment = .data$`y[1]`) %>%
    arrange(desc(.data$enrichment))
}


#' @title Convert object of class \code{pareg} to class \code{enrichResult}.
#'
#' @description The resulting object can be passed to any method from the
#'  enrichplot package and thus allows for nice visualizations of the
#'  enrichment results. Note: term similarities are included if available.
#'
#' @export
#'
#' @param x An object of class \code{pareg}.
#' @param pvalue_threshold Treshold to select genes for count statistics.
#'
#' @return Object of class \code{enrichResult}.
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
#' fit <- pareg(df_genes, df_terms, max_iterations = 10)
#' as_enrichplot_object(fit)
#' @importFrom rlang .data
#' @importFrom methods new
#' @importFrom dplyr mutate filter rename inner_join rowwise tally n
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom purrr pluck
#' @importFrom stringr str_c
#' @importClassesFrom DOSE enrichResult
as_enrichplot_object <- function(x, pvalue_threshold = 0.05) {
  # find genes most relevant to enrichment
  sig_genes <- x$df_genes %>%
    filter(.data$pvalue <= pvalue_threshold) %>%
    pull(.data$gene)

  # helper functions
  get_term_genes <- function(df, variable) {
    df %>%
      filter(variable == .data$term) %>%
      pull(.data$gene) %>%
      str_c(collapse = "/")
  }

  sig_member_count <- function(df, variable) {
    df %>%
      filter(
        variable == .data$term &
          .data$gene %in% sig_genes
      ) %>%
      pluck(dim, 1)
  }

  # create class
  cls <- new(
    "enrichResult",
    result = x %>%
      as.data.frame() %>%
      inner_join(
        x$df_terms %>%
          group_by(.data$term) %>%
          tally(),
        by = "term"
      ) %>%
      rename(term_size = n) %>%
      rowwise() %>%
      mutate(
        ID = .data$term,
        Description = .data$term,
        p.adjust = .data$enrichment,
        Count = sig_member_count(x$df_terms, .data$ID),
        GeneRatio = paste0(.data$Count, "/", .data$term_size),
        geneID = get_term_genes(x$df_terms, .data$ID)
      ) %>%
      column_to_rownames("term"),
    geneSets = x$df_terms %>%
      pipe_split("term", "gene"),
    gene = as.character(sig_genes),
    readable = FALSE,
    keytype = "UNKNOWN",
    ontology = "UNKNOWN"
  )

  if (!is.null(x$term_network)) {
    cls@method <- "pareg"

    cls@termsim <- x$term_network
    cls@termsim[lower.tri(cls@termsim, diag = TRUE)] <- NA
  }

  cls
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
#' @importFrom progress progress_bar
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

  pb <- progress_bar$new(total = size)
  pb$tick(0)

  # draw samples
  for (i in seq_len(size)) {
    pb$tick()

    prob <- (1 - similarity_factor) *
      rep(1, sample_num) / sample_num +
      similarity_factor * sim_mat[current_sample, ]
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
#' @importFrom stats runif
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


#' @noRd
pipe_split <- function(df, group_col, value_col) {
  df <- as.data.frame(df)
  split(df[, value_col], df[, group_col])
}


#' @title Convert matrices.
#'
#' @description Convert sparse similarity matrix from package data to a
#' dense version with 1 on its diagonal.
#' This matrix can then be used by \code{pareg}.
#'
#' @export
#' @param mat_sparse Sparse matrix.
#'
#' @return Dense matrix
#'
#' @examples
#' transform_y(c(0, 0.5, 1))
#' @importFrom Matrix forceSymmetric
as_dense_sim <- function(mat_sparse) {
  mat_sparse <- forceSymmetric(mat_sparse, uplo = "U")
  mat_dense <- as.matrix(mat_sparse)
  diag(mat_dense) <- 1
  mat_dense
}
