#' @title Plot result of enrichment computation.
#'
#' @description Visualize pathway enrichments as network.
#'
#' @export
#'
#' @param x An object of class \code{pareg}.
#' @param show_term_names Whether to plot node labels.
#' @param min_similarity Don't plot edges for similarities below this value.
#'
#' @return ggplot object.
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
#' plot(fit)
#' @import tidyverse ggraph
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarize distinct pull
#' @importFrom magrittr %<>% extract2
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom ggplot2 aes theme element_rect coord_fixed
plot.pareg <- function(x, show_term_names = TRUE, min_similarity = 0) {
  # prepare data
  df_enr <- as.data.frame(x)
  df_terms <- x$df_terms
  term_network <- x$term_network

  term_sizes <- df_terms %>%
    group_by(term) %>%
    summarize(size = n())

  if (is.null(term_network)) {
    # pareg was run without network regularization
    term_list <- df_terms %>%
      distinct(term) %>%
      pull(term)
    term_network <- matrix(0, length(term_list), length(term_list))
    rownames(term_network) <- colnames(term_network) <- term_list
  }

  # subset term network
  term_network[term_network < min_similarity] <- 0

  # create plot
  term_graph <- as_tbl_graph(igraph::graph_from_adjacency_matrix( # nolint
    term_network,
    weighted = TRUE
  )) %>%
    activate(nodes) %>%
    mutate(
      enrichment = data.frame(term = .data$name) %>%
        left_join(df_enr, by = "term") %>%
        pull(enrichment),
      term_size = data.frame(term = .data$name) %>%
        left_join(term_sizes, by = "term") %>%
        pull(size),
    )

  edge_count <- term_graph %>%
    activate(edges) %>%
    as_tibble() %>%
    dim() %>%
    extract2(1)

  if (edge_count > 0) {
    term_graph %<>%
      activate(edges) %>%
      mutate(
        term_similarity = .data$weight
      )
  }

  p <- term_graph %>%
    ggraph(layout = "fr")

  if (edge_count > 0) {
    p <- p + geom_edge_link(aes(alpha = .data$term_similarity))
  }

  p <- p +
    geom_node_point(
      aes(size = .data$term_size, color = .data$enrichment)
    ) +
    scale_size(range = c(2, 10)) +
    scale_color_gradient2(
      low = "red", mid = "grey", high = "blue", midpoint = 0,
      na.value = "black"
    ) +
    scale_edge_alpha() +
    coord_fixed() +
    theme(
      panel.background = element_rect(fill = "white")
    )

  if (show_term_names) {
    p <- p + geom_node_text(aes(label = .data$name))
  }

  return(p)
}
