#' @export
plot.pareg <- function(x, show_term_names = TRUE) {
  # prepare data
  df_enr <- as.data.frame(x)
  df_terms <- x$df_terms
  term_network <- x$term_network

  term_sizes <- df_terms %>%
    group_by(name) %>%
    summarize(size = n())

  if (is.null(term_network)) {
    # pareg was run without network regularization
    term_list <- df_terms %>% distinct(term) %>% pull(term)
    term_network <- matrix(0, length(term_list), length(term_list))
    rownames(term_network) <- colnames(term_network) <- term_list
  }

  # create plot
  # TODO: make sure that node order checks out
  p <- as_tbl_graph(igraph::graph_from_adjacency_matrix(  # nolint
    term_network, weighted = TRUE
  )) %>%
    activate(nodes) %>%
    mutate(
      enrichment = df_enr$enrichment,
      label = df_enr$name,
      termsize = term_sizes$size,
    ) %>%
    activate(edges) %>%
    mutate(

    ) %>%
    ggraph(layout = "mds") +
      geom_edge_link() +
      geom_node_point(
        aes(size = .data$termsize, color = .data$enrichment)
      ) +
      scale_size(range = c(2, 10)) +
      scale_color_gradient2(
        low = "red", mid = "grey", high = "blue", midpoint = 0,
        na.value = "black"
      ) +
      coord_fixed() +
      theme(
        panel.background = element_rect(fill = "white")
      )

  if (show_term_names) {
    p <- p + geom_node_text(aes(label = .data$label))
  }

  return(p)
}
