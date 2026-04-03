project_root <- function() {
  if (requireNamespace("r4projects", quietly = TRUE)) {
    return(r4projects::get_project_wd())
  }
  base::getwd()
}

init_public_plotting <- function() {
  base::setwd(project_root())
  base::invisible(project_root())
}

dir_for <- function(path) {
  base::dir.create(base::dirname(path), recursive = TRUE, showWarnings = FALSE)
}

choose_feature_label <- function(df) {
  if ("true_name" %in% base::names(df)) {
    label <- df$true_name
    label[base::is.na(label) | label == ""] <- df$data_set2[base::is.na(label) | label == ""]
    return(label)
  }
  if ("Metabolite" %in% base::names(df)) {
    label <- df$Metabolite
    label[base::is.na(label) | label == ""] <- df$data_set2[base::is.na(label) | label == ""]
    return(label)
  }
  df$data_set2
}

make_correlation_plot <- function(df, title = NULL, top_n = 60) {
  base::stopifnot(base::all(c("data_set1", "cor") %in% base::names(df)))

  df$feature_label <- choose_feature_label(df)

  if ("p_adjust" %in% base::names(df)) {
    df <- df[base::order(df$p_adjust, -base::abs(df$cor)), , drop = FALSE]
  } else {
    df <- df[base::order(-base::abs(df$cor)), , drop = FALSE]
  }

  if (base::nrow(df) > top_n) {
    df <- df[base::seq_len(top_n), , drop = FALSE]
  }

  df$data_set1 <- base::factor(df$data_set1, levels = base::unique(df$data_set1))
  df$feature_label <- base::factor(df$feature_label, levels = base::rev(base::unique(df$feature_label)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = data_set1, y = feature_label, fill = cor)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title, fill = "cor")

  if ("class" %in% base::names(df)) {
    p <- p + ggplot2::facet_wrap(~class)
  }

  p
}

save_correlation_plot <- function(input_csv, output_pdf, title = NULL, top_n = 60) {
  df <- readr::read_csv(input_csv, show_col_types = FALSE)
  p <- make_correlation_plot(df, title = title, top_n = top_n)
  dir_for(output_pdf)
  ggplot2::ggsave(output_pdf, p, width = 10, height = 8)
}

read_significant_edges <- function(path, node1_class, node2_class, node1_label = "data_set1") {
  df <- readr::read_csv(path, show_col_types = FALSE)
  df$node1 <- df$data_set1
  df$node2 <- df$data_set2
  df$class1 <- node1_class
  df$class2 <- node2_class
  df$node_name1 <- df[[node1_label]]
  df$node_name2 <- choose_feature_label(df)
  df
}

save_network_plot <- function(edges, output_path, title = NULL) {
  nodes <- base::unique(base::rbind(
    base::data.frame(node = edges$node1, node_name = edges$node_name1, class = edges$class1),
    base::data.frame(node = edges$node2, node_name = edges$node_name2, class = edges$class2)
  ))

  graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = base::data.frame(
      from = base::match(edges$node1, nodes$node),
      to = base::match(edges$node2, nodes$node),
      cor = edges$cor
    ),
    directed = FALSE
  ) |>
    dplyr::mutate(Degree = tidygraph::centrality_degree(mode = "all"))

  p <- ggraph::ggraph(graph, layout = "fr") +
    ggraph::geom_edge_link(color = "#8a8a8a", alpha = 0.5, show.legend = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(fill = class, size = Degree), shape = 21, color = "black") +
    ggraph::geom_node_text(ggplot2::aes(label = node_name), repel = TRUE, size = 2.5) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title)

  dir_for(output_path)
  ggplot2::ggsave(output_path, p, width = 8, height = 8)
}
