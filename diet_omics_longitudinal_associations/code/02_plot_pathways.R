source("code/utils.R")
init_public_plotting()

ratio_data <- readr::read_csv("data_analysis/plot_data/pathways/pathway_metabolite_ratio.csv", show_col_types = FALSE) %>%
  tidyr::pivot_longer(-Pathway, names_to = "group", values_to = "ratio")

p1 <- ratio_data %>%
  ggplot2::ggplot(ggplot2::aes(x = ratio, y = stats::reorder(Pathway, ratio), fill = group)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Pathway metabolite ratio", y = NULL, fill = NULL)
p1
dir_for("data_analysis/pathways/permutation_pathway_comparison.pdf")
ggplot2::ggsave(filename = "data_analysis/pathways/permutation_pathway_comparison.pdf", plot = p1, width = 10, height = 8)
ggplot2::ggsave(filename = "data_analysis/pathways/permutation_pathway_comparison.png", plot = p1, width = 10, height = 8)

top10_pathways <- ratio_data %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(max_ratio = base::max(ratio), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(max_ratio)) %>%
  dplyr::slice_head(n = 10)

p2 <- ratio_data %>%
  dplyr::filter(Pathway %in% top10_pathways[["Pathway"]]) %>%
  dplyr::mutate(max_ratio = top10_pathways$max_ratio[base::match(Pathway, top10_pathways$Pathway)]) %>%
  ggplot2::ggplot(ggplot2::aes(x = ratio, y = stats::reorder(Pathway, max_ratio), fill = group)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Pathway metabolite ratio", y = NULL, fill = NULL)
p2
dir_for("data_analysis/pathways/top10_permutation_pathway_comparison.pdf")
ggplot2::ggsave(filename = "data_analysis/pathways/top10_permutation_pathway_comparison.pdf", plot = p2, width = 10, height = 7)
ggplot2::ggsave(filename = "data_analysis/pathways/top10_permutation_pathway_comparison.png", plot = p2, width = 10, height = 7)
