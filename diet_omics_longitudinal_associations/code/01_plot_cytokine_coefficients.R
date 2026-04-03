source("code/utils.R")
init_public_plotting()

data <- readr::read_csv("data_analysis/plot_data/coefficients.csv", show_col_types = FALSE) %>%
  dplyr::select(Y_Col, `Diet Pattern`, SSPG) %>%
  tidyr::pivot_longer(-Y_Col, names_to = "metric", values_to = "value")

p <- ggplot2::ggplot(data, ggplot2::aes(x = metric, y = Y_Col, fill = value)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.2) +
  ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank()
  )

dir_for("data_analysis/cytokine_coefficients/plot.pdf")
ggplot2::ggsave("data_analysis/cytokine_coefficients/plot.pdf", p, width = 6, height = 8)
