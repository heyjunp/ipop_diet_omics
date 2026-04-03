source("code/utils.R")
init_public_plotting()

sync_files <- base::list.files(
  "data_analysis/plot_data/8_feature_shape_intepretation/1_synchronized_response/1_nutrition_vs_microbiome",
  pattern = "cor_data\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

delay_files <- base::list.files(
  "data_analysis/plot_data/8_feature_shape_intepretation/2_delayed_response/1_nutrition_vs_microbiome",
  pattern = "cor_data\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

plot_feature_summary <- function(files, output) {
  df <- purrr::map_dfr(files, function(path) {
    dat <- readr::read_csv(path, show_col_types = FALSE)
    base::data.frame(
      comparison = base::basename(base::dirname(path)),
      n_sig = base::sum(dat$p_adjust < 0.05, na.rm = TRUE)
    )
  })

  p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(comparison, n_sig), y = n_sig)) +
    ggplot2::geom_col(fill = "#4c78a8") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Significant correlations")

  ggplot2::ggsave(output, p, width = 8, height = 6)
}

plot_feature_summary(
  sync_files,
  "data_analysis/8_feature_shape_intepretation/1_synchronized_response/1_nutrition_vs_microbiome/summary_counts.pdf"
)

plot_feature_summary(
  delay_files,
  "data_analysis/8_feature_shape_intepretation/2_delayed_response/1_nutrition_vs_microbiome/summary_counts.pdf"
)
