source("code/utils.R")
init_public_plotting()

food_microbiome <- read_significant_edges(
  "data_analysis/plot_data/1_food_group_vs_microbiome_t1_t4/based_on_sspg/significant_cor_IR_IS.csv",
  node1_class = "food",
  node2_class = "microbiome"
)

food_metabolomics <- read_significant_edges(
  "data_analysis/plot_data/3_food_group_vs_metabolomics_t1_t4/based_on_sspg/significant_cor_IR_IS.csv",
  node1_class = "food",
  node2_class = "metabolomics"
)

save_network_plot(
  dplyr::bind_rows(food_microbiome, food_metabolomics),
  "data_analysis/4_food_microbiome_metabolomics_network/network_is.pdf",
  "Food-Microbiome-Metabolomics Network"
)

nutrition_microbiome <- read_significant_edges(
  "data_analysis/plot_data/nutrition_vs_microbiome_t1_t4/based_on_sspg/significant_cor_IR_IS.csv",
  node1_class = "nutrition",
  node2_class = "microbiome"
)

nutrition_metabolomics <- read_significant_edges(
  "data_analysis/plot_data/nutrition_vs_metabolomics_t1_t4/based_on_sspg/significant_cor_IR_IS.csv",
  node1_class = "nutrition",
  node2_class = "metabolomics"
)

save_network_plot(
  dplyr::bind_rows(nutrition_microbiome, nutrition_metabolomics),
  "data_analysis/7_nutrition_microbiome_metabolomics_network/network_is.pdf",
  "Nutrition-Microbiome-Metabolomics Network"
)
