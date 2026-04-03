source("code/utils.R")
init_public_plotting()

plot_jobs <- tibble::tribble(
  ~input, ~output, ~title,
  "data_analysis/plot_data/nutrition_vs_microbiome_t1_t4/based_on_sspg/cor_data.csv", "data_analysis/nutrition_vs_microbiome_t1_t4/based_on_sspg/all_cor.pdf", "Nutrition vs Microbiome (All)",
  "data_analysis/plot_data/nutrition_vs_microbiome_t1_t4/based_on_sspg/cor_data_IR.csv", "data_analysis/nutrition_vs_microbiome_t1_t4/based_on_sspg/IR_cor.pdf", "Nutrition vs Microbiome (IR)",
  "data_analysis/plot_data/nutrition_vs_microbiome_t1_t4/based_on_sspg/cor_data_IS.csv", "data_analysis/nutrition_vs_microbiome_t1_t4/based_on_sspg/IS_cor.pdf", "Nutrition vs Microbiome (IS)",
  "data_analysis/plot_data/nutrition_vs_metabolomics_t1_t4/based_on_sspg/cor_data.csv", "data_analysis/nutrition_vs_metabolomics_t1_t4/based_on_sspg/all_cor.pdf", "Nutrition vs Metabolomics (All)",
  "data_analysis/plot_data/nutrition_vs_metabolomics_t1_t4/based_on_sspg/cor_data_IR.csv", "data_analysis/nutrition_vs_metabolomics_t1_t4/based_on_sspg/IR_cor.pdf", "Nutrition vs Metabolomics (IR)",
  "data_analysis/plot_data/nutrition_vs_metabolomics_t1_t4/based_on_sspg/cor_data_IS.csv", "data_analysis/nutrition_vs_metabolomics_t1_t4/based_on_sspg/IS_cor.pdf", "Nutrition vs Metabolomics (IS)",
  "data_analysis/plot_data/1_food_group_vs_microbiome_t1_t4/based_on_sspg/significant_cor_IR_IS.csv", "data_analysis/1_food_group_vs_microbiome_t1_t4/based_on_sspg/all_cor.pdf", "Food Group vs Microbiome",
  "data_analysis/plot_data/3_food_group_vs_metabolomics_t1_t4/based_on_sspg/significant_cor_IR_IS.csv", "data_analysis/3_food_group_vs_metabolomics_t1_t4/based_on_sspg/all_cor.pdf", "Food Group vs Metabolomics",
  "data_analysis/plot_data/nutrition_vs_microbiome/based_on_sspg/cor_data_IS.csv", "data_analysis/nutrition_vs_microbiome/based_on_sspg/IS_cor.pdf", "Nutrition vs Microbiome (IS)",
  "data_analysis/plot_data/nutrition_vs_metabolomics/based_on_a1c/cor_data_normal.csv", "data_analysis/nutrition_vs_metabolomics/based_on_a1c/Normal_cor.pdf", "Nutrition vs Metabolomics (Normal)"
)

for (i in base::seq_len(base::nrow(plot_jobs))) {
  save_correlation_plot(
    input_csv = plot_jobs$input[i],
    output_pdf = plot_jobs$output[i],
    title = plot_jobs$title[i]
  )
}
