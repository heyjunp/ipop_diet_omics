library(dplyr); library(readr); library(readxl)
input_dir  <- "../data"
output_dir <- "."
metab_file <- file.path(input_dir, "metabolites_kw_test.csv")
diet_file  <- file.path(input_dir, "merged_diet_demographics.csv")
supp_file  <- file.path(input_dir, "Supplementary_Table_12.xlsx")

bh1 <- 0.1; bh2 <- 0.2
hm_pal  <- c("#142621", "#71d9bc", "#FFFFFF", "#e39d5f", "#4a331f")
vrange  <- c(-0.55, 0.55)
excl_diet <- c("Fol_DFE", "VitD_mcg", "VitB_NE", "VitA_RAE", "Water", "OCarb",
               "Fib", "SolFib", "FatCals", "SatCals", "Sugar", "Disacc", "Caroten", "Chrom")

#partial spearman
pspear <- function(x, y, covariates) {
  valid <- complete.cases(x, y, covariates)
  if (sum(valid) < 10) return(c(rho = NA_real_, pvalue = NA_real_))
  x_v <- rank(x[valid])
  y_v <- rank(y[valid])
  cov_mat <- cbind(1, as.matrix(covariates[valid, , drop = FALSE]))
  x_resid <- residuals(lm.fit(cov_mat, x_v))
  y_resid <- residuals(lm.fit(cov_mat, y_v))
  ct <- cor.test(x_resid, y_resid, method = "pearson")
  c(rho = unname(ct$estimate), pvalue = ct$p.value)
}


metabolites <- read_csv(metab_file, show_col_types = FALSE)
diet_demo   <- read_csv(diet_file, show_col_types = FALSE)
supp_table  <- read_excel(supp_file)

metabolites$SubjectID <- sub("^(\\d+-\\d+).*", "\\1", metabolites$SampleID)
metab_cols <- grep("^(pHILIC|nHILIC|pRPLC|nRPLC)", names(metabolites), value = TRUE)

mmed <- metabolites %>%
  group_by(SubjectID) %>%
  summarise(across(all_of(metab_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop")

dietary_cols <- setdiff(names(diet_demo), c("SubjectID", "Sex", "Age"))
merged <- inner_join(mmed, diet_demo, by = "SubjectID")

##correlation loop
n_tests <- length(dietary_cols) * length(metab_cols)
merged$Sex_numeric <- as.integer(merged$Sex == "M")
covars <- merged[, c("Age", "Sex_numeric")]

results_list <- vector("list", n_tests)
idx <- 0L
for (di in seq_along(dietary_cols)) {
  dc <- dietary_cols[di]
  x_vals <- merged[[dc]]
  for (mc in metab_cols) {
    idx <- idx + 1L
    res <- pspear(x_vals, merged[[mc]], covars)
    results_list[[idx]] <- data.frame(
      dietary_component = dc, metabolite = mc,
      correlation = res["rho"], pvalue = res["pvalue"],
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
}
results_df <- bind_rows(results_list)

#BH correction
valid_mask <- !is.na(results_df$pvalue)
results_df$pvalue_adj <- NA_real_
results_df$pvalue_adj[valid_mask] <- p.adjust(results_df$pvalue[valid_mask], method = "BH")

sig_01 <- results_df[!is.na(results_df$pvalue_adj) & results_df$pvalue_adj < bh1, ]
sig_02 <- results_df[!is.na(results_df$pvalue_adj) & results_df$pvalue_adj < bh2, ]

##### annotate metabolite info
metab_map <- setNames(supp_table$true_name, supp_table$variable_id)
super_map <- setNames(supp_table$Super.pathway, supp_table$variable_id)
sub_map   <- setNames(supp_table$Sub.pathway, supp_table$variable_id)

results_df$metabolite_name <- metab_map[results_df$metabolite]
results_df$super_pathway   <- super_map[results_df$metabolite]
results_df$sub_pathway     <- sub_map[results_df$metabolite]
sig_01$metabolite_name <- metab_map[sig_01$metabolite]
sig_01$super_pathway   <- super_map[sig_01$metabolite]
sig_01$sub_pathway     <- sub_map[sig_01$metabolite]
sig_02$metabolite_name <- metab_map[sig_02$metabolite]
sig_02$super_pathway   <- super_map[sig_02$metabolite]
sig_02$sub_pathway     <- sub_map[sig_02$metabolite]

write_csv(results_df, file.path(output_dir, "all_correlations.csv"))

corr_matrix <- results_df %>%
  select(dietary_component, metabolite, correlation) %>%
  tidyr::pivot_wider(names_from = metabolite, values_from = correlation) %>%
  tibble::column_to_rownames("dietary_component")

pval_raw <- results_df %>%
  select(dietary_component, metabolite, pvalue) %>%
  tidyr::pivot_wider(names_from = metabolite, values_from = pvalue) %>%
  tibble::column_to_rownames("dietary_component")

pval_bh <- results_df %>%
  select(dietary_component, metabolite, pvalue_adj) %>%
  tidyr::pivot_wider(names_from = metabolite, values_from = pvalue_adj) %>%
  tibble::column_to_rownames("dietary_component")

write.csv(corr_matrix, file.path(output_dir, "correlation_matrix.csv"))
write.csv(pval_raw, file.path(output_dir, "pvalue_matrix_raw.csv"))
write.csv(pval_bh, file.path(output_dir, "pvalue_matrix_BH_adjusted.csv"))
write_csv(sig_01, file.path(output_dir, "significant_correlations_BH_p01.csv"))
write_csv(sig_02, file.path(output_dir, "significant_correlations_BH_p02.csv"))

summary_df <- data.frame(
  Analysis = "Partial Spearman Correlation", Covariates = "Age, Sex",
  N_subjects = nrow(merged), N_dietary = length(dietary_cols),
  N_metabolites = length(metab_cols), N_tests = nrow(results_df),
  N_sig_BH01 = nrow(sig_01), N_sig_BH02 = nrow(sig_02),
  N_diet_BH01 = n_distinct(sig_01$dietary_component),
  N_metab_BH01 = n_distinct(sig_01$metabolite),
  N_diet_BH02 = n_distinct(sig_02$dietary_component),
  N_metab_BH02 = n_distinct(sig_02$metabolite),
  stringsAsFactors = FALSE
)
write_csv(summary_df, file.path(output_dir, "analysis_summary.csv"))

#############################
######## heatmap ########
###########################

corr_filt <- corr_matrix[!rownames(corr_matrix) %in% excl_diet, , drop = FALSE]
pval_filt <- pval_bh[!rownames(pval_bh) %in% excl_diet, , drop = FALSE]

sig_mask <- pval_filt < bh2
sig_mask[is.na(sig_mask)] <- FALSE
corr_filt <- corr_filt[rowSums(sig_mask) > 0, colSums(sig_mask) > 0, drop = FALSE]
pval_filt <- pval_filt[rownames(corr_filt), colnames(corr_filt), drop = FALSE]

# ordering
nutrient_order <- c(
  "percent_carb_cal", "percent_pro_cal", "Pro", "Carb", "SugAdd",
  "TransFat", "Chol", "Retinol", "VitC", "Fluor", "Iodine", "Mn",
  "Moly", "Se", "W3", "Caff", "Chln"
)
food_order <- c("Coffee", "Meat_Poultry", "Sauces, Spices, and Herbs",
                "Seafood", "Legume_Products")
row_ord <- c(nutrient_order, food_order)
row_ord <- row_ord[row_ord %in% rownames(corr_filt)]
row_ord <- c(row_ord, setdiff(rownames(corr_filt), row_ord))

corr_filt <- corr_filt[row_ord, , drop = FALSE]
pval_filt <- pval_filt[row_ord, , drop = FALSE]
sep_idx <- which(rownames(corr_filt) %in% food_order)[1]

###
new_cols <- ifelse(is.na(metab_map[colnames(corr_filt)]),
                   colnames(corr_filt), metab_map[colnames(corr_filt)])
colnames(corr_filt) <- new_cols
colnames(pval_filt) <- new_cols

vmin <- vrange[1]; vmax <- vrange[2]
cmat <- pmax(pmin(as.matrix(corr_filt), vmax), vmin)
cmap <- colorRampPalette(hm_pal)(256)
nr <- nrow(cmat); nc <- ncol(cmat)
fig_w <- max(22, nc * 0.5 + 4)
fig_h <- max(10, nr * 0.4 + 4)

plot_hm <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(15, 1))
  par(mar = c(10, 7, 4, 1))

  val_sc <- (cmat - vmin) / (vmax - vmin)
  val_sc[val_sc < 0] <- 0; val_sc[val_sc > 1] <- 1
  ci <- pmax(1, pmin(256, round(val_sc * 255) + 1))
  cm <- matrix(cmap[ci], nrow = nr, ncol = nc)

  plot(NA, xlim = c(0, nc), ylim = c(nr, 0),
       xlab = "", ylab = "", axes = FALSE, asp = NA)

  for (i in 1:nr) for (j in 1:nc)
    rect(j - 1, i - 1, j, i, col = cm[i, j], border = NA)

  # sig --> asterisks
  for (i in 1:nr) for (j in 1:nc) {
    pv <- pval_filt[i, j]
    if (!is.na(pv) && pv < bh2) {
      tc <- if (abs(cmat[i, j]) > (vmax - vmin) * 0.3) "white" else "black"
      text(j - 0.5, i - 0.5, "*", cex = 1.2, col = tc, font = 2)
    }
  }

  if (!is.na(sep_idx)) abline(h = sep_idx - 1, col = "black", lwd = 2.5)

  axis(1, at = seq(0.5, nc - 0.5), labels = colnames(cmat),
       las = 2, cex.axis = 0.7, tick = FALSE, line = -0.5)
  axis(2, at = seq(0.5, nr - 0.5), labels = rownames(cmat),
       las = 1, cex.axis = 0.85, tick = FALSE, line = -0.5)
  mtext("Metabolites", side = 1, line = 8.5, cex = 1.2)
  mtext("Dietary Components", side = 2, line = 5.5, cex = 1.2)
  title(main = "Partial Spearman Correlation: Diet vs Metabolites\n(Adjusted for Age and Sex, * = BH-adjusted p < 0.2)",
        cex.main = 1.3, font.main = 2)

  if (!is.na(sep_idx)) {
    mtext("Nutrients", side = 4, at = (sep_idx - 1) / 2,
          las = 0, line = 0.3, cex = 0.9, font = 2)
    mtext("Food\nGroups", side = 4,
          at = sep_idx - 1 + (nr - sep_idx + 1) / 2,
          las = 0, line = 0.3, cex = 0.9, font = 2)
  }

  #color
  par(mar = c(10, 0.5, 4, 3))
  bv <- seq(vmin, vmax, length.out = 256)
  image(1, bv, matrix(bv, nrow = 1), col = cmap,
        xlab = "", ylab = "", axes = FALSE)
  axis(4, las = 1, cex.axis = 0.9)
  mtext("Partial Spearman Correlation (r)", side = 4, line = 2.5, cex = 0.9)
}

png(file.path(output_dir, "heatmap_larger_fonts_final.png"),
    width = fig_w, height = fig_h, units = "in", res = 300, bg = "white")
plot_hm(); dev.off()

pdf(file.path(output_dir, "heatmap_larger_fonts_final.pdf"),
    width = fig_w, height = fig_h)
plot_hm(); dev.off()
