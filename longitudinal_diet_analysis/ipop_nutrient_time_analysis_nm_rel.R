library(dplyr); library(tidyverse); library(lme4); library(lmerTest); library(emmeans); library(broom.mixed)

data_dir <- "../../Diet_Data/Clean_Data"
setwd(data_dir)

df_a_nutrients <- read.csv("df_a_nutrients.csv", stringsAsFactors = F)
df_a_demographics <- read.csv("df_a_demographics.csv", stringsAsFactors = F)

dim(df_a_nutrients) #151 64
dim(df_a_demographics) #151 25
df_a_demographics$X <- NULL

intersect(names(df_a_nutrients), names(df_a_demographics)) #SubjectID, Time

df_nutrients_a <- df_a_nutrients %>%
  left_join(df_a_demographics, by = c("SubjectID", "Time"))

dim(df_nutrients_a) #151 86
names(df_nutrients_a)

######################################
df_nutrients_a <- df_nutrients_a %>%
  mutate(
    Time = factor(Time, levels = c("T1","T2","T3","T4")),
    sspg_status = factor(sspg_status, levels = c("IS","IR")))

str(df_nutrients_a)

nutrient_cols <- c("Energy", "percent_carb_cal", "percent_pro_cal", "percent_fat_cal", "FatCals", "SatCals", "Pro", "Carb", "TotFib", "TotSolFib", "Fib", "SolFib",
                    "Sugar", "SugAdd", "MonSac", "Disacc", "OCarb", "Fat", "SatFat", "MonoFat", "PolyFat", "TransFat", "Chol", "Water", "VitA_IU", "VitA_RAE",
                    "Caroten", "Retinol", "BetaCaro", "VitB1", "VitB2", "VitB3", "VitB_NE", "VitB6", "VitB12", "Biot", "VitC", "VitD_IU", "VitD_mcg",
                    "VitE_a_Toco", "Folate", "Fol_DFE", "VitK", "Panto", "Ca", "Chrom", "Cu", "Fluor", "Iodine", "Fe", "Mg", "Mn", "Moly", "Phos", "K",
                    "Se", "Na", "Zn", "W3", "W6", "Caff", "Chln")

setdiff(nutrient_cols, names(df_nutrients_a))

##############################################################
### Mixed models per nutrient (Cases 1-3 + plotting means) ###
##############################################################
run_mixed_for_nutrient <- function(nut_name, dat) {

  dat_all <- dat %>% dplyr::filter(!is.na(.data[[nut_name]]))

  sample_size_table <- dat_all %>%
    dplyr::group_by(Time) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_IS    = sum(sspg_status == "IS", na.rm = TRUE),
      n_IR    = sum(sspg_status == "IR", na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(nutrient = nut_name, .before = 1)
  print(sample_size_table)

  #### 1. OVERALL MODEL (Case 1)

  f_overall <- as.formula(paste0("log1p(", nut_name, ") ~ Time + (1|SubjectID)"))
  fit_overall <- lme4::lmer(f_overall, data = dat_all, REML = FALSE)
  emm_overall <- emmeans::emmeans(fit_overall, ~ Time)
  cont_overall <- emmeans::contrast(emm_overall,
                                    method = "trt.vs.ctrl",
                                    ref = "T1") %>%
    summary(infer = c(TRUE, TRUE)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      nutrient = nut_name,
      group    = "Overall",
      .before  = 1
    ) %>%
    dplyr::rename(
      estimate = estimate,
      se       = SE,
      df       = df,
      t        = t.ratio,
      p_raw    = p.value
    ) %>%
    dplyr::select(nutrient, group, contrast, estimate, se, df, t, p_raw)

  means_overall <- summary(emm_overall, infer = c(TRUE, TRUE)) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      nutrient = nut_name,
      Time,
      group    = "Overall",
      mean     = expm1(emmean),
      lower_ci = expm1(lower.CL),
      upper_ci = expm1(upper.CL)
    )



  #### 2. IS & IR MODEL (Case 2 & 3)
  dat_status <- dat_all %>%
    dplyr::filter(!is.na(sspg_status),
                  sspg_status %in% c("IS","IR"))
  if (nrow(dat_status) > 0 && dplyr::n_distinct(dat_status$sspg_status) == 2) {

    f_group <- as.formula(paste0("log1p(", nut_name, ") ~ Time * sspg_status + (1|SubjectID)"))
    fit_group <- lme4::lmer(f_group, data = dat_status, REML = FALSE)

    emm_group <- emmeans::emmeans(fit_group, ~ Time | sspg_status)
    cont_group <- emmeans::contrast(emm_group,
                                    method = "trt.vs.ctrl",
                                    ref = "T1") %>%
      summary(infer = c(TRUE, TRUE)) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        nutrient = nut_name,
        group    = as.character(sspg_status),
        .before  = 1
      ) %>%
      dplyr::rename(
        estimate = estimate,
        se       = SE,
        df       = df,
        t        = t.ratio,
        p_raw    = p.value
      ) %>%
      dplyr::select(nutrient, group, contrast, estimate, se, df, t, p_raw)

    means_group <- summary(emm_group, infer = c(TRUE, TRUE)) %>%
      tibble::as_tibble() %>%
      dplyr::transmute(
        nutrient = nut_name,
        Time,
        group    = as.character(sspg_status),
        mean     = expm1(emmean),
        lower_ci = expm1(lower.CL),
        upper_ci = expm1(upper.CL)
      )

  } else {
    cont_group  <- tibble::tibble()
    means_group <- tibble::tibble()
  }

  # Combine
  contrasts_tbl <- dplyr::bind_rows(cont_overall, cont_group)
  means_tbl     <- dplyr::bind_rows(means_overall, means_group)

  list(
    contrasts   = contrasts_tbl,
    means       = means_tbl,
    sample_size = sample_size_table
  )
}

res_list <- purrr::map(
  nutrient_cols,
  run_mixed_for_nutrient,
  dat = df_nutrients_a
)

mixed_contrasts_raw <- purrr::map_dfr(res_list, "contrasts") %>%
  dplyr::mutate(
    percent_change = (exp(estimate) - 1) * 100) %>%
  dplyr::group_by(nutrient, group) %>%
  dplyr::mutate(p_BH = p.adjust(p_raw, method = "BH")) %>%
  dplyr::ungroup()

plot_means_all <- purrr::map_dfr(res_list, "means") %>%
  dplyr::mutate(group = factor(group, levels = c("Overall","IS","IR")))

sample_sizes_all <- purrr::map_dfr(res_list, "sample_size")

head(mixed_contrasts_raw)
head(plot_means_all)

#################################################################################
##################### Plot Comparisons (All participants) ####################
################################################################################

library(ggplot2)
p_all <- ggplot(plot_means_all,
                aes(x = Time, y = mean,
                    group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(values = c("Overall" = "gray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(x = "Time", y = "Value", color = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

p_all

if (!dir.exists("Figures")) {
  dir.create("Figures")
}

ggplot2::ggsave(
  filename = "Figures/plot_all_participants.png", plot = p_all, width = 14, height = 10, dpi= 300)

ggplot2::ggsave(filename = "Figures/p_all.pdf", plot = p_all, width = 14, height = 10)

## Mann–Whitney tests per nutrient per time (IS vs IR)
run_mw_for_nutrient <- function(nut_name, dat) {
  dat_n <- dat %>%
    filter(!is.na(.data[[nut_name]]),
           !is.na(sspg_status),
           sspg_status %in% c("IS","IR"))

  map_dfr(levels(dat_n$Time), function(tlev) {
    subdat <- dat_n %>% filter(Time == tlev)

    if (n_distinct(subdat$sspg_status) < 2) {
      tibble(
        nutrient = nut_name,
        Time     = tlev,
        n_IS     = sum(subdat$sspg_status == "IS"),
        n_IR     = sum(subdat$sspg_status == "IR"),
        p_raw    = NA_real_,
        W        = NA_real_
      )
    } else {
      w <- wilcox.test(
        subdat[[nut_name]] ~ subdat$sspg_status,
        exact = FALSE
      )
      tibble(
        nutrient = nut_name,
        Time     = tlev,
        n_IS     = sum(subdat$sspg_status == "IS"),
        n_IR     = sum(subdat$sspg_status == "IR"),
        p_raw    = w$p.value,
        W        = w$statistic
      )
    }
  })
}

mw_results_raw <- map_dfr(nutrient_cols, run_mw_for_nutrient, dat = df_nutrients_a)

# BH adjust within each nutrient over the four time points
mw_results <- mw_results_raw %>%
  group_by(nutrient) %>%
  mutate(p_BH = p.adjust(p_raw, method = "BH")) %>%
  ungroup()

write.csv(mw_results, file = "mann_whitney_results.csv", row.names = FALSE)

mw_results %>%
  group_by(nutrient, Time) %>%
  summarise(n_IS = unique(n_IS), n_IR = unique(n_IR), .groups = "drop") %>%
  arrange(n_IS + n_IR)

#################
## Significant nutrients by RAW p val < 0.05

# Case 1–3 (Overall, IS, IR)
sig_mixed_raw <- mixed_contrasts_raw %>%
  dplyr::filter(p_raw < 0.05) %>%
  dplyr::pull(nutrient) %>%
  unique()

# Case 4 (IS vs IR at any time)
sig_mw_raw <- mw_results %>%
  dplyr::filter(p_raw < 0.05) %>%
  dplyr::pull(nutrient) %>%
  unique()

# Combine (at least one significant effect)
sig_nutrients_raw <- union(sig_mixed_raw, sig_mw_raw)

length(sig_nutrients_raw)
sig_nutrients_raw

#####
plot_means_sig_raw <- plot_means_all %>%
  dplyr::filter(nutrient %in% sig_nutrients_raw)

p_sig_raw <- ggplot(plot_means_sig_raw,
                    aes(x = Time, y = mean,
                        group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(values = c("Overall" = "gray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(x = "Time", y = "Value", color = NULL,
       title = "Nutrients with ≥1 significant effect (raw p < 0.05)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

p_sig_raw

#####
## Significant nutrients by BH-adjusted p < 0.20

# Case 1–3
sig_mixed_BH <- mixed_contrasts_raw %>%
  dplyr::filter(p_BH < 0.20) %>%
  dplyr::pull(nutrient) %>%
  unique()

# Case 4
sig_mw_BH <- mw_results %>%
  dplyr::filter(p_BH < 0.20) %>%
  dplyr::pull(nutrient) %>%
  unique()

sig_nutrients_BH <- union(sig_mixed_BH, sig_mw_BH)
sig_mixed_BH

length(sig_nutrients_BH) #27 nutrients
sig_nutrients_BH

plot_means_sig_BH <- plot_means_all %>%
  dplyr::filter(nutrient %in% sig_nutrients_BH)

p_sig_BH <- ggplot(plot_means_sig_BH,
                   aes(x = Time, y = mean,
                       group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(values = c("Overall" = "darkgray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(x = "Time", y = "Value", color = NULL,
       title = "Nutrients with ≥1 effect (BH-adjusted p < 0.20)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

p_sig_BH

if (!dir.exists("Figures")) {
  dir.create("Figures")
}

ggplot2::ggsave(
  filename = "Figures/p_sig_BH.png",
  plot     = p_sig_BH,
  width    = 14,
  height   = 10,
  dpi      = 300
)

ggplot2::ggsave(
  filename = "Figures/p_sig_BH.pdf",
  plot     = p_sig_BH,
  width    = 14,
  height   = 10
)

###############################################
## Plot v3: selected nutrients in given order
###############################################

nutrients_v3 <- c(
  "Fib", "SolFib", "Sugar", "Fat", "Chol", "Na",
  "VitA_IU", "Caroten", "VitB_NE", "VitC",
  "VitE_a_Toco", "VitK", "Chln", "Biot",
  "K", "Iodine")

nutrient_labels <- c(
  Fib          = "Fiber",
  SolFib       = "Soluble Fiber",
  Sugar        = "Total Sugars",
  Fat          = "Fat",
  Chol         = "Cholesterol",
  Na           = "Sodium",
  VitA_IU      = "Vitamin A (IU)",
  Caroten      = "Carotene",
  VitB_NE      = "Vitamin B3",
  VitC         = "Vitamin C",
  VitE_a_Toco  = "Vitamin E",
  VitK         = "Vitamin K",
  Chln         = "Choline",
  Biot         = "Biotin",
  K            = "Potassium",
  Iodine       = "Iodine"
)

plot_means_v3 <- plot_means_all %>%
  dplyr::filter(nutrient %in% nutrients_v3) %>%
  dplyr::mutate(
    nutrient       = factor(nutrient, levels = nutrients_v3),
    nutrient_label = factor(nutrient_labels[nutrient],
                            levels = unname(nutrient_labels))
  )


p_v3 <- ggplot(plot_means_v3,
               aes(x = Time, y = mean,
                   color = group, group = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient_label, scales = "free_y",
             nrow = 4, ncol = 4) +
  scale_color_manual(values = c("Overall" = "darkgray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(
    x = "Time",
    y = "Dietary Intake",
    color = NULL,
    title = " "
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold", size = 17),   # facet banner
    axis.title.x     = element_text(size = 16),                  # x-axis title
    axis.title.y     = element_text(size = 16),                  # y-axis title
    axis.text.x      = element_text(size = 14),                  # x tick labels
    axis.text.y      = element_text(size = 14),                  # y tick labels
    legend.text      = element_text(size = 15),                  # legend labels
    legend.title     = element_text(size = 17),                  # legend title (if any)
    legend.position  = "right"
  )

p_v3

if (!dir.exists("Figures")) dir.create("Figures")

ggplot2::ggsave("Figures/plot_v4_selected_nutrients.png",
                p_v3, width = 14, height = 12, dpi = 300)

ggplot2::ggsave("Figures/plot_v4_selected_nutrients.pdf",
                p_v3, width = 14, height = 12)

# high Vitamin A intakes
vit_a <- df_nutrients_a %>%
  dplyr::filter(VitA_IU > 10000) %>%
  dplyr::select(SubjectID, Time, VitA_IU, dplyr::everything())



###############################################################
####################### Food Group ############################
###############################################################

df_a_food_energy <- read.csv("df_a_energy.csv", stringsAsFactors = F)
df_a_demographics <- read.csv("df_a_demographics.csv", stringsAsFactors = F)

dim(df_a_food_energy) #151 22
dim(df_a_demographics) #151 25
df_a_demographics$X <- NULL

names(df_a_food_energy)

intersect(names(df_a_food_energy), names(df_a_demographics)) #SubjectID, Time

df_a_food_energy <- df_a_food_energy %>%
  left_join(df_a_demographics, by = c("SubjectID", "Time"))

dim(df_a_food_energy) #151 44
names(df_a_food_energy)

######################################
df_a_food_energy <- df_a_food_energy %>%
  mutate(
    Time = factor(Time, levels = c("T1","T2","T3","T4")),
    sspg_status = factor(sspg_status, levels = c("IS","IR")))

str(df_a_food_energy)
#paste(names(df_a_food_energy),collapse=", ")

food_cols <- c("Bread_Products", "Coffee", "Dairy_Products", "Fats.and.Oils", "Fruits", "Meat_Poultry", "Non_Starchy_Vegetables", "Nuts_Seeds", "Salad",
                    "Starchy_Vegetables", "Sweets_Sweeteners", "Cereals_Grains", "Sauces..Spices..and.Herbs", "Snacks", "Seafood", "Tea", "Eggs", "Legume_Products",
                    "Mixed_Meal", "Rice_Pasta_Noodles")

setdiff(food_cols, names(df_a_food_energy))

##############################################################
### Mixed models per food group ##############################
##############################################################

res_list <- purrr::map(
  food_cols,
  run_mixed_for_nutrient,
  dat = df_a_food_energy)

mixed_contrasts_raw <- purrr::map_dfr(res_list, "contrasts") %>%
  dplyr::mutate(
    percent_change = (exp(estimate) - 1) * 100
  ) %>%
  dplyr::group_by(nutrient, group) %>%
  dplyr::mutate(p_BH = p.adjust(p_raw, method = "BH")) %>%
  dplyr::ungroup()

plot_means_all <- purrr::map_dfr(res_list, "means") %>%
  dplyr::mutate(group = factor(group, levels = c("Overall","IS","IR")))

sample_sizes_all <- purrr::map_dfr(res_list, "sample_size")

head(mixed_contrasts_raw)
head(plot_means_all)

#################################################################################
##################### Plot Comparisons (All participants) ######################
################################################################################

library(ggplot2)
p_all <- ggplot(plot_means_all,
                aes(x = Time, y = mean,
                    group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(values = c("Overall" = "gray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(x = "Time", y = "Value", color = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

p_all

if (!dir.exists("Figures")) {
  dir.create("Figures")
}

ggplot2::ggsave(
  filename = "Figures/plot_all_participants_food_energy.png", plot = p_all, width = 14, height = 10, dpi= 300)

ggplot2::ggsave(filename = "Figures/p_all_food_energy.pdf", plot = p_all, width = 14, height = 10)

## Mann–Whitney tests per nutrient per time (IS vs IR)
########################################################

run_mw_for_nutrient <- function(nut_name, dat) {
  dat_n <- dat %>%
    filter(!is.na(.data[[nut_name]]),
           !is.na(sspg_status),
           sspg_status %in% c("IS","IR"))

  map_dfr(levels(dat_n$Time), function(tlev) {
    subdat <- dat_n %>% filter(Time == tlev)

    if (n_distinct(subdat$sspg_status) < 2) {
      tibble(
        nutrient = nut_name,
        Time     = tlev,
        n_IS     = sum(subdat$sspg_status == "IS"),
        n_IR     = sum(subdat$sspg_status == "IR"),
        p_raw    = NA_real_,
        W        = NA_real_
      )
    } else {
      w <- wilcox.test(
        subdat[[nut_name]] ~ subdat$sspg_status,
        exact = FALSE
      )
      tibble(
        nutrient = nut_name,
        Time     = tlev,
        n_IS     = sum(subdat$sspg_status == "IS"),
        n_IR     = sum(subdat$sspg_status == "IR"),
        p_raw    = w$p.value,
        W        = w$statistic
      )
    }
  })
}

mw_results_raw <- map_dfr(food_cols, run_mw_for_nutrient, dat = df_a_food_energy)

#BH adjust
mw_results <- mw_results_raw %>%
  group_by(nutrient) %>%
  mutate(p_BH = p.adjust(p_raw, method = "BH")) %>%
  ungroup()

write.csv(mw_results, file = "mann_whitney_results_food_energy.csv", row.names = FALSE)

mw_results %>%
  group_by(nutrient, Time) %>%
  summarise(n_IS = unique(n_IS), n_IR = unique(n_IR), .groups = "drop") %>%
  arrange(n_IS + n_IR)


#################
## ---- Significant nutrients by RAW p val < 0.05 ----

# Case 1–3 (Overall, IS, IR)
sig_mixed_raw <- mixed_contrasts_raw %>%
  dplyr::filter(p_raw < 0.05) %>%
  dplyr::pull(nutrient) %>%
  unique()

# Case 4 (IS vs IR at any time)
sig_mw_raw <- mw_results %>%
  dplyr::filter(p_raw < 0.05) %>%
  dplyr::pull(nutrient) %>%
  unique()

# Combine: at least one significant effect
sig_nutrients_raw <- union(sig_mixed_raw, sig_mw_raw)

length(sig_nutrients_raw)
sig_nutrients_raw
#####
plot_means_sig_raw <- plot_means_all %>%
  dplyr::filter(nutrient %in% sig_nutrients_raw)

p_sig_raw <- ggplot(plot_means_sig_raw,
                    aes(x = Time, y = mean,
                        group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(values = c("Overall" = "gray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(x = "Time", y = "Value", color = NULL,
       title = "Nutrients with ≥1 significant effect (raw p < 0.05)") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

p_sig_raw

#####
## ---- Significant nutrients by BH-adjusted p < 0.20 ----

# Case 1–3
sig_mixed_BH <- mixed_contrasts_raw %>%
  dplyr::filter(p_BH < 0.20) %>%
  dplyr::pull(nutrient) %>%
  unique()

# Case 4
sig_mw_BH <- mw_results %>%
  dplyr::filter(p_BH < 0.20) %>%
  dplyr::pull(nutrient) %>%
  unique()

sig_nutrients_BH <- union(sig_mixed_BH, sig_mw_BH)
sig_mixed_BH

length(sig_nutrients_BH)
sig_nutrients_BH

library(forcats)

plot_means_sig_BH <- plot_means_all %>%
  dplyr::filter(nutrient %in% sig_nutrients_BH) %>%
  dplyr::mutate(
    nutrient = fct_recode(
      nutrient,
      "Sauces_Spices_Herbs" = "Sauces..Spices..and.Herbs",
      "Pizza_Sandwiches_Wraps" = "Mixed_Meal"
    )
  )

p_sig_BH <- ggplot(plot_means_sig_BH,
                   aes(x = Time, y = mean,
                       group = group, color = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.1) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(values = c("Overall" = "darkgray",
                                "IS"      = "#298C8C",
                                "IR"      = "#F1A226")) +
  labs(x = "Time", y = "Value", color = NULL,
       title = "") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

p_sig_BH

if (!dir.exists("Figures")) {
  dir.create("Figures")
}

ggplot2::ggsave(
  filename = "Figures/p_sig_BH_food_energy.png",
  plot     = p_sig_BH,
  width    = 14,
  height   = 10,
  dpi      = 300
)

ggplot2::ggsave(
  filename = "Figures/p_sig_BH_food_energy.pdf",
  plot     = p_sig_BH,
  width    = 14,
  height   = 10
)


############### Dietary Pattern Clusters ############################
#####################################################################

df_b_nutrients <- read.csv("df_b_nutrients.csv", stringsAsFactors = F)
df_b_demographics <- read.csv("df_b_demographics.csv", stringsAsFactors = F)
df_cluster <- read.csv("diet_clusters_ae801.csv", stringsAsFactors = F)

dim(df_b_nutrients) #71 63
dim(df_b_demographics) #71 24
dim(df_cluster) #71 2
df_b_demographics$X <- NULL

intersect(names(df_b_nutrients), names(df_b_demographics)) #SubjectID
df_nutrients_b <- df_b_nutrients %>%
  left_join(df_b_demographics, by = c("SubjectID"))

intersect(names(df_nutrients_b), names(df_cluster)) #SubjectID

df_nutrients_b <- df_nutrients_b %>%
  left_join(df_cluster, by = c("SubjectID"))

dim(df_nutrients_b) #71 86
names(df_nutrients_b)

df_b_food_energy <- read.csv("df_b_energy.csv", stringsAsFactors = F)

dim(df_b_food_energy) #71 21

intersect(names(df_b_food_energy), names(df_nutrients_b)) #SubjectID

df_diet_b <- df_nutrients_b %>%
  left_join(df_b_food_energy, by = c("SubjectID"))

dim(df_diet_b) #71 106
names(df_diet_b)
str(df_diet_b)

######################################
df_diet_b <- df_diet_b %>%
  mutate(
    diet_cluster = factor(diet_cluster, levels = c("0", "1")),
    sspg_status = factor(sspg_status, levels = c("IS","IR")))

str(df_diet_b)
#paste(names(df_nutrients_a),collapse=", ")

diet_cols <- c("Energy", "percent_carb_cal", "percent_pro_cal", "percent_fat_cal", "FatCals", "SatCals", "Pro", "Carb", "TotFib", "TotSolFib", "Fib", "SolFib",
                   "Sugar", "SugAdd", "MonSac", "Disacc", "OCarb", "Fat", "SatFat", "MonoFat", "PolyFat", "TransFat", "Chol", "Water", "VitA_IU", "VitA_RAE",
                   "Caroten", "Retinol", "BetaCaro", "VitB1", "VitB2", "VitB3", "VitB_NE", "VitB6", "VitB12", "Biot", "VitC", "VitD_IU", "VitD_mcg",
                   "VitE_a_Toco", "Folate", "Fol_DFE", "VitK", "Panto", "Ca", "Chrom", "Cu", "Fluor", "Iodine", "Fe", "Mg", "Mn", "Moly", "Phos", "K",
                   "Se", "Na", "Zn", "W3", "W6", "Caff", "Chln",
                   "Bread_Products", "Coffee", "Dairy_Products", "Fats.and.Oils", "Fruits", "Meat_Poultry", "Non_Starchy_Vegetables", "Nuts_Seeds", "Salad",
                   "Starchy_Vegetables", "Sweets_Sweeteners", "Cereals_Grains", "Sauces..Spices..and.Herbs", "Snacks", "Seafood", "Tea", "Eggs", "Legume_Products",
                   "Mixed_Meal", "Rice_Pasta_Noodles")

setdiff(diet_cols, names(df_diet_b))

run_mw_for_dietcluster <- function(var_name, dat) {
  dat_n <- dat %>%
    filter(
      !is.na(.data[[var_name]]),
      !is.na(diet_cluster),
      diet_cluster %in% c("0", "1")
    )

  if (n_distinct(dat_n$diet_cluster) < 2L) {
    return(
      tibble(
        variable  = var_name,
        n_0       = sum(dat_n$diet_cluster == "0"),
        n_1       = sum(dat_n$diet_cluster == "1"),
        median_0  = NA_real_,
        median_1  = NA_real_,
        p_raw     = NA_real_,
        W         = NA_real_
      )
    )
  }

  w <- wilcox.test(
    dat_n[[var_name]] ~ dat_n$diet_cluster,
    exact = FALSE
  )

  tibble(
    variable  = var_name,
    n_0       = sum(dat_n$diet_cluster == "0"),
    n_1       = sum(dat_n$diet_cluster == "1"),
    median_0  = median(dat_n[[var_name]][dat_n$diet_cluster == "0"], na.rm = TRUE),
    median_1  = median(dat_n[[var_name]][dat_n$diet_cluster == "1"], na.rm = TRUE),
    p_raw     = w$p.value,
    W         = unname(w$statistic)
  )
}

diet_mw_results_raw <- map_dfr(diet_cols, run_mw_for_dietcluster, dat = df_diet_b)

diet_mw_results <- diet_mw_results_raw %>%
  mutate(p_BH = p.adjust(p_raw, method = "BH")) %>%
  arrange(p_BH)

write.csv(
  diet_mw_results,
  file = "mann_whitney_results_diet_cluster_temp.csv",
  row.names = FALSE
)

sessionInfo()
