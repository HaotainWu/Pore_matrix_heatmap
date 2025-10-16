# ------------------------------------------------------------
# Script: 03_plot_Pearson_heatmap.R
# Purpose: Calculate Pearson correlations between soil properties and theta differences for each pore size range
# Author: Haotian Wu
# Date: 14.10.2025
# ------------------------------------------------------------
# Required data:
#   - theta_diff_total_df.csv derived from '02_build_pore_matrix'
#   - soil_properties.csv (inclduing SOC content, and clay/silt/sand content,and climate information)

# Output:
#   - pearson_corr.csv : Pearson's correlation results between soil properties and pore size range matrix
# ------------------------------------------------------------

# Load required packages
library(dplyr)
library(ppcor)

# load soil properties dataset
theta_diff_total_df <- read.csv("Data/theta_diff_total_df.csv")
Soil_properties <- read.csv("Data/Soil_properties.csv")

# Merge dataset
Soil_properties$layer_id <- as.character(Soil_properties$layer_id)
theta_diff_total_df <- left_join(Soil_properties, theta_diff_total_df, by = "layer_id")

# Grouping the data by each pair of r_min and r_max and calculating correlation coefficients and p-values
pearson_corr <- theta_diff_total_df %>%
  group_by(r_min, r_max) %>%
  filter(n() > 1) %>%  # Filter out groups with only one row
  summarise(
    oc_corr = cor(oc, theta_diff),
    oc_p_value = cor.test(oc, theta_diff)$p.value,
    oc_lm_slope = coef(lm(theta_diff ~ oc))[2],
    sand_tot_psa_corr = cor(sand_tot_psa, theta_diff),
    sand_tot_psa_p_value = cor.test(sand_tot_psa, theta_diff)$p.value,
    sand_lm_slope = coef(lm(theta_diff ~ sand_tot_psa))[2],
    clay_tot_psa_corr = cor(clay_tot_psa, theta_diff),
    clay_tot_psa_p_value = cor.test(clay_tot_psa, theta_diff)$p.value,
    clay_lm_slope = coef(lm(theta_diff ~ clay_tot_psa))[2],
    silt_tot_psa_corr = cor(silt_tot_psa, theta_diff),
    silt_tot_psa_p_value = cor.test(silt_tot_psa, theta_diff)$p.value,
    silt_lm_slope = coef(lm(theta_diff ~ silt_tot_psa))[2],
    .groups = "drop"
  )

# Export dataset
write.csv(
  pearson_corr,
  "Data/pearson_corr.csv",
  row.names = FALSE
)
