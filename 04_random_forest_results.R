# ------------------------------------------------------------
# Script: 04_compute_RF_heatmap.R
# Purpose: Compute Random Forest model results for theta differences across soil layers
# Author: Haotian Wu
# Date: 14.10.2025
# ------------------------------------------------------------
# Required data:
#   - theta_diff_total_df.csv derived from '03_plot_Pearson_heatmap.R'

# Output:
#   - RF_single_results_df.csv : Contains Random Forest model performance metrics 
#                                (RMSE and R²) for each pore size range (r_min, r_max)
#   - RF_variable_importance_df.csv : Contains variable importance and contribution values 
#                                     for each predictor (OC, clay, silt, sand, climate_classes)
#                                     across all pore size ranges
# ------------------------------------------------------------

# Load required packages
install.packages("party")
install.packages("caret")

library(party)
library(caret) 
library(dplyr)


# Initialize result storage
RF_single_results_heatmap <- list()

# Get all unique r_min and r_max combinations
r_combinations <- unique(theta_diff_total_df[, c("r_min", "r_max")])

# Loop through each unique r_min and r_max pair
for (i in 1:nrow(r_combinations)) {
  
  r_min_val <- r_combinations$r_min[i]
  r_max_val <- r_combinations$r_max[i]
  
  # Subset data for the current combination
  subset_data <- theta_diff_total_df[theta_diff_total_df$r_min == r_min_val & theta_diff_total_df$r_max == r_max_val, ]
  
  # Define model formula
  my_formula <- as.formula(paste("theta_diff ~ oc + clay_tot_psa + silt_tot_psa + sand_tot_psa + climate_classes"))
  
  # Build Random Forest model
  set.seed(1234)
  RF_model <- cforest(my_formula, data = subset_data,
                      controls = cforest_unbiased(mtry = 2, ntree = 700, minbucket = 1, minsplit = 5))
  
  # Model evaluation
  fit_single <- predict(RF_model, newdata = NULL, OOB = TRUE)
  residuals_single <- subset_data$theta_diff - fit_single
  MSE_single <- sum(residuals_single^2) / length(subset_data$theta_diff)
  RMSE_single <- sqrt(MSE_single)
  R2_single <- 1 - MSE_single / var(subset_data$theta_diff)
  
  # Variable importance
  set.seed(1234)
  var_imp <- caret::varImp(RF_model)
  var_imp$variable <- as.factor(rownames(var_imp))
  var_imp <- var_imp[order(var_imp$Overall, decreasing = TRUE), ]
  
  # Relative importance
  var_imp <- var_imp %>%
    mutate(Overall = if_else(Overall < 0, 0, Overall)) %>%
    mutate(rel = Overall / sum(Overall) * 100)
  
  # Save results
  RF_single_results_heatmap[[paste0("r_min_", r_min_val, "_r_max_", r_max_val)]] <- list(
    RMSE = RMSE_single,
    R2 = R2_single,
    variable_importance = var_imp
  )
  
  # Clean temporary variables
  rm(fit_single, residuals_single, MSE_single, RMSE_single, R2_single, var_imp, subset_data, RF_model)
}

# Combine model summary results
RF_single_results_df <- do.call(rbind, lapply(names(RF_single_results_heatmap), function(name) {
  data.frame(
    r_min = as.numeric(sub("r_min_", "", strsplit(name, "_r_max_")[[1]][1])),
    r_max = as.numeric(sub("r_max_", "", strsplit(name, "_r_max_")[[1]][2])),
    RMSE = RF_single_results_heatmap[[name]]$RMSE,
    R2 = RF_single_results_heatmap[[name]]$R2,
    stringsAsFactors = FALSE
  )
}))

# Combine variable importance results
RF_variable_importance_df <- do.call(rbind, lapply(names(RF_single_results_heatmap), function(name) {
  var_imp <- RF_single_results_heatmap[[name]]$variable_importance
  var_imp$r_min <- as.numeric(sub("r_min_", "", strsplit(name, "_r_max_")[[1]][1]))
  var_imp$r_max <- as.numeric(sub("r_max_", "", strsplit(name, "_r_max_")[[1]][2]))
  var_imp
}))

# Calculate contribution of each feature
RF_variable_importance_df <- RF_variable_importance_df %>%
  left_join(RF_single_results_df, by = c("r_min", "r_max")) %>%
  mutate(
    variable = as.character(variable),
    contribution = as.numeric(rel * R2)
  )

# Export dataset
write.csv(RF_single_results_df, "data/RF_single_results_df.csv", row.names = FALSE)
write.csv(RF_variable_importance_df, "data/RF_variable_importance_df.csv", row.names = FALSE)




