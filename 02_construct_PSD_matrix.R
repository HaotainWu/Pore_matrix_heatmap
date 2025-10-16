# ------------------------------------------------------------
# Script: 02_build_pore_matrix.R
# Purpose: conduct pore size range matrix for each sample
# Author: Haotian Wu
# Date: 14.10.2025
# ------------------------------------------------------------
# Required data:
  #   - Export_data.csv derived from '01_fit_vanGenuchten.R'

# Output:
  #   - theta_diff_total_df.csv : combined results including layer_id, r_min, r_max, theta_diff
# ------------------------------------------------------------

# Load required packages
library(dplyr)

# load dataset
Export_data <- read.csv("Data/Export_data.csv")

# ---- Step 1: Define matric potential sequence ----
h_map <- 10^seq(-1, 3, by = 0.02)  # h in meters

# ---- Step 2: Create all possible combinations of h_min and h_max ----
h_combinations <- expand.grid(h_min = h_map, h_max = h_map)

# Calculate corresponding pore radii (r_min, r_max) 
h_combinations$r_min <- 30 / h_combinations$h_max
h_combinations$r_max <- 30 / h_combinations$h_min

# Keep only valid combinations where r_max > r_min
h_valid <- h_combinations[h_combinations$r_max > h_combinations$r_min, ]
h_valid$r_min <- round(h_valid$r_min, 7)
h_valid$r_max <- round(h_valid$r_max, 7)

# ---- Step 3: Initialize list to store results for each sample ----
theta_diff_total_list <- list()

# ---- Step 4: Calculate theta differences for each layer_id ----
for (layer in unique(Export_data$layer_id)) {
  
  # Filter data for current layer
  layer_data <- Export_data %>% filter(layer_id == layer)
  
  # Initialize list to store theta differences for this layer
  theta_diff_list <- list()
  
  # Loop over all valid pore size combinations
  for (i in 1:nrow(h_valid)) {
    r_min <- h_valid$r_min[i]
    r_max <- h_valid$r_max[i] 
    
    # Extract theta values corresponding to r_min and r_max
    theta_r_min <- layer_data$theta[layer_data$r_um == r_min]
    theta_r_max <- layer_data$theta[layer_data$r_um == r_max]  
    
    # Compute theta difference if both radii exist
    if (length(theta_r_min) > 0 && length(theta_r_max) > 0) {
      theta_diff <- theta_r_max - theta_r_min
    } else {
      theta_diff <- NA  # Assign NA if either radius is missing
    }
    
    # If theta_diff is a vector, take the first value
    if (length(theta_diff) > 1) {
      theta_diff <- theta_diff[1]
    }
    
    # Store results for this combination
    theta_diff_list[[i]] <- data.frame(
      layer_id = layer,
      h_min = h_valid$h_min[i],
      h_max = h_valid$h_max[i],
      r_min = r_min,
      r_max = r_max,
      theta_diff = theta_diff
    )
  }
  
  # Combine all combinations for each sample and add to total list
  theta_diff_total_list[[layer]] <- do.call(rbind, theta_diff_list)
}

# ---- Step 5: Combine results from all samples ----
theta_diff_total_df <- do.call(rbind, theta_diff_total_list)

# ---- Save results ----
write.csv(
  theta_diff_total_df,
  "Data/theta_diff_total_df.csv",
  row.names = FALSE
)


