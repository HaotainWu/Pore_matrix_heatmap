# ------------------------------------------------------------
# Script: 01_fit_vanGenuchten.R
# Purpose: Fit van Genuchten water retention model for each soil layer and generate pore size distribution curve
# Author: Haotian Wu
# Date: 14.10.2025
# ------------------------------------------------------------
# Description:
# This script performs the following tasks:
# 1. Loads the water retention dataset (Data_wrc.csv)
# 2. Estimates van Genuchten parameters for each soil layer using 
#    'soilhypfit::fit_wrc_hcc()' (maximum likelihood fitting)
# 3. Merges fitted parameters with original water retention data
# 4. Computes predicted water retention curves (theta vs h) for each sample
# 5. Evaluates model performance (R² and RMSE)
# 6. Converts matric potential (h) to equivalent pore radius (r, µm)
# 7. Exports combined results to CSV for further analysis

# Required data:
#   - Data_wrc: dataframe containing columns
#       layer_id, lab_head_m (lab measured matric potential, m), lab_wrc (lab measured volumetric water content, m3/m3)
#
# Output:
#   - Export_data.csv : combined results including layer_id, h, theta, R², RMSE, and r (µm)
# ------------------------------------------------------------
# Load required packages 
install.packages("soilhypfit")
library(soilhypfit)
library(dplyr)

#Load data
Data_wrc <- read.csv("Data/Data_wrc.csv")

# calculate parameters of VG models using 'soilhypfit' package
system.time(r1 <- fit_wrc_hcc(
  wrc_formula = lab_wrc ~ lab_head_m | layer_id,
  data = Data_wrc,
  control = control_fit_wrc_hcc(
    min_nobs_wc = 4,
    keep_empty_fits = FALSE,
    nloptr = control_nloptr(maxeval = 250),
    param_bound = param_boundf(
      alpha = c(1.490116e-07 , 100.),
      n = c(1., 7.)
    )
  )
))

Parameters<- coef(r1,gof = TRUE)

Parameters <- cbind(layer_id = rownames(Parameters), Parameters)
rownames(Parameters) <- NULL

# combine parameters data and original WRC data
Parameters$layer_id <- as.character(Parameters$layer_id)
Data_wrc$layer_id <- as.character(Data_wrc$layer_id)
Merged_data <- left_join(Data_wrc, Parameters, by = "layer_id")

# Initialize empty dataframe for storing results 
Export_data <- data.frame()

# Define range of matric potentials (h, in m)
h <- 10^seq(-2, 3, by = 0.02)

# Fit van Genuchten model for each sample 
for (id in unique(Merged_data$layer_id)) {
  
  # Subset data for the current layer_id
  params <- subset(Merged_data, layer_id == id)
  
  # Calculate predicted theta values over defined h range
  theta <- wc_model(
    h,
    nlp = c(alpha = params$alpha[1], n = params$n[1]),
    lp  = c(thetar = params$thetar[1], thetas = params$thetas[1])
  )
  
  # Model performance evaluation 
  # Compute predicted theta for observed heads
  actual_theta    <- params$lab_wrc
  predicted_theta <- wc_model(
    params$lab_head_m,
    nlp = c(alpha = params$alpha[1], n = params$n[1]),
    lp  = c(thetar = params$thetar[1], thetas = params$thetas[1])
  )
  
  # Calculate R² and RMSE
  ss_res    <- sum((actual_theta - predicted_theta)^2)
  ss_tot    <- sum((actual_theta - mean(actual_theta))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  rmse      <- sqrt(mean((actual_theta - predicted_theta)^2))
  
  # Combine results for export 
  df_plot <- data.frame(
    h         = h,
    theta     = theta,
    layer_id  = rep(id, length(h)),
    r_squared = r_squared,
    rmse      = rmse
  )
  
  # Convert matric potential (h) to equivalent pore radius (r, µm)
  df_plot$r_um <- round(30 / df_plot$h,7)
  
  # Append to export dataframe
  Export_data <- rbind(Export_data, df_plot)
}

# ---- Save results ----
write.csv(
  Export_data,
  theta_diff_total_df,
  "Data/Export_data.csv",
  row.names = FALSE
)
