# Pore_matrix_heatmap

This repository contains R scripts for analyzing soil pore size distributions, 
calculating correlations with soil properties, and evaluating variable importance using Random Forest models. 
The scripts cover the full workflow from fitting water retention curves to generating heatmap results for the pore size range matrix.

---

## Repository Structure
``` 
├── data/                          # Folder containing input and output data files
│   ├── Data_wrc.csv               # Measured water retention data
│   ├── soil_properties.csv        # Soil texture and organic carbon data
│   ├── Export_data.csv            # Fitted van Genuchten model output (θ–h and θ–r)
│   ├── Pearson_corr.csv           # Pearson's correlation coefficients results
├── scripts/                       # Folder containing analysis scripts
│   ├── 01_fit_vanGenuchten.R      # Fit van Genuchten model and generate PSD curve
│   ├── 02_construct_PSD_matrix.R  # Construct pore size range matrix and compute θ differences
│   ├── 03_Pearson_correlation.R   # Calculate correlations between soil properties and θ differences
│   └── 04_random_forest_results.R # Compute Random Forest results and variable importance
└── README.md                      # This file
``` 

---

## Scripts Overview

### 01_fit_vanGenuchten.R
**Purpose:** Fit the van Genuchten model to soil water retention data and generate corresponding pore size distribution curves.  
**Input:** `data/Data_wrc.csv`  
**Output:** `data/Export_data.csv` (θ–h and θ–r data with R² and RMSE for each sample)

### 02_construct_PSD_matrix.R
**Purpose:** Construct the pore size range matrix and compute θ differences (`theta_diff`) for each pore range and sample.  
**Input:** `data/Export_data.csv`  
**Output:** `theta_diff_total_df` (θ differences with soil texture properties)

### 03_Pearson_correlation.R
**Purpose:** Calculate Pearson correlations between soil properties (OC, clay, silt, sand) and θ differences across pore size ranges.  
**Input:** `theta_diff_total_df`, `data/soil_properties.csv`  
**Output:** `pearson_corr.csv` (correlation coefficients, p-values, and slopes for each pore range)

### 04_random_forest_results.R
**Purpose:** Perform Random Forest analysis to evaluate variable importance and contributions for pore volume of each range.  
**Input:** `theta_diff_total_df`  
**Output:**  
- `RF_single_results_df.csv` : Random Forest model performance metrics (RMSE and R²) for each pore range (`r_min`, `r_max`)  
- `RF_variable_importance_df.csv` : Variable importance and contribution values for each predictor across all pore ranges

---

## Requirements

- **R (≥ 4.3)**
- Packages:  
```r
library(soilhypfit)
library(dplyr)
library(caret)
library(party)
```
---

## License
This repository is licensed under the MIT License.
You are free to use, modify, and distribute the code with attribution.

## Author
Haotian Wu  
haotian.wu@tum.de

Technical University of Munich  
TUM School of Life Sciences  
Chair of Soil Science  






