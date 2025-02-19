# Thesis_Armenia

## Vegetation and Lepidoptera Analysis

 This project conducts a detailed analysis of vegetation and Lepidoptera species 
 using various diversity indices, beta diversity partitioning, and Generalized Additive Models (GAMs). 
 The environmental and species data are extracted from an Excel file, 
 and multiple statistical analyses are performed to understand biodiversity patterns.

## Requirements

# R packages necessary for the analysis:
required_packages <- c(
  "readxl",    # For reading data from Excel files
  "vegan",     # For calculating diversity indices and beta diversity
  "ggplot2",   # For creating graphs and visualizations
  "stats",     # For performing statistical analyses
  "dplyr",     # For data manipulation
  "patchwork", # For combining plots
  "adespatial",# For calculating Renyi's diversity profile
  "betapart",  # For computing beta diversity and partitioning turnover/nestedness
  "mgcv",      # For building Generalized Additive Models (GAMs)
  "tidyr",     # For transforming data into long format
  "ggdendro"   # For creating dendrograms with ggplot2
)

# Install missing packages
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)
if(length(missing_packages)) install.packages(missing_packages)

# Load libraries
lapply(required_packages, library, character.only = TRUE)

## Analysis Workflow

# 1. Alpha Diversity Analysis
 - Calculation of Diversity Indices:
   - Shannon Index
   - Inverse Simpson Index
   - Species Richness
   - Total Abundance
   - Margalef and Menhinick indices
   - Absolute Effective Diversity (AED)
 - Visualization: Boxplots of calculated indices

# 2. Renyi's Diversity Profile
 - Computation of Renyi's diversity profile at different equitability scales
 - Enhanced graphical representation

# 3. Non-metric Multidimensional Scaling (NMDS)
 - Calculation of Bray-Curtis dissimilarity metric
 - NMDS application for visualizing vegetation community structure
 - Graphical representation with labeled sample points

# 4. Permutational Multivariate Analysis of Variance (PERMANOVA)
 - Bray-Curtis PERMANOVA to test the influence of altitude on vegetation composition
 - Jaccard PERMANOVA to assess significance in species composition differences

# 5. Beta Diversity Analysis (betapart package)
 - Calculation of beta diversity using Sorensen and Jaccard metrics
 - Partitioning beta diversity into turnover and nestedness components
 - Multi-site beta diversity estimation

# 6. Generalized Additive Models (GAMs)
 - Analysis of altitude’s influence on species richness and Shannon diversity
 - GAM visualization for interpreting biodiversity-environment relationships

# 7. Vegetation and Lepidoptera Analysis
 - Vegetation Gamma Diversity:
   - Observed and estimated species richness (Chao1 estimator)
   - GAM modeling of richness and Shannon index
   - Renyi diversity profile across altitudinal groups
 - Lepidoptera Analysis:
   - Calculation of Lepidoptera species richness
   - Regression models for altitude and road distance influence
   - NMDS visualization of Lepidoptera community structure
 - Vegetation-Lepidoptera Relationship:
   - Correlation analysis between vegetation diversity and Lepidoptera richness
   - GAM models assessing biodiversity interaction

## Input Files

# Data is loaded from an Excel file
file_path <- "Final_merged_file_Italian_and_Armenian_data.xlsx"

# Sheets in the Excel file:
sheets <- list(
  spec_veg = "Vegetation composition data",
  env_veg_partial = "Environmental data (altitude, road distance, etc.)",
  spec_lep = "Lepidoptera composition data",
  env_lep = "Environmental data for Lepidoptera analysis"
)

# Read data
veg_data <- read_excel(file_path, sheet = "spec_veg")
env_veg  <- read_excel(file_path, sheet = "env_veg_partial")
lep_data <- read_excel(file_path, sheet = "spec_lep")
env_lep  <- read_excel(file_path, sheet = "env_lep")

 This script provides a comprehensive framework for biodiversity analysis in Armenia,
 integrating ecological statistics, modeling, and visualization tools.
