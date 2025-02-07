# Thesis_Armenia
## Vegetation Analysis
This project performs a detailed analysis of vegetation using various diversity indices, diversity profiles, and generalized additive models (GAM). The vegetation and environmental data are read from an Excel file, and several diversity indices and beta diversity components are calculated.

## Requirements
R packages:
- `readxl`
- `vegan`
- `ggplot2`
- `stats`
- `dplyr`
- `patchwork`
- `adespatial`
- `betapart`
- `mgcv`
- `tidyr`

## Analysis
The script file analisi_vegetazione.R performs the following analyses:

Calculation of Diversity Indices:
- Shannon Index
- Simpson Index
- Species Richness
- Total Abundance
- Margalef Index
- Menhinick Index
- Average Evenness Diversity (AED)

Renyi Diversity Profile:
- Calculation of the Renyi diversity profile
- Visualization of the diversity profile

Non-metric Multidimensional Scaling (NMDS):
- Calculation of the Bray-Curtis dissimilarity metric
- Visualization of NMDS

Permutational Analysis of Variance (PERMANOVA):
- PERMANOVA with Bray-Curtis index
- PERMANOVA with Jaccard index

Beta Diversity with the betapart package:
- Calculation of beta diversity
- Extraction of turnover and nestedness components

Generalized Additive Models (GAM):
- Study of the influence of elevation on species richness and Shannon diversity
- Visualization of GAM results
