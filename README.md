# Script_Armenia - Analisi di Diversità

## Descrizione
Analisi riassuntiva della biodiversità (alfa, beta, gamma) per vegetazione e lepidotteri in Armenia, 
utilizzando modelli GAM, NMDS, PERMANOVA e profili di Rényi.

## Requisiti
required_packages <- c("readxl", "vegan", "ggplot2", "stats", "dplyr", "patchwork",
                       "adespatial", "betapart", "mgcv", "tidyr", "ggdendro", "indicspecies", "ggrepel")

## Workflow
1. **Dati:**  
   - Caricamento dei file Excel (vegetazione, lepidotteri, specie aliene/autoctone).  
   - Pre-elaborazione: creazione di matrici presenza/assenza, trasposizione e calcolo delle ricchezze.

2. **Diversità Alfa:**  
   - Calcolo di indici (Shannon, Simpson, Richness, AED).  
   - Modellizzazione con GLM e GAM e creazione di grafici in funzione di altitudine, distanza dalla strada e dagli insediamenti.

3. **Diversità Beta:**  
   - Calcolo di dissimilarità (Bray-Curtis, Jaccard) e partizionamento (turnover, nestedness).  
   - NMDS, dendrogrammi e analisi pairwise con variabili ambientali.

4. **Diversità Gamma:**  
   - Stima della diversità totale (osservata e tramite metodi come Chao1).  
   - Profilo di Diversità di Rényi:**  
     - Calcolo del profilo a scale diverse (α).  
     - Visualizzazione grafica e raggruppamento per altitudine.

5. **Lepidotteri:**  
   - Calcolo della ricchezza delle specie di lepidotteri.  
   - Analisi GAM e correlazione (Spearman) tra diversità vegetale (AED) e ricchezza dei lepidotteri.
