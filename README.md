# Script_Armenia - Analisi di Alfa, Beta, Gamma Diversità e Lepidotteri

## Descrizione
Questo progetto effettua un'analisi approfondita sulla diversità di specie vegetali e Lepidotteri in Armenia. Utilizza vari indici di diversità, decomposizione della diversità beta e modelli additivi generalizzati (GAM) per studiare i pattern di biodiversità.

## Requisiti

### Pacchetti R richiesti:
required_packages:
-  "readxl",       # Lettura di file Excel
-  "vegan",        # Calcolo di indici di diversità e beta diversità
-  "ggplot2",      # Creazione di grafici e visualizzazioni
-  "stats",        # Funzioni statistiche di base
-  "dplyr",        # Manipolazione dei dati
-  "patchwork",    # Combinazione di grafici
-  "adespatial",   # Analisi spaziali e profili di diversità di Rényi
-  "betapart",     # Calcolo e partizionamento della diversità beta
-  "mgcv",         # Costruzione di modelli additivi generalizzati (GAM)
-  "tidyr",        # Trasformazione dei dati in formato lungo
-  "ggdendro",     # Creazione di dendrogrammi
-  "indicspecies"  # Analisi di specie indicatrici

## Installazione automatica dei pacchetti mancanti
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)
if(length(missing_packages)) install.packages(missing_packages)

### Caricamento librerie
lapply(required_packages, library, character.only = TRUE)
Workflow di Analisi
### 1. Caricamento dati e pre-elaborazione
r
Copia
Modifica
### Definizione dei percorsi dei file
file_path <- "C:/Users/Martina/Desktop/Tesi/Final_merged_file_Italian_and_Armenian_data.xlsx"
file_path_alien <- "C:/Users/Martina/Desktop/Tesi/Armenian_sp_alien.xlsx"
### 2. Analisi della diversità alfa
Calcolo degli indici di diversità:
- Indice di Shannon
- Indice di Simpson inverso
- Ricchezza specifica
- Abbondanza totale
- Indici di Margalef e Menhinick
- Diversità Effettiva Assoluta (AED)
Visualizzazione:
- Grafici delle metriche in funzione dell'altitudine e della distanza dalle strade
### 3. Profilo di diversità di Rényi
Calcolo e rappresentazione grafica del profilo a diverse scale di equità
### 4. Analisi NMDS (Non-metric Multidimensional Scaling)
Calcolo della dissimilarità di Bray-Curtis
Rappresentazione grafica della struttura della comunità
### 5. Analisi PERMANOVA
Test per verificare l'effetto dell'altitudine sulla composizione vegetale
Test Jaccard per le differenze nella composizione delle specie
### 6. Analisi della diversità beta
Calcolo di turnover e nestedness usando metriche di Sørensen e Jaccard
Analisi NMDS per visualizzare le differenze tra siti
### 7. Modelli additivi generalizzati (GAMs)
Modelli GAM per analizzare l'effetto di altitudine e distanza su:
Diversità alfa
Diversità gamma
Composizione dei Lepidotteri
### 8. Analisi Lepidotteri e relazione con la vegetazione
Calcolo della ricchezza delle specie di Lepidotteri
Analisi della correlazione tra diversità vegetale e Lepidotteri
