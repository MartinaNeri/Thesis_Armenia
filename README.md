# Thesis_Armenia
  ## Vegetation Analysis
  Questo progetto esegue un'analisi dettagliata della vegetazione utilizzando diversi indici di diversità, profili di diversità e modelli additivi generalizzati (GAM). I dati sulla vegetazione e sull'ambiente vengono letti da un file Excel e vengono calcolati vari indici di diversità e componenti della diversità beta.

  ## Requirements
    R packages necessari per l'analisi:
    -readxl – Per leggere i dati da file Excel
    -vegan – Per calcolare gli indici di diversità e la diversità beta
    -ggplot2 – Per creare grafici e visualizzazioni
    -stats – Per eseguire analisi statistiche
    -dplyr – Per manipolare i dati
    -patchwork – Per unire i grafici
    -adespatial – Per calcolare il profilo di diversità di Renyi
    -betapart – Per calcolare la diversità beta e partizionare turnover/nestedness
    -mgcv – Per costruire modelli additivi generalizzati (GAM)
    -tidyr – Per trasformare i dati in formato long
    -Analysis
    
  Il file di script analisi_vegetazione.R esegue le seguenti analisi:
  ## 1. Calcolo degli Indici di Diversità:
    -Indice di Shannon
    -Indice di Simpson (Inverse Simpson)
    -Ricchezza specifica
    -Abbondanza totale
    -Indice di Margalef
    -Indice di Menhinick
    -Average Evenness Diversity (AED)
    -Viene inoltre creato un boxplot per visualizzare la distribuzione degli indici calcolati.
	
 ## 2. Profilo di Diversità di Renyi:
 	Calcolo del profilo di diversità di Renyi a diverse scale di equità
	Visualizzazione del profilo di diversità con grafico migliorato
 
 ## 3. Non-metric Multidimensional Scaling (NMDS):
 	Calcolo della metrica di dissimilarità di Bray-Curtis	
	Applicazione della tecnica NMDS per rappresentare la struttura della comunità vegetale
 	Visualizzazione dei risultati in un grafico bidimensionale con etichette numeriche per i campioni
	
 ## 4. Analisi della Varianza Permutazionale (PERMANOVA):
 	PERMANOVA utilizzando l'indice di Bray-Curtis per testare l'influenza dell'altitudine sulla composizione della vegetazione
	PERMANOVA utilizzando l'indice di Jaccard per valutare la significatività delle differenze nella composizione delle specie
 
 ## 5. Analisi della Diversità Beta con il pacchetto betapart:
 	Calcolo della diversità beta utilizzando le metriche di Sorensen e Jaccard
	Partizionamento della diversità beta nelle componenti di turnover e nestedness
 	Calcolo della diversità beta multi-sito per valutare il grado di differenziazione tra le comunità vegetali
	
 ## 6. Modelli Additivi Generalizzati (GAM):
 	Studio dell'influenza dell'altitudine sulla ricchezza delle specie e sulla diversità di Shannon
	Visualizzazione dei risultati dei modelli GAM per interpretare meglio le relazioni tra variabili ambientali e biodiversità
 
 ## File di Input
 	I dati vengono caricati da un file Excel (Final_merged_file_Italian_and_Armenian_data.xlsx), che contiene:
	-Foglio "spec_veg": Dati sulla vegetazione (composizione delle specie)
 	-Foglio "env_veg_partial": Dati ambientali (altitudine, esposizione, distanza dalle strade, ecc.)
