# Caricamento dei pacchetti
library(readxl)    # Per leggere i dati da file Excel
library(vegan)     # Per calcolare gli indici di diversità
library(ggplot2)   # Per creare grafici e visualizzazioni
library(stats)     # Per eseguire analisi statistiche
library(dplyr)     # Per manipolare i dati
library(patchwork) # Per unire i grafici
library(adespatial) # Per calcolare il profilo di diversità di Renyi
library(betapart)  # Per calcolare la diversità beta
library(mgcv)      # Per i modelli additivi generalizzati
library(tidyr)     # Per trasformare i dati in formato long

# Caricamento dei dati da file Excel
file_path <- "C:/Users/Martina/Desktop/Tirocinio/Tirocinio Armenia/Final_merged_file_Italian_and_Armenian_data.xlsx"

############## 1. VEGETATION ANALYSIS 

# Lettura dei dati di vegetazione e ambiente
veg_data <- read_excel(file_path, sheet = "spec_veg")
env_veg <- read_excel(file_path, sheet = "env_veg_partial")
veg_altitude <- as.numeric(as.matrix(env_veg[1, 2:84]))  # Converte in vettore numerico
veg_aspect <- as.numeric(as.matrix(env_veg[2, 2:84]))
veg_street_d <- as.numeric(as.matrix(env_veg[3, 2:84]))

# Mostra la struttura dei dati caricati
str(veg_data)
str(env_veg)

# Trasforma il dataset in una matrice (escludendo la prima colonna)
veg_matrix <- as.matrix(veg_data[,-1])

# Trasponi la matrice per calcolare gli indici per ogni colonna
veg_matrix_t <- t(veg_matrix)

# Converte la matrice in una matrice di presenza/assenza
veg_matrix_pa <- ifelse(veg_matrix_t > 0, 1, 0)

############### 1.1 INDICES AND BOXPLOT ###############

# Calcolo dell'indice di Shannon per ogni colonna
shannon_diversity <- diversity(veg_matrix_t, index = "shannon")
print(shannon_diversity)

# Calcolo dell'indice di Simpson per ogni colonna
simpson_diversity <- diversity(veg_matrix_t, index = "invsimpson")
print(simpson_diversity)

# Calcolo della ricchezza delle specie e dell'abbondanza totale per ogni colonna
veg_richness <- colSums(veg_matrix > 0)
print(veg_richness)
total_abundance <- colSums(veg_matrix)
print(total_abundance)

# Calcolo dell'indice di Margalef per ogni colonna
margalef <- (veg_richness - 1) / log(total_abundance)
print(margalef)

# Calcolo dell'indice di Menhinick per ogni colonna
menhinick <- veg_richness / sqrt(total_abundance)
print(menhinick)

# Calcolo dell'AED per ogni colonna
H0 <- veg_richness
H1 <- exp(shannon_diversity)
H2 <- simpson_diversity
aed <- H0 + ((H1^2)/(2*H2))
print(aed)

# Crea una matrice dagli indici calcolati
index_matrix <- rbind(veg_altitude, veg_richness, simpson_diversity, shannon_diversity, aed, margalef, menhinick)
print(index_matrix)

# Trasforma la matrice in un data frame per la visualizzazione
index_df <- as.data.frame(t(index_matrix))
colnames(index_df) <- c("Altitude", "Richness", "Simpson", "Shannon", "AED", "Margalef", "Menhinick")

# Controllo dei dati
str(veg_altitude)
str(shannon_diversity)
str(simpson_diversity)

# Creazione di un data frame con altitudine, Shannon, Simpson e altri indici
env_index <- data.frame(
  Altitude = veg_altitude,
  Distance = veg_street_d,
  Aspect = veg_aspect,
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  AED = aed
)
print(env_index)

############### 2. RENYI DIVERSITY PROFILE ###############

# Valutare la diversità considerando diverse scale di equità
# Il profilo mostra come la diversità cambia in base all’enfasi sulle specie rare o dominanti.
# Calcolo del profilo di diversità Renyi
renyi_profile <- renyi(veg_matrix_t, scales = seq(0, 3, by = 0.5))
print(renyi_profile)

#𝛼= 0 → Ricchezza specifica (S): Conta il numero di specie senza considerare l'abbondanza relativa.
#𝛼= 1 → Indice di Shannon (𝐻′): Tiene conto della ricchezza e dell'equità.
#𝛼= 2 → Indice di Simpson (1/𝜆): Pone maggiore enfasi sulle specie dominanti.
#𝛼→ ∞ → L'indice è dominato dalla specie più abbondante.
# In generale quando α aumenta, l'indice diventa sempre più influenzato dalle specie più abbondanti nella comunità.
# Valori elevati di diversità indicano una distribuzione più equa delle abbondanze tra le specie.

# In molte righe, i valori tendono a diminuire all'aumentare di α, indicando che
# le comunità sono dominate da poche specie abbondanti.
# Ad esempio, la prima riga (1.0) passa da 1.0986 (𝛼=0) a 0.7772 (α=3), segnalando che
#alcune specie dominano la comunità.

#Alcune comunità hanno valori alti per tutti i valori di α (es. righe 4.0, 5.0, 7.1, 8.5), 
#suggerendo una maggiore equità tra le specie presenti.
#Altre comunità mostrano un rapido calo (es. righe 5.5, 6.1, 10.2), suggerendo che poche specie dominano.

#Alcune righe (es. 13.0, 14.1, 14.2) mostrano valori costanti per tutti i valori di α, 
#indicando che la comunità ha un'uniformità perfetta o pochissime specie.

# Trasforma il profilo di diversità Renyi in un formato long
renyi_df <- pivot_longer(as.data.frame(renyi_profile), cols = everything(), names_to = "scale", values_to = "diversity")
renyi_df$scale <- as.numeric(renyi_df$scale)
print(renyi_df)

# Plot del profilo di diversità Renyi migliorato
renyi_plot <- ggplot(renyi_df, aes(x = scale, y = diversity)) +
  geom_line(color = "blue", linewidth = 1) +  # Linea blu più spessa
  geom_point(color = "red", size = 2) +  # Punti rossi
  labs(title = "Renyi Diversity Profile", x = "Scale parameter", y = "Diversity") +
  theme_minimal() +  # Tema minimalista per migliorare la leggibilità
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Centrare il titolo e renderlo più grande
    axis.title = element_text(size = 12),  # Dimensione dei titoli degli assi
    axis.text = element_text(size = 10)  # Dimensione dei testi degli assi
  ) +
  geom_vline(xintercept = seq(0, 3, by = 0.5), linetype = "dashed", color = "grey", size = 0.5) +  # Linee verticali tratteggiate
  geom_hline(yintercept = seq(0, max(renyi_df$diversity), by = 0.5), linetype = "dashed", color = "grey", size = 0.5)  # Linee orizzontali tratteggiate
print(renyi_plot)

#La linea blu rappresenta il profilo di diversità di Renyi per le tue campionature di vegetazione. 
#Ogni punto lungo la linea mostra la diversità calcolata per un particolare 
#parametro di scala (indicato sull'asse delle X). La scala di Renyi è un continuum di indici
#di diversità che va da 0 a 3 (in questo caso), e ogni valore sulla scala riflette 
#una diversa sensibilità alla rarità delle specie.
#Linee tratteggiate grigie: Le linee tratteggiate grigie sull'asse X e Y aiutano a visualizzare meglio i valori 
#specifici di diversità per ciascun parametro di scala e per ciascuna campionatura di vegetazione.
#I punti rossi rappresentano i valori specifici di diversità calcolati per ciascun parametro di scala. 
#Questi punti sono effettivamente i dati calcolati dalla funzione renyi per ogni campionatura di vegetazione.


############### 3. NON-METRIC MULTIDIMENSIONAL SCALING (NMDS) ###############

# Calcolo della metrica di dissimilarità di Bray-Curtis
veg_bray <- vegdist(veg_matrix_t, method = "bray")

#NMDS
#NMDS è una tecnica di ordinamento che cerca di rappresentare le relazioni di dissimilarità 
#tra le osservazioni in uno spazio a bassa dimensione (in questo caso, 2 dimensioni).
#k = 2 indica che stiamo cercando di rappresentare i dati in uno spazio bidimensionale.
#Stress: Lo "stress" è una misura di quanto bene l'ordinamento a bassa dimensione 
#rappresenta le relazioni di dissimilarità originali. Un valore di stress più basso indica una migliore rappresentazione.
#Procrustes: La trasformazione di Procrustes è utilizzata per confrontare diverse 
#soluzioni di ordinamento. rmse (Root Mean Square Error) e max resid (residuo massimo) 
#sono misure di quanto bene le soluzioni si allineano. Un valore di rmse più basso indica un migliore allineamento.

nmds <- metaMDS(veg_bray, k = 2)
print(nmds)

#global Multidimensional Scaling using monoMDS

#Data:     veg_bray 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.124161 
#Stress type 1, weak ties
#Best solution was not repeated after 20 tries
#The best solution was from try 15 (random start)
#Scaling: centring, PC rotation, halfchange scaling 
#Species: scores missing

#Il valore di stress è 0.1242231. Lo stress è una misura di quanto bene l'ordinamento 
#bidimensionale rappresenta le relazioni di dissimilarità originali tra i campioni. 
#Un valore di stress più basso indica una migliore rappresentazione (è accettabile)

# Ottenere le coordinate NMDS e aggiungere i numeri dei campioni
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)


# Plot NMDS con numeri dei campioni come etichette
nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, label = Sample)) +
  geom_point() +
  geom_text(vjust = -0.5, hjust = 0.5) +
  labs(title = "Non-metric Multidimensional Scaling (NMDS)", x = "NMDS1", y = "NMDS2") +
  theme_minimal()
print(nmds_plot)

#Ogni punto nel grafico rappresenta un campione di vegetazione.
#La posizione dei punti è determinata dalle relazioni di dissimilarità tra i campioni. 
#Campioni che sono più simili tra loro saranno rappresentati vicini nel grafico, 
#mentre campioni più dissimili saranno rappresentati più lontani.

############### 4. ANALISI DELLA VARIANZA PERMUTAZIONALE (PERMANOVA) ###############

# PERMANOVA con indice di Bray-Curtis
permanova_bray <- adonis2(veg_matrix_t ~ veg_altitude, method = "bray")
print(permanova_bray)

# PERMANOVA con indice di Jaccard
permanova_jaccard <- adonis2(veg_matrix_t ~ veg_altitude, method = "jaccard")
print(permanova_jaccard)

#Entrambe le analisi PERMANOVA con gli indici di Bray-Curtis e Jaccard indicano 
#che l'altitudine ha un effetto significativo sulla composizione delle specie di 
#vegetazione nel tuo dataset. I valori p (0.001) sono molto bassi, suggerendo una 
#forte evidenza contro l'ipotesi nulla (che l'altitudine non ha effetto sulla composizione delle specie).
#Bray-Curtis: Spiega circa il 7.887% della variazione nella composizione delle specie.
#Jaccard: Spiega circa il 5.594% della variazione nella composizione delle specie.
#Questi risultati possono essere utilizzati per supportare l'ipotesi che l'altitudine 
#influenzi la composizione delle comunità vegetali.

#I bassi valori di ( R^2 ) indicano che l'altitudine spiega solo una porzione limitata 
#della variazione totale nella composizione delle specie.
#Tuttavia, il fatto che il valore ( p ) sia molto basso (( p < 0.001 )) significa che l'effetto dell'altitudine è comunque statisticamente significativo. Questo suggerisce che, anche se l'effetto è piccolo, è reale e non dovuto al caso.

############### 5. DIVERSITÀ BETA CON BETAPART PACKAGE

# Calcolo della diversità beta utilizzando la matrice di presenza/assenza
beta_pair <- beta.pair(veg_matrix_pa, index.family = "sorensen")
print(beta_pair)

# Estrazione delle componenti turnover e nestedness
beta_turnover <- beta_pair$beta.sim
beta_nestedness <- beta_pair$beta.sne
print(beta_turnover)
print(beta_nestedness)

# Converti le matrici in data frame per la scrittura nei file CSV
#beta_turnover_df <- as.data.frame(as.matrix(beta_turnover))
#beta_nestedness_df <- as.data.frame(as.matrix(beta_nestedness))

# Scrittura dei dati nei file CSV
#write.csv(beta_turnover_df, "Beta_Turnover.csv", row.names = TRUE)
#write.csv(beta_nestedness_df, "Beta_Nestedness.csv", row.names = TRUE)


#########5.2 CALCOLO DELLA DIVERSITÀ BETA USANDO IL METODO DI BASELGA##########

# Calcolo della diversità beta utilizzando il metodo di Baselga

# Calcola la diversità beta totale
beta_total <- beta.pair(veg_matrix_pa, index.family = "jaccard")
print(beta_total)

# Partiziona la diversità beta nelle sue componenti di turnover e nestedness
beta_turnover_jaccard <- beta_total$beta.jtu
beta_nestedness_jaccard <- beta_total$beta.jne
print(beta_turnover_jaccard)
print(beta_nestedness_jaccard)

# Calcolo delle misure multiple-site
beta_multi <- beta.multi(veg_matrix_pa, index.family = "sorensen")
print(beta_multi)

#$beta.SIM
#[1] 0.98795
#rappresenta il turnover delle specie, cioè la proporzione di differenze nella composizione 
#delle specie dovuta alla sostituzione di specie tra i siti. Un valore vicino a 1 
#indica un alto turnover, suggerendo che la maggior parte della dissimilarità è dovuta 
#al fatto che le specie presenti in un sito non sono presenti in un altro.

#$beta.SNE
#[1] 0.006268038
#rappresenta il nestedness della diversità beta, cioè la proporzione di differenze 
#nella composizione delle specie dovuta alla perdita o guadagno di specie tra i siti. 
#Un valore vicino a 0 indica che il nestedness contribuisce poco alla dissimilarità complessiva.

#$beta.SOR
#[1] 0.9942181
#rappresenta la dissimilarità totale di Sørensen, combinando sia il turnover che il nestedness. 
#Un valore vicino a 1 indica un'alta dissimilarità totale tra i siti.


# Misure multiple-site per la famiglia di Jaccard
beta_multi_jaccard <- beta.multi(veg_matrix_pa, index.family = "jaccard")
print(beta_multi_jaccard)

#$beta.JTU
#[1] 0.9939385
#rappresenta il turnover delle specie nella famiglia di Jaccard. Analogamente a 
#βSIM, un valore vicino a 1 indica un alto turnover.

#$beta.JNE
#[1] 0.003162158
#appresenta il nestedness della diversità beta nella famiglia di Jaccard. 
#Analogamente a βSNE, un valore vicino a 0 indica che il nestedness contribuisce 
#poco alla dissimilarità complessiva.

#$beta.JAC
#[1] 0.9971007
#rappresenta la dissimilarità totale di Jaccard, combinando sia il turnover che 
#il nestedness. Un valore vicino a 1 indica un'alta dissimilarità totale tra i siti.


############### 5. DIVERSITÀ BETA CON BETAPART PACKAGE

# Calcolo della diversità beta utilizzando la matrice di presenza/assenza
beta_pair <- beta.pair(veg_matrix_pa, index.family = "jaccard")
print(beta_pair)

# Estrazione delle componenti turnover e nestedness
beta_turnover <- beta_pair$beta.jtu
beta_nestedness <- beta_pair$beta.jne
print(beta_turnover)
print(beta_nestedness)

# Converti le matrici in data frame per la scrittura nei file CSV
#beta_turnover_df <- as.data.frame(as.matrix(beta_turnover))
#beta_nestedness_df <- as.data.frame(as.matrix(beta_nestedness))

# Scrittura dei dati nei file CSV
#write.csv(beta_turnover_df, "Beta_Turnover.csv", row.names = TRUE)
#write.csv(beta_nestedness_df, "Beta_Nestedness.csv", row.names = TRUE)

# Aggiunta delle analisi proposte da Baselga

# Calcola la diversità beta totale
beta_total <- beta.pair(veg_matrix_pa, index.family = "jaccard")
print(beta_total)

# Partiziona la diversità beta nelle sue componenti di turnover e nestedness
beta_turnover_jaccard <- beta_total$beta.jtu
beta_nestedness_jaccard <- beta_total$beta.jne
print(beta_turnover_jaccard)
print(beta_nestedness_jaccard)

# Calcolo delle misure multiple-site
beta_multi_jaccard <- beta.multi(veg_matrix_pa, index.family = "jaccard")
print(beta_multi_jaccard)

#$beta.JTU
#[1] 0.9849822
#rappresenta il turnover delle specie nella famiglia di Jaccard. Un valore vicino a 1 
#indica un alto turnover, suggerendo che la maggior parte della dissimilarità tra 
#i siti è dovuta alla sostituzione delle specie. In altre parole, i siti tendono ad avere 
#specie differenti tra loro.

#$beta.JNE
#[1] 0.005062547
#rappresenta il nestedness della diversità beta nella famiglia di Jaccard. 
#Un valore vicino a 0 indica che il nestedness contribuisce poco alla dissimilarità complessiva, 
#suggerendo che le specie presenti in un sito non sono semplicemente un sottoinsieme 
#di quelle presenti in un altro sito. Il nestedness descrive quanto i siti con meno 
#specie sono sottoinsiemi di siti più ricchi di specie.


#$beta.JAC
#[1] 0.9900448
#rappresenta la dissimilarità totale di Jaccard, combinando sia il turnover che il nestedness.
#Un valore vicino a 1 indica un'alta dissimilarità totale tra i siti. Questo significa che 
#i siti tendono ad avere composizioni di specie molto diverse tra loro.


############### 6. GENERALIZED ADDITIVE MODEL (GAM) ###############

# Utilizziamo i modelli additivi generalizzati (GAM) per studiare l'effetto dell'altitudine 
# sulla ricchezza delle specie e sulla diversità di Shannon.

# GAM per studiare l'influenza dell'elevazione sulla ricchezza delle specie
# ho aumentato k con maggiore complessità
gam_richness <- gam(veg_richness ~ s(veg_altitude, k = 30), data = env_index)
summary(gam_richness)

#R-sq.(adj): 0.798 il coefficiente di determinazione aggiustato, che rappresenta la proporzione 
#della variabilità nella ricchezza delle specie spiegata dal modello. 
#Un valore di 0.798 indica che il modello spiega circa il 60.9% della variabilità totale.

# GAM per studiare l'influenza dell'elevazione sulla diversità di Shannon
gam_shannon <- gam(Shannon ~ s(Altitude), data = env_index)
summary(gam_shannon)

#R-sq.(adj): 0.651 il coefficiente di determinazione aggiustato, che rappresenta la proporzione
#della variabilità nella diversità di Shannon spiegata dal modello. Un valore di 0.651 indica 
#che il modello spiega circa il 65.1% della variabilità totale.

# Diagnostica del modello
par(mfrow = c(2, 2))
plot(gam_richness, residuals = TRUE, pch = 19, cex = 0.5)
gam.check(gam_richness)

#k-index: 1.08
#p-value: 0.71
#L'indice k è superiore a 1 e il p-value è alto (0.71), indicando che 
#la dimensione della base k scelta è adeguata e non c'è necessità di aumentarla ulteriormente.

# Previsione del modello GAM
predict_richness <- data.frame(Altitude = veg_altitude, Predicted = predict(gam_richness, type = "response"))
print(predict_richness)


#Il modello GAM suggerisce che l'altitudine ha un effetto non lineare sulla ricchezza delle specie. 
#La ricchezza delle specie sembra aumentare fino a un certo punto con l'altitudine, 
#raggiungendo un picco a altitudini intermedie, e poi diminuire a altitudini molto elevate. 
#Questi risultati possono essere utilizzati per comprendere meglio come la variazione 
#altitudinale influenzi la


