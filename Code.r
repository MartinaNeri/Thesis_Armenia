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

############### 1.1 INDICES AND BOXPLOT

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

############### 2. RENYI DIVERSITY PROFILE

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

# Plot del profilo di diversità Renyi
renyi_plot <- ggplot(renyi_df, aes(x = scale, y = diversity)) +
  geom_line() +
  labs(title = "Renyi Diversity Profile", x = "Scale parameter", y = "Diversity")
print(renyi_plot)


ggplot(renyi_df, aes(x = scale, y = diversity, group = interaction(row_number()), color = factor(row_number()))) +
  geom_line() +
  labs(title = "Renyi Diversity Profile", x = "Scale parameter", y = "Diversity") +
  theme_minimal()

############### 3. NON-METRIC MULTIDIMENSIONAL SCALING (NMDS)

# Calcolo della metrica di dissimilarità di Bray-Curtis
veg_bray <- vegdist(veg_matrix_t, method = "bray")

# NMDS
nmds <- metaMDS(veg_bray, k = 2)
print(nmds)

# Plot NMDS
nmds_plot <- ggplot(as.data.frame(scores(nmds)), aes(x = NMDS1, y = NMDS2)) +
  geom_point() +
  labs(title = "Non-metric Multidimensional Scaling (NMDS)", x = "NMDS1", y = "NMDS2")
print(nmds_plot)

############### 4. ANALISI DELLA VARIANZA PERMUTAZIONALE (PERMANOVA)

# PERMANOVA con indice di Bray-Curtis
permanova_bray <- adonis2(veg_matrix_t ~ veg_altitude, method = "bray")
print(permanova_bray)

# PERMANOVA con indice di Jaccard
permanova_jaccard <- adonis2(veg_matrix_t ~ veg_altitude, method = "jaccard")
print(permanova_jaccard)

############### 5. DIVERSITÀ BETA CON BETAPART PACKAGE

# Calcolo della diversità beta utilizzando la matrice di presenza/assenza
beta_pair <- beta.pair(veg_matrix_pa, index.family = "sorensen")
print(beta_pair)

# Estrazione delle componenti turnover e nestedness
beta_turnover <- beta_pair$beta.sim
beta_nestedness <- beta_pair$beta.sne
print(beta_turnover)
print(beta_nestedness)

############### 6. GENERALIZED ADDITIVE MODEL

# GAM per studiare l'influenza dell'elevazione sulla ricchezza delle specie
gam_richness <- gam(veg_richness ~ s(veg_altitude), data = env_index)
summary(gam_richness)

# GAM per studiare l'influenza dell'elevazione sulla diversità di Shannon
gam_shannon <- gam(Shannon ~ s(Altitude), data = env_index)
summary(gam_shannon)

# Previsione del modello GAM per la ricchezza delle specie
predict_richness <- data.frame(Altitude = veg_altitude, Predicted = predict(gam_richness, type = "response"))
print(predict_richness)

# Plot GAM per ricchezza delle specie
gam_richness_plot <- ggplot(predict_richness, aes(x = Altitude, y = Predicted)) +
  geom_line() +
  labs(title = "GAM - Species Richness vs Altitude", x = "Altitude", y = "Species Richness")
print(gam_richness_plot)

# Previsione del modello GAM per la diversità di Shannon
predict_shannon <- data.frame(Altitude = veg_altitude, Predicted = predict(gam_shannon, type = "response"))
print(predict_shannon)

# Plot GAM per diversità di Shannon
gam_shannon_plot <- ggplot(predict_shannon, aes(x = Altitude, y = Predicted)) +
  geom_line() +
  labs(title = "GAM - Shannon Diversity vs Altitude", x = "Altitude", y = "Shannon Diversity")
print(gam_shannon_plot)

