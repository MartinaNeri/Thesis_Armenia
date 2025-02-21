### Script_Armenia - Analisi di Alfa, Beta, Gamma Diversità e Lepidotteri ####
#### 1. Caricamento librerie e file path #### 
library(readxl)       # Per leggere file Excel
library(vegan)        # Per il calcolo degli indici di diversità
library(ggplot2)      # Per la creazione dei grafici
library(stats)        # Per funzioni statistiche di base
library(dplyr)        # Per la manipolazione dei dati
library(patchwork)    # Per combinare grafici
library(adespatial)   # Per analisi spaziali e profili di diversità
library(betapart)     # Per l’analisi della diversità beta
library(mgcv)         # Per i modelli additivi generalizzati (GAM)
library(tidyr)        # Per la trasformazione dei dati (formato long)
library(ggdendro)     # Per la creazione di dendrogrammi
library(indicspecies) # (Opzionale) per analisi di specie indicatrici

# Definizione dei percorsi dei file
file_path <- "C:/Users/Martina/Desktop/Tesi/Final_merged_file_Italian_and_Armenian_data.xlsx"
file_path_alien <- "C:/Users/Martina/Desktop/Tesi/Armenian_sp_alien.xlsx"

#### 2. Caricamento e Pre-elaborazione dei Dati ####
###### 2.1 Dati dei Subplot (Vegetazione) ######
veg_data_subplot <- read_excel(file_path, sheet = "spec_veg")
env_veg_subplot  <- read_excel(file_path, sheet = "env_veg_partial")
veg_altitude_subplot <- as.numeric(as.matrix(env_veg_subplot[1, 2:84]))
veg_street_d_subplot <- as.numeric(as.matrix(env_veg_subplot[3, 2:84]))
veg_matrix_subplot   <- as.matrix(veg_data_subplot[,-1])
veg_matrix_subplot_t <- t(veg_matrix_subplot)                # righe = siti, colonne = specie
veg_matrix_subplot_pa<- ifelse(veg_matrix_subplot_t > 0, 1, 0)  # matrice presenza/assenza

###### 2.2 Dati a Livello di Plot (Vegetazione) ######
veg_plot <- read_excel(file_path, sheet = "spec_veg_plot")
env_plot  <- read_excel(file_path, sheet = "env_veg_plot")
veg_altitude_plot <- as.numeric(as.matrix(env_plot[1, 2:15]))
veg_street_d_plot <- as.numeric(as.matrix(env_plot[3, 2:15]))
veg_matrix_plot   <- as.matrix(veg_plot[,-1])
veg_matrix_plot_t <- t(veg_matrix_plot)
veg_matrix_plot_pa<- ifelse(veg_matrix_plot_t > 0, 1, 0)

###### 2.3 Dati sui Lepidotteri ######
lep_data_plot <- read_excel(file_path, sheet = "spec_lep")
env_lep_plot  <- read_excel(file_path, sheet = "env_lep")
lep_altitude_plot <- as.numeric(as.matrix(env_lep_plot[1, 2:15]))
lep_dist_plot     <- as.numeric(as.matrix(env_lep_plot[2, 2:15]))
lep_matrix_plot   <- as.matrix(lep_data_plot[,-1])
lep_matrix_plot_t <- t(lep_matrix_plot)
lep_richness_plot <- colSums(lep_matrix_plot > 0)
lep_df <- as.data.frame(t(rbind(
  Altitude = lep_altitude_plot,
  Distance = lep_dist_plot,
  Richness = lep_richness_plot
)))

###### 2.4 Dati su Specie Aliene e Autoctone ######
sp_alien  <- read_excel(file_path_alien, sheet = "abb_aliene")
sp_autoct <- read_excel(file_path_alien, sheet = "abb_autoctone")
abb_aliene    <- as.matrix(sp_alien[,-1])
abb_autoctone <- as.matrix(sp_autoct[,-1])
abb_aliene_t    <- t(abb_aliene)
abb_autoctone_t <- t(abb_autoctone)
# Calcolo del rapporto specie aliene su totale
endem_ratio <- rowSums(abb_aliene_t) / (rowSums(abb_autoctone_t) + rowSums(abb_aliene_t))
endem_df <- data.frame(
  Altitude = veg_altitude_subplot,
  Distance = veg_street_d_subplot,
  EndemRatio = endem_ratio
)

#### 3. DIVERSITÀ ALFA ####

###### 3.1 Calcolo degli Indici di Diversità Alfa ######
# Calcolo della ricchezza di specie, indice di Shannon, Simpson e indice composito AED
shannon_diversity_subplot <- diversity(veg_matrix_subplot_t, index = "shannon")
simpson_diversity_subplot <- diversity(veg_matrix_subplot_t, index = "invsimpson")
veg_richness_subplot      <- colSums(veg_matrix_subplot > 0)
total_abundance_subplot   <- colSums(veg_matrix_subplot)
margalef_subplot          <- (veg_richness_subplot - 1) / log(total_abundance_subplot)
menhinick_subplot         <- veg_richness_subplot / sqrt(total_abundance_subplot)
H0_subplot <- veg_richness_subplot         # Numero di specie
H1_subplot <- exp(shannon_diversity_subplot)  # Diversità effettiva (trasformazione di Shannon)
H2_subplot <- simpson_diversity_subplot       # Simpson (inverso)
aed_subplot <- H0_subplot + ((H1_subplot^2)/(2*H2_subplot))

###### 3.2 Creazione del Dataframe per l'Analisi Alfa ######
env_index_subplot <- data.frame(
  Altitude = veg_altitude_subplot,
  Distance = veg_street_d_subplot,
  Shannon  = shannon_diversity_subplot,
  Simpson  = simpson_diversity_subplot,
  AED      = aed_subplot,
  Richness = veg_richness_subplot
)

###### 3.3 Modellizzazione della Relazione tra Diversità Alfa e Variabili Ambientali ######
######### 3.3.1 Utilizzo di GLM #########
glm_shannon <- glm(Shannon ~ Altitude + Distance, data = env_index_subplot, family = gaussian())
glm_simpson <- glm(Simpson ~ Altitude + Distance, data = env_index_subplot, family = gaussian())
glm_aed     <- glm(AED ~ Altitude + Distance, data = env_index_subplot, family = gaussian())

######### 3.3.2 Utilizzo di GAM (per relazioni non lineari) #########
gam_shannon <- gam(Shannon ~ s(Altitude) + s(Distance), data = env_index_subplot)
gam_simpson <- gam(Simpson ~ s(Altitude) + s(Distance), data = env_index_subplot)
gam_aed     <- gam(AED ~ s(Altitude) + s(Distance), data = env_index_subplot)

summary(glm_shannon)
summary(glm_simpson)
summary(glm_aed)

summary(gam_shannon)
# Altitudine: p < 2e-16 → Effetto altamente significativo e non lineare.
# Distanza dalle strade: p = 0.0952 → Effetto quasi significativo.
# Conclusione: L'altitudine ha un forte effetto non lineare sulla diversità di 
#Shannon. L'effetto della distanza è borderline, ma potrebbe avere un impatto in alcune zone.
summary(gam_simpson)
# Altitudine: p < 2e-16 → Effetto altamente significativo e non lineare.
# Distanza dalle strade: p = 0.0947 → Effetto quasi significativo.
# Conclusione: Come per Shannon, anche qui l'altitudine è un predittore significativo con un effetto non lineare. 
summary(gam_aed)
# Altitudine: p < 2e-16 → Effetto significativo e non lineare.
# Distanza dalle strade: p = 0.0201 → Effetto significativo.
# Conclusione: Sia l'altitudine che la distanza influenzano significativamente l'AED, suggerendo che ci sono pattern complessi nelle relazioni tra la diversità effettiva e le variabili ambientali.

# i modelli GAM mostrano un chiaro effetto non lineare su tutti gli indici di diversità. 
# Probabilmente ci sono fasce altitudinali specifiche in cui la diversità aumenta 
# o diminuisce bruscamente (es. diversità maggiore a medie altitudini).

###### 3.4 Visualizzazione Grafica per la Diversità Alfa ######
######### 3.4.1 Grafici in funzione dell'Altitudine #########
p1 <- ggplot(env_index_subplot, aes(x = Altitude, y = Shannon)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "Shannon Index vs Altitude", x = "Altitude (m)", y = "Shannon Index")
p2 <- ggplot(env_index_subplot, aes(x = Altitude, y = Simpson)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "Simpson Index vs Altitude", x = "Altitude (m)", y = "Simpson Index")
p3 <- ggplot(env_index_subplot, aes(x = Altitude, y = AED)) +
  geom_point(color = 'purple') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'pink') +
  labs(title = "AED vs Altitude", x = "Altitude (m)", y = "AED")
######### 3.4.2 Grafici in funzione della Distanza dalle Strade #########
p4 <- ggplot(env_index_subplot, aes(x = Distance, y = Shannon)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "Shannon Index vs Distance", x = "Distance (m)", y = "Shannon Index")
p5 <- ggplot(env_index_subplot, aes(x = Distance, y = Simpson)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "Simpson Index vs Distance", x = "Distance (m)", y = "Simpson Index")
p6 <- ggplot(env_index_subplot, aes(x = Distance, y = AED)) +
  geom_point(color = 'purple') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'pink') +
  labs(title = "AED vs Distance", x = "Distance (m)", y = "AED")
# Visualizzazione combinata
(p1 | p2 | p3) / (p4 | p5 | p6)


###### 3.5 aliene e autoctone ######
# Calcolo indici di diversità per specie aliene e autoctone
shannon_aliene <- diversity(abb_aliene_t, index = "shannon")
simpson_aliene <- diversity(abb_aliene_t, index = "invsimpson")

shannon_autoctone <- diversity(abb_autoctone_t, index = "shannon")
simpson_autoctone <- diversity(abb_autoctone_t, index = "invsimpson")

# Dataframe di confronto
diversity_df <- data.frame(
  Altitude = veg_altitude_subplot,
  Distance = veg_street_d_subplot,
  Shannon_Aliene = shannon_aliene,
  Simpson_Aliene = simpson_aliene,
  Shannon_Autoctone = shannon_autoctone,
  Simpson_Autoctone = simpson_autoctone,
)

# Controllo se ci sono ancora NA
any(is.na(diversity_df$Richness_Aliene))
any(is.na(diversity_df$Richness_Autoctone))

# Rimuovi eventuali NA
diversity_df <- na.omit(diversity_df)

######### 3.5.1 Analisi Statistica e grafici #########
# Modelli GAM
# Specie aliene
gam_shannon_aliene <- gam(Shannon_Aliene ~ s(Altitude) + s(Distance), data = diversity_df)
gam_simpson_aliene <- gam(Simpson_Aliene ~ s(Altitude) + s(Distance), data = diversity_df)

# Specie autoctone
gam_shannon_autoctone <- gam(Shannon_Autoctone ~ s(Altitude) + s(Distance), data = diversity_df)
gam_simpson_autoctone <- gam(Simpson_Autoctone ~ s(Altitude) + s(Distance), data = diversity_df)

print(summary(gam_shannon_aliene))
print(summary(gam_simpson_aliene))
print(summary(gam_shannon_autoctone))
print(summary(gam_simpson_autoctone))

# Confronto tra aliene e autoctone
shannon_ttest <- t.test(diversity_df$Shannon_Aliene, diversity_df$Shannon_Autoctone)
simpson_ttest <- t.test(diversity_df$Simpson_Aliene, diversity_df$Simpson_Autoctone)

# Visualizza i risultati
print(shannon_ttest)
print(simpson_ttest)


# Visualizzazione Grafica 
# Grafici per aliene
paa1 <- ggplot(diversity_df, aes(x = Altitude, y = Shannon_Aliene)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "Shannon (Aliene) vs Altitudine", x = "Altitudine (m)", y = "Indice di Shannon")

# Grafici per autoctone
paa2 <- ggplot(diversity_df, aes(x = Altitude, y = Shannon_Autoctone)) +
  geom_point(color = 'green') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkgreen') +
  labs(title = "Shannon (Autoctone) vs Altitudine", x = "Altitudine (m)", y = "Indice di Shannon")

# Visualizzazione combinata
(paa1 | paa2)

#### 4. DIVERSITÀ BETA ####

######  4.1 Calcolo delle Dissimilarità tra Siti ######
#########  4.1.1 Dissimilarità basata sulle abbondanze (Bray-Curtis) a livello di Plot #########
veg_bray_beta <- vegdist(veg_matrix_plot_t, method = "bray")
#########  4.1.2 Dissimilarità binarie (Sørensen e Jaccard) sui dati dei Subplot #########
dissim_sorensen <- vegdist(veg_matrix_subplot_pa, method = "bray", binary = TRUE)
dissim_jaccard   <- vegdist(veg_matrix_subplot_pa, method = "jaccard")

###### 4.2 Decomposizione della Diversità Beta ###### 
beta_pair       <- beta.pair(veg_matrix_subplot_pa, index.family = "sorensen")
beta_turnover   <- beta_pair$beta.sim      # Turnover (sostituzione specie)
beta_nestedness <- beta_pair$beta.sne      # Nestedness (perdita/guadagno specie)
beta_total      <- beta_pair$beta.sor      # Dissimilarità totale
beta_multi_sorensen <- beta.multi(veg_matrix_subplot_pa, index.family = "sorensen")
print(beta_multi_sorensen)

###### 4.3 Analisi PERMANOVA ###### 
permanova_bray_beta   <- adonis2(veg_matrix_plot_t ~ veg_altitude_plot, method = "bray")
print(permanova_bray_beta)
permanova_jaccard_beta<- adonis2(veg_matrix_plot_pa ~ veg_altitude_plot, method = "jaccard")
print(permanova_jaccard_beta)

## 4.4 NMDS per Visualizzare le Differenze nella Composizione delle Comunità
nmds_beta <- metaMDS(veg_bray_beta, k = 2, trymax = 100)
nmds_scores_beta <- as.data.frame(scores(nmds_beta))
nmds_scores_beta$Sample <- rownames(nmds_scores_beta)
nmds_plot_beta <- ggplot(nmds_scores_beta, aes(x = NMDS1, y = NMDS2, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5) +
  labs(title = "NMDS: Composizione delle Comunità", x = "NMDS1", y = "NMDS2") +
  theme_minimal()
print(nmds_plot_beta)

## 4.5 Analisi della Relazione tra Dissimilarità e Altitudine
bray_altitude_df <- data.frame(
  Altitude = veg_altitude_plot,
  Dissimilarity = rowMeans(as.matrix(veg_bray_beta))
)
plot_br_alt <- ggplot(bray_altitude_df, aes(x = Altitude, y = Dissimilarity)) +
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "gam", se = TRUE, color = "black") +
  labs(title = "Altitudine vs Dissimilarità (Bray-Curtis)",
       x = "Altitude", y = "Dissimilarità Media") +
  theme_minimal()
print(plot_br_alt)
summary(gam(Dissimilarity ~ Altitude, data = bray_altitude_df))

## 4.6 Creazione di un Dendrogramma Basato sulla Dissimilarità di Bray-Curtis
bray_clustering_beta <- hclust(veg_bray_beta, method = "average")
dendro_data_beta <- as.dendrogram(bray_clustering_beta)
dendro_plot_beta <- ggdendrogram(dendro_data_beta, theme_dendro = FALSE) +
  ggtitle("Dendrogramma: Bray-Curtis") +
  theme_minimal()
print(dendro_plot_beta)

## 4.7 Analisi Pairwise e Relazione con Variabili Ambientali
beta_bray_pairwise <- vegdist(veg_matrix_plot_t, method = "bray")
beta_bray_df <- as.data.frame(as.matrix(beta_bray_pairwise))
beta_bray_df$Site1 <- rownames(beta_bray_df)
beta_bray_long <- beta_bray_df %>%
  pivot_longer(cols = -Site1, names_to = "Site2", values_to = "Beta_Diversity") %>%
  filter(Site1 != Site2)

# Assegnazione dei nomi dei siti (presupponendo che la prima colonna di veg_plot contenga i nomi)
site_names <- veg_plot[[1]]
rownames(veg_matrix_plot) <- site_names
colnames(veg_matrix_plot_t) <- site_names

env_data <- data.frame(
  Site = colnames(veg_matrix_plot),
  Altitude = veg_altitude_plot,
  Street_Distance = veg_street_d_plot
)

beta_bray_long <- beta_bray_long %>%
  left_join(env_data, by = c("Site1" = "Site")) %>%
  rename(Altitude1 = Altitude, Street_Distance1 = Street_Distance) %>%
  left_join(env_data, by = c("Site2" = "Site")) %>%
  rename(Altitude2 = Altitude, Street_Distance2 = Street_Distance) %>%
  mutate(
    Altitude_Diff = abs(Altitude1 - Altitude2),
    Street_Distance_Diff = abs(Street_Distance1 - Street_Distance2)
  )
beta_bray_long_clean <- beta_bray_long %>%
  mutate(Pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-Pair)

plot_altitude_beta <- ggplot(beta_bray_long_clean, aes(x = Altitude_Diff, y = Beta_Diversity)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "gam", se = TRUE, color = "black") +
  labs(title = "Differenza di Altitudine vs Beta Diversità",
       x = "Differenza di Altitudine",
       y = "Beta Diversità (Bray-Curtis)") +
  theme_minimal()
print(plot_altitude_beta)
summary(gam(Beta_Diversity ~ Altitude_Diff, data = beta_bray_long_clean))

plot_street_beta <- ggplot(beta_bray_long_clean, aes(x = Street_Distance_Diff, y = Beta_Diversity)) +
  geom_point(alpha = 0.5, color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "Differenza di Distanza vs Beta Diversità",
       x = "Differenza di Distanza",
       y = "Beta Diversità (Bray-Curtis)") +
  theme_minimal()
print(plot_street_beta)
summary(gam(Beta_Diversity ~ Street_Distance_Diff, data = beta_bray_long_clean))

## 4.8 Analisi Specie Endemiche/Autoctone (Opzionale)
gam_shannon_aliene    <- gam(diversity(abb_aliene_t, index = "shannon") ~ Altitude + Distance, data = endem_df, family = gaussian())
print(summary(gam_shannon_aliene))
gam_shannon_autoctone <- gam(diversity(abb_autoctone_t, index = "shannon") ~ Altitude + Distance, data = endem_df, family = gaussian())
print(summary(gam_shannon_autoctone))

p_aliene <- ggplot(endem_df, aes(x = Altitude, y = diversity(abb_aliene_t, index = "shannon"))) +
  geom_point(color = 'darkgreen') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'green') +
  labs(title = "Specie Aliene: Shannon vs Altitude", x = "Altitude (m)", y = "Shannon Index")
p_autoctone <- ggplot(endem_df, aes(x = Altitude, y = diversity(abb_autoctone_t, index = "shannon"))) +
  geom_point(color = 'darkblue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'blue') +
  labs(title = "Specie Autoctone: Shannon vs Altitude", x = "Altitude (m)", y = "Shannon Index")
p_aliene + p_autoctone

#### 5. DIVERSITÀ GAMMA E ANALISI DEI LEPIDOTTERI ####

###### 5.1 Stima della Diversità Gamma ###### 
gamma_diversity_obs <- sum(rowSums(veg_matrix_subplot > 0) > 0)
cat("Gamma diversity osservata:", gamma_diversity_obs, "\n")
chao_estimates <- specpool(veg_matrix_subplot_t)
cat("Stima della ricchezza (Chao1, ACE, etc.):\n")
print(chao_estimates)

###### 5.2 Modellizzazione tramite GAM per la Vegetazione (Gamma Diversità) ###### 
#########  5.2.1 GAM per la Richness #########
gam_richness <- gam(Richness ~ s(Altitude, k = 30) + s(Distance), data = env_index_subplot)
cat("Risultati GAM per Richness:\n")
print(summary(gam_richness))
plot(gam_richness, residuals = TRUE, pch = 19, cex = 0.5,
     main = "GAM: Richness ~ s(Altitude, k=30) + s(Distance)")
#########  5.2.2 GAM per lo Shannon Index #########
gam_shannon_gamma <- gam(Shannon ~ s(Altitude) + s(Distance), data = env_index_subplot)
cat("Risultati GAM per Shannon (Gamma):\n")
print(summary(gam_shannon_gamma))
plot(gam_shannon_gamma, residuals = TRUE, pch = 19, cex = 0.5,
     main = "GAM: Shannon ~ s(Altitude) + s(Distance)")

###### 5.3 Profili di Diversità di Rényi ###### 
renyi_profile <- renyi(veg_matrix_subplot_t, scales = c(0, 0.5, 1, 1.5, 2, 2.5, 3, Inf))
cat("Profilo di Rényi:\n")
print(renyi_profile)
renyi_df <- as.data.frame(renyi_profile)
renyi_df$Sample <- rownames(renyi_df)
renyi_df <- pivot_longer(renyi_df, cols = -Sample, names_to = "scale", values_to = "diversity")
renyi_df <- renyi_df %>% filter(!is.na(as.numeric(scale))) %>% mutate(scale = as.numeric(scale))
renyi_plot <- ggplot(renyi_df, aes(x = scale, y = diversity, group = Sample, color = Sample)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "Profili di Rényi (tutti i subplot)",
       x = "Scale parameter (α)",
       y = "Diversity") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
print(renyi_plot)

### Raggruppamento per Altitudine e Confronto dei Profili
alt_breaks <- quantile(veg_altitude_subplot, probs = c(0, 0.33, 0.66, 1))
subplot_groups <- cut(veg_altitude_subplot, breaks = alt_breaks,
                      include.lowest = TRUE, labels = c("Low", "Medium", "High"))
renyi_profiles_group <- lapply(split(1:length(subplot_groups), subplot_groups), function(idx) {
  mat_group <- veg_matrix_subplot_t[idx, ]
  renyi(mat_group, scales = c(0, 0.5, 1, 1.5, 2, 2.5, 3, Inf))
})
renyi_group_list <- lapply(names(renyi_profiles_group), function(group) {
  rp <- as.data.frame(renyi_profiles_group[[group]])
  rp$Sample <- rownames(rp)
  rp <- pivot_longer(rp, cols = -Sample, names_to = "scale", values_to = "diversity")
  rp$scale <- as.numeric(rp$scale)
  rp$Group <- group
  return(rp)
})
renyi_group_df <- do.call(rbind, renyi_group_list)
renyi_group_plot <- ggplot(renyi_group_df, aes(x = scale, y = diversity,
                                               color = Group, group = interaction(Sample, Group))) +
  geom_line(alpha = 0.5) +
  geom_point() +
  labs(title = "Profili di Rényi per Gruppo di Altitudine",
       x = "Scale parameter (α)",
       y = "Diversity") +
  theme_minimal()
print(renyi_group_plot)

######  5.4 Analisi dei Lepidotteri ###### 
#########  5.4.1 Grafico: Lepidoptera Richness vs Altitude #########
lep_alt_plot <- ggplot(lep_df, aes(x = Altitude, y = Richness)) +
  geom_point(color = "slateblue1") +
  geom_smooth(method = "gam", se = FALSE, color = "slateblue3") +
  labs(title = "Lepidoptera Richness vs Altitude",
       x = "Altitude (m)", y = "Species Richness") +
  theme_minimal()
print(lep_alt_plot)
######### 5.4.2 Modello GAM: Richness ~ Altitude nei Lepidotteri #########
leprich_alt <- gam(Richness ~ Altitude, data = lep_df)
cat("Risultati GAM per Lepidotteri (Richness ~ Altitude):\n")
print(summary(leprich_alt))
######### 5.4.3 Grafico: Lepidoptera Richness vs Distance ###### 
lep_dist_plot <- ggplot(lep_df, aes(x = Distance, y = Richness)) +
  geom_point(color = "orange1") +
  geom_smooth(method = "gam", se = FALSE, color = "orange3") +
  labs(title = "Lepidoptera Richness vs Distance",
       x = "Distance (m)", y = "Species Richness") +
  theme_minimal()
print(lep_dist_plot)
####### 5.4.4 Modello GAM: Richness ~ Distance nei Lepidotteri #########
leprich_dist <- gam(Richness ~ Distance, data = lep_df)
cat("Risultati GAM per Lepidotteri (Richness ~ Distance):\n")
print(summary(leprich_dist))

###### 5.5 Relazione tra Vegetazione e Lepidotteri ###### 
# Calcolo degli indici di diversità per la vegetazione a livello di Plot
shannon_diversity_plot <- diversity(veg_matrix_plot_t, index = "shannon")
simpson_diversity_plot <- diversity(veg_matrix_plot_t, index = "invsimpson")
veg_richness_plot      <- colSums(veg_matrix_plot > 0)
total_abundance_plot   <- colSums(veg_matrix_plot_t)
margalef_plot          <- (veg_richness_plot - 1) / log(total_abundance_plot)
menhinick_plot         <- veg_richness_plot / sqrt(total_abundance_plot)
H0_plot <- veg_richness_plot         
H1_plot <- exp(shannon_diversity_plot) 
H2_plot <- simpson_diversity_plot     
aed_plot <- H0_plot + ((H1_plot^2)/(2*H2_plot))
# Regressione e correlazione tra AED e Lepidotteri
veg_lep_df <- data.frame(
  AED = aed_plot,
  LepRichness = lep_richness_plot
)
veg_lep_lm <- lm(LepRichness ~ AED, data = veg_lep_df)
cat("Risultati della regressione (LepRichness ~ AED):\n")
print(summary(veg_lep_lm))
cor_aed_lep <- cor.test(aed_plot, lep_richness_plot, method = "spearman", exact = FALSE)
cat("Correlazione Spearman tra AED e Lepidoptera Richness:\n")
print(cor_aed_lep)
