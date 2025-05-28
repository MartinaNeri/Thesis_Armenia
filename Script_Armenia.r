### Script per l'analisi della diversità (alfa, beta, gamma) e dei lepidotteri ###

#### 1. Caricamento delle librerie e definizione dei percorsi dei file ####
library(readxl)       # Lettura di file Excel
library(vegan)        # Calcolo degli indici di diversità
library(ggplot2)      # Creazione di grafici
library(stats)        # Funzioni statistiche di base
library(dplyr)        # Manipolazione dei dati
library(patchwork)    # Combinazione di grafici
library(adespatial)   # Analisi spaziali e profili di diversità
library(betapart)     # Analisi della diversità beta
library(mgcv)         # Modelli additivi generalizzati (GAM)
library(tidyr)        # Trasformazione dei dati in formato long
library(ggdendro)     # Creazione di dendrogrammi
library(indicspecies) # Analisi delle specie indicatrici (opzionale)
library(ggrepel)      # Etichettatura dei punti nei grafici

# Impostazione della directory di lavoro
setwd("C:/Users/Martina/Desktop/Tesi")

# Definizione dei percorsi dei file di input
file_path <- "C:/Users/Martina/Desktop/Tesi/Final_merged_file_Italian_and_Armenian_data.xlsx"
file_path_alien <- "C:/Users/Martina/Desktop/Tesi/Armenian_sp_alien.xlsx"

#### 2. Caricamento e pre-elaborazione dei dati ####

##### 2.1 Dati dei subplot (vegetazione) #####
veg_data_subplot <- read_excel(file_path, sheet = "spec_veg")
env_veg_subplot  <- read_excel(file_path, sheet = "env_veg_partial")
veg_altitude_subplot <- as.numeric(as.matrix(env_veg_subplot[1, 2:84]))
veg_street_d_subplot <- as.numeric(as.matrix(env_veg_subplot[3, 2:84]))
veg_h_s_distance_subplot <- as.numeric(as.matrix(env_veg_subplot[4, 2:84]))
veg_matrix_subplot   <- as.matrix(veg_data_subplot[,-1])
veg_matrix_subplot_t <- t(veg_matrix_subplot)  # Trasposizione: righe = siti, colonne = specie
veg_matrix_subplot_pa<- ifelse(veg_matrix_subplot_t > 0, 1, 0)  # Creazione della matrice presenza/assenza

##### 2.2 Dati a livello di plot (vegetazione) #####
veg_plot <- read_excel(file_path, sheet = "spec_veg_plot")
env_plot  <- read_excel(file_path, sheet = "env_veg_plot")
veg_altitude_plot <- as.numeric(as.matrix(env_plot[1, 2:15]))
veg_street_d_plot <- as.numeric(as.matrix(env_plot[3, 2:15]))
veg_h_s_distance_plot <- as.numeric(as.matrix(env_plot[4, 2:15]))
veg_matrix_plot   <- as.matrix(veg_plot[,-1])
veg_matrix_plot_t <- t(veg_matrix_plot)
veg_matrix_plot_pa<- ifelse(veg_matrix_plot_t > 0, 1, 0)

##### 2.3 Dati relativi ai lepidotteri #####
lep_data_plot <- read_excel(file_path, sheet = "spec_lep")
env_lep_plot  <- read_excel(file_path, sheet = "env_lep")
lep_altitude_plot <- as.numeric(as.matrix(env_lep_plot[1, 2:15]))
lep_dist_plot     <- as.numeric(as.matrix(env_lep_plot[2, 2:15]))
lep_h_s_distance_plot <- as.numeric(as.matrix(env_plot[4, 2:15]))
lep_matrix_plot   <- as.matrix(lep_data_plot[,-1])
lep_matrix_plot_t <- t(lep_matrix_plot)
lep_richness_plot <- colSums(lep_matrix_plot > 0)
lep_df <- as.data.frame(t(rbind(
  Altitude = lep_altitude_plot,
  Distance = lep_dist_plot,
  Richness = lep_richness_plot,
  H_S_Distance = lep_h_s_distance_plot
)))

##### 2.4 Dati relativi alle specie aliene e autoctone #####
sp_alien  <- read_excel(file_path_alien, sheet = "abb_aliene")
sp_autoct <- read_excel(file_path_alien, sheet = "abb_autoctone")
abb_aliene    <- as.matrix(sp_alien[,-1])
abb_autoctone <- as.matrix(sp_autoct[,-1])
abb_aliene_t    <- t(abb_aliene)
abb_autoctone_t <- t(abb_autoctone)
# Calcolo del rapporto tra specie aliene e totale
endem_ratio <- rowSums(abb_aliene_t) / (rowSums(abb_autoctone_t) + rowSums(abb_aliene_t))
endem_df <- data.frame(
  Altitude = veg_altitude_subplot,
  Distance = veg_street_d_subplot,
  HS_Distance = veg_h_s_distance_subplot,
  EndemRatio = endem_ratio
)

#### 3. Diversità Alfa ####

##### 3.1 Calcolo degli indici di diversità alfa #####
shannon_diversity_subplot <- diversity(veg_matrix_subplot_t, index = "shannon")
simpson_diversity_subplot <- diversity(veg_matrix_subplot_t, index = "invsimpson")
veg_richness_subplot      <- colSums(veg_matrix_subplot > 0)  # Numero di specie presenti
total_abundance_subplot   <- colSums(veg_matrix_subplot)
H0_subplot <- veg_richness_subplot         # Ricchezza (numero di specie)
H1_subplot <- exp(shannon_diversity_subplot)  # Diversità effettiva (indice di Shannon esponenziato)
H2_subplot <- simpson_diversity_subplot       
aed_subplot <- H0_subplot + ((H1_subplot^2)/(2*H2_subplot))  # Calcolo dell'indice composito AED

##### 3.2 Creazione del dataframe per l'analisi della diversità alfa #####
env_index_subplot <- data.frame(
  Altitude = veg_altitude_subplot,
  Distance = veg_street_d_subplot,
  Human_Settlement_Distance = veg_h_s_distance_subplot,
  Shannon  = shannon_diversity_subplot,
  Simpson  = simpson_diversity_subplot,
  AED      = aed_subplot,
  Richness = veg_richness_subplot
)

##### 3.3 Modellizzazione della relazione tra diversità alfa e variabili ambientali #####
# Modelli lineari generalizzati (GLM)
glm_shannon <- glm(Shannon ~ Altitude + Distance + Human_Settlement_Distance, data = env_index_subplot, family = gaussian())
glm_simpson <- glm(Simpson ~ Altitude + Distance + Human_Settlement_Distance, data = env_index_subplot, family = gaussian())
glm_aed     <- glm(AED ~ Altitude + Distance + Human_Settlement_Distance, data = env_index_subplot, family = gaussian())

# Modelli additivi generalizzati (GAM)
gam_shannon <- gam(Shannon ~ s(Altitude) + s(Distance) + s(Human_Settlement_Distance), data = env_index_subplot)
gam_simpson <- gam(Simpson ~ s(Altitude) + s(Distance) + s(Human_Settlement_Distance), data = env_index_subplot)
gam_aed     <- gam(AED ~ s(Altitude) + s(Distance) + s(Human_Settlement_Distance), data = env_index_subplot)

summary(glm_shannon)
summary(glm_simpson)
summary(glm_aed)
summary(gam_shannon)
summary(gam_simpson)
summary(gam_aed)

##### 3.4 Visualizzazione grafica della diversità alfa #####

###### 3.4.1 Grafici in funzione dell'altitudine ######
p1 <- ggplot(env_index_subplot, aes(x = Altitude, y = Shannon)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "a. Shannon vs Altitudine", x = "Altitudine (m)", y = "Indice di Shannon")

p2 <- ggplot(env_index_subplot, aes(x = Altitude, y = Simpson)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "b. Simpson vs Altitudine", x = "Altitudine (m)", y = "Indice di Simpson")

p3 <- ggplot(env_index_subplot, aes(x = Altitude, y = AED)) +
  geom_point(color = 'purple') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'pink') +
  labs(title = "c. AED vs Altitudine", x = "Altitudine (m)", y = "AED")

###### 3.4.2 Grafici in funzione della Dist. dalle strade ######
p4 <- ggplot(env_index_subplot, aes(x = Distance, y = Shannon)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "d. Shannon vs Dist. dalla strada", x = "Distanza (m)", y = "Indice di Shannon")

p5 <- ggplot(env_index_subplot, aes(x = Distance, y = Simpson)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "e. Simpson vs Dist. dalla strada", x = "Distanza (m)", y = "Indice di Simpson")

p6 <- ggplot(env_index_subplot, aes(x = Distance, y = AED)) +
  geom_point(color = 'purple') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'pink') +
  labs(title = "f. AED vs Dist. dalla strada", x = "Distance (m)", y = "AED")

###### 3.4.3 Grafici in funzione della distanza dagli insediamenti umani ######
p7 <- ggplot(env_index_subplot, aes(x = Human_Settlement_Distance, y = Shannon)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "g. Shannon vs Dist. dagli insediamenti", x = "Distanza (m)", y = "Indice di Shannon")

p8 <- ggplot(env_index_subplot, aes(x = Human_Settlement_Distance, y = Simpson)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "h. Simpson vs Dist. dagli insediamenti", x = "Distanza (m)", y = "Indice di Simpson")

p9 <- ggplot(env_index_subplot, aes(x = Human_Settlement_Distance, y = AED)) +
  geom_point(color = 'purple') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'pink') +
  labs(title = "i. AED vs Dist. dagli insediamenti", x = "Distanza (m)", y = "AED")

# Combinazione dei grafici mediante patchwork
combined_alpha <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)
png("combined_alpha.png", width = 12*300, height = 12*300, res = 300)
print(combined_alpha)
dev.off()

##### 3.5 Analisi della diversità per specie aliene e autoctone #####
shannon_aliene <- diversity(abb_aliene_t, index = "shannon")
simpson_aliene <- diversity(abb_aliene_t, index = "invsimpson")
shannon_autoctone <- diversity(abb_autoctone_t, index = "shannon")
simpson_autoctone <- diversity(abb_autoctone_t, index = "invsimpson")

diversity_df <- data.frame(
  Altitude = veg_altitude_subplot,
  Distance = veg_street_d_subplot,
  Human_Settlement_Distance = veg_h_s_distance_subplot,
  Shannon_Aliene = shannon_aliene,
  Simpson_Aliene = simpson_aliene,
  Shannon_Autoctone = shannon_autoctone,
  Simpson_Autoctone = simpson_autoctone
)
simpson_aliene <- ifelse(rowSums(abb_aliene_t) == 0, NA, diversity(abb_aliene_t, index = "invsimpson"))
diversity_df$Simpson_Aliene <- simpson_aliene
diversity_df <- na.omit(diversity_df)

summary(diversity_df$Simpson_Aliene)

###### 3.5.1 Analisi statistica e grafica per specie aliene e autoctone ######

# Modellizzazione GAM per specie aliene e autoctone
#Shannon
gam_shannon_aliene    <- gam(diversity(abb_aliene_t, index = "shannon") ~ Altitude + Distance + Human_Settlement_Distance, 
                             data = env_index_subplot, family = gaussian())
print(summary(gam_shannon_aliene))
gam_shannon_autoctone <- gam(diversity(abb_autoctone_t, index = "shannon") ~ Altitude + Distance + Human_Settlement_Distance, 
                             data = env_index_subplot, family = gaussian())
print(summary(gam_shannon_autoctone))

#Simpson
gam_simpson_aliene    <- gam(diversity(abb_aliene_t, index = "simpson") ~ Altitude + Distance + Human_Settlement_Distance, 
                             data = env_index_subplot, family = gaussian())
print(summary(gam_simpson_aliene))
gam_simpson_autoctone <- gam(diversity(abb_autoctone_t, index = "simpson") ~ Altitude + Distance + Human_Settlement_Distance, 
                             data = env_index_subplot, family = gaussian())
print(summary(gam_simpson_autoctone))

# Grafici per specie aliene
paa1 <- ggplot(diversity_df, aes(x = Altitude, y = Shannon_Aliene)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "a. Shannon Al. vs Altitudine", x = "Altitudine (m)", y = "Indice di Shannon")

paa2 <- ggplot(diversity_df, aes(x = Distance, y = Shannon_Aliene)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "b. Shannon Al. vs Dist. dalle strade", x = "Distanza (m)", y = "Indice di Shannon")

paa3 <- ggplot(diversity_df, aes(x = Human_Settlement_Distance, y = Shannon_Aliene)) +
  geom_point(color = 'green') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkgreen') +
  labs(title = "c. Shannon Al. vs Dist. Insediamenti", x = "Distanza (m)", y = "Indice di Shannon")

# Grafici per specie autoctone
paa4 <- ggplot(diversity_df, aes(x = Altitude, y = Shannon_Autoctone)) +
  geom_point(color = 'blue') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkblue') +
  labs(title = "d. Shannon Aut. vs Altitudine", x = "Altitudine (m)", y = "Indice di Shannon")

paa5 <- ggplot(diversity_df, aes(x = Distance, y = Shannon_Autoctone)) +
  geom_point(color = 'red') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkred') +
  labs(title = "e. Shannon Aut. vs Dist. dalla strada", x = "Distanza (m)", y = "Indice di Shannon")

paa6 <- ggplot(diversity_df, aes(x = Human_Settlement_Distance, y = Shannon_Autoctone)) +
  geom_point(color = 'green') +
  stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE, color = 'darkgreen') +
  labs(title = "f. Shannon Aut. vs Dist. Insediamenti", x = "Distanza (m)", y = "Indice di Shannon")

combined_endem <- paa1 + paa2 + paa3 + paa4 + paa5 + paa6
png("combined_endem.png", width = 12*300, height = 5*300, res = 300)
print(combined_endem)
dev.off()

#### 4. Diversità Beta ####

##### 4.1 Calcolo delle dissimilarità tra siti #####
veg_bray_beta <- vegdist(veg_matrix_plot_t, method = "bray")
dissim_sorensen <- vegdist(veg_matrix_subplot_pa, method = "bray", binary = TRUE)
dissim_jaccard   <- vegdist(veg_matrix_subplot_pa, method = "jaccard")

print(dissim_sorensen)

##### 4.2 Decomposizione della diversità beta #####
beta_pair       <- beta.pair(veg_matrix_subplot_pa, index.family = "sorensen")
beta_turnover   <- beta_pair$beta.sim      
beta_nestedness <- beta_pair$beta.sne      
beta_total      <- beta_pair$beta.sor   
beta_multi_sorensen <- beta.multi(veg_matrix_subplot_pa, index.family = "sorensen")
print(beta_multi_sorensen)

##### 4.3 Analisi PERMANOVA #####
permanova_bray_beta   <- adonis2(veg_matrix_plot_t ~ veg_altitude_plot + veg_h_s_distance_plot, method = "bray")
print(permanova_bray_beta)
permanova_jaccard_beta<- adonis2(veg_matrix_plot_pa ~ veg_altitude_plot + veg_h_s_distance_plot, method = "jaccard")
print(permanova_jaccard_beta)

#### 4.4 NMDS per visualizzare le differenze nella composizione delle comunità ####
nmds_beta <- metaMDS(veg_bray_beta, k = 2, trymax = 100)
nmds_scores_beta <- as.data.frame(scores(nmds_beta))
nmds_scores_beta$Sample <- rownames(nmds_scores_beta)
nmds_plot_beta <- ggplot(nmds_scores_beta, aes(x = NMDS1, y = NMDS2, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5) +
  labs(title = "NMDS: Composizione delle Comunità", x = "NMDS1", y = "NMDS2") +
  theme_minimal()
png("nmds_plot_beta.png", width = 8*300, height = 6*300, res = 300)
print(nmds_plot_beta)
dev.off()

#### 4.5 Dendrogramma basato sulla dissimilarità di Bray-Curtis ####
bray_clustering_beta <- hclust(veg_bray_beta, method = "average")
dendro_data_beta <- as.dendrogram(bray_clustering_beta)
dendro_plot_beta <- ggdendrogram(dendro_data_beta, theme_dendro = FALSE) +
  ggtitle("Dendrogramma: Bray-Curtis") +
  theme_minimal()
png("dendrogramma.png", width = 8*300, height = 6*300, res = 300)
print(dendro_plot_beta)
dev.off()

#### 4.6 Analisi pairwise e relazione con variabili ambientali ####
beta_bray_pairwise <- vegdist(veg_matrix_plot_t, method = "bray")
beta_bray_df <- as.data.frame(as.matrix(beta_bray_pairwise))
beta_bray_df$Site1 <- rownames(beta_bray_df)
beta_bray_long <- beta_bray_df %>%
  pivot_longer(cols = -Site1, names_to = "Site2", values_to = "Beta_Diversity") %>%
  filter(Site1 != Site2)

# Assegnazione dei nomi dei siti
site_names <- veg_plot[[1]]
rownames(veg_matrix_plot) <- site_names
colnames(veg_matrix_plot_t) <- site_names

env_data <- data.frame(
  Site = colnames(veg_matrix_plot),
  Altitude = veg_altitude_plot,
  Street_Distance = veg_street_d_plot,
  Human_Settlement_Distance = veg_h_s_distance_plot
)

beta_bray_long <- beta_bray_long %>%
  left_join(env_data, by = c("Site1" = "Site")) %>%
  rename(Altitude1 = Altitude, Street_Distance1 = Street_Distance, Human_Settlement_Distance1 = Human_Settlement_Distance) %>%
  left_join(env_data, by = c("Site2" = "Site")) %>%
  rename(Altitude2 = Altitude, Street_Distance2 = Street_Distance, Human_Settlement_Distance2 = Human_Settlement_Distance) %>%
  mutate(
    Altitude_Diff = abs(Altitude1 - Altitude2),
    Street_Distance_Diff = abs(Street_Distance1 - Street_Distance2),
    Human_Settlement_Distance_Diff = abs(Human_Settlement_Distance1 - Human_Settlement_Distance2)
  )
beta_bray_long_clean <- beta_bray_long %>%
  mutate(Pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = "_")) %>%
  distinct(Pair, .keep_all = TRUE) %>%
  select(-Pair)

plot_altitude_beta <- ggplot(beta_bray_long_clean, aes(x = Altitude_Diff, y = Beta_Diversity)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "gam", se = TRUE, color = "black") +
  labs(title = "a. Diff. di Altitudine vs Beta Diversità",
       x = "Differenza di Altitudine",
       y = "Beta Diversità (Bray-Curtis)") +
  theme_minimal()

summary(gam(Beta_Diversity ~ Altitude_Diff, data = beta_bray_long_clean))

plot_street_beta <- ggplot(beta_bray_long_clean, aes(x = Street_Distance_Diff, y = Beta_Diversity)) +
  geom_point(alpha = 0.5, color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "b. Diff. di Dist. dalla strada vs Beta Diversità",
       x = "Differenza di Distanza",
       y = "Beta Diversità (Bray-Curtis)") +
  theme_minimal()

summary(gam(Beta_Diversity ~ Street_Distance_Diff, data = beta_bray_long_clean))

plot_hs_beta <- ggplot(beta_bray_long_clean, aes(x = Human_Settlement_Distance_Diff, y = Beta_Diversity)) +
  geom_point(alpha = 0.5, color = "green") +
  geom_smooth(method = "gam", se = TRUE, color = "black") +
  labs(title = "c. Diff. di Dist. dagli Insediamenti vs Beta Diversità",
       x = "Differenza di Distanza",
       y = "Beta Diversità (Bray-Curtis)") +
  theme_minimal()

summary(gam(Beta_Diversity ~ Human_Settlement_Distance_Diff, data = beta_bray_long_clean))

#stampo in png tutti e tre i grafici insieme
combined_beta <- plot_altitude_beta + plot_street_beta + plot_hs_beta
png("combined_beta.png", width = 14*300, height = 5*300, res = 300)
print(combined_beta)
dev.off()

#### 5. Diversità Gamma e analisi dei lepidotteri ####
##### 5.1 Stima della diversità gamma #####
gamma_diversity_obs <- sum(rowSums(veg_matrix_subplot > 0) > 0)
chao_estimates <- specpool(veg_matrix_subplot_t)
print(chao_estimates)

##### 5.1 Profili di diversità di Rényi #####
renyi_profile <- renyi(veg_matrix_plot_t, scales = c(0, 0.5, 1, 1.5, 2, 2.5, 3, Inf))
print(renyi_profile)
renyi_df <- as.data.frame(renyi_profile)
renyi_df$Sample <- rownames(renyi_df)
renyi_df <- pivot_longer(renyi_df, cols = -Sample, names_to = "scale", values_to = "diversity")
renyi_df <- renyi_df %>% filter(!is.na(as.numeric(scale))) %>% mutate(scale = as.numeric(scale))

# Selezione del primo punto di ogni profilo
first_points <- renyi_df %>%
  group_by(Sample) %>%
  slice(1) %>%  # Seleziona il primo punto
  ungroup()

# Creazione del grafico con le etichette posizionate all'inizio
renyi_plot <- ggplot(renyi_df, aes(x = scale, y = diversity, group = Sample, color = Sample)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_label_repel(data = first_points, aes(label = Sample), size = 3, box.padding = 0.3, point.padding = 0.2, show.legend = FALSE) + 
  labs(
    title = "Profili di Rényi (tutti i plot)",
    x = "Parametro scala (α)",
    y = "Diversity"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Salvataggio del grafico
png("renyi_plot_with_start_labels.png", width = 12*300, height = 8*300, res = 300)
print(renyi_plot)
dev.off()


# Raggruppamento per Altitudine e Confronto dei Profili
alt_breaks <- quantile(veg_altitude_plot, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE)
plot_groups <- cut(veg_altitude_plot, breaks = alt_breaks, include.lowest = TRUE, labels = c("Low", "Medium", "High"))

# Verifica che gli indici siano validi prima di accedere alla matrice
renyi_profiles_group <- lapply(split(1:nrow(veg_matrix_plot_t), plot_groups), function(idx) {
  if (all(idx <= nrow(veg_matrix_plot_t))) {
    mat_group <- veg_matrix_plot_t[idx, , drop = FALSE]
    return(renyi(mat_group, scales = c(0, 0.5, 1, 1.5, 2, 2.5, 3, Inf)))
  }
})

# Rimuovi gruppi nulli
renyi_profiles_group <- Filter(Negate(is.null), renyi_profiles_group)

# Creazione del dataframe per i profili di Renyi
renyi_group_list <- lapply(names(renyi_profiles_group), function(group) {
  rp <- as.data.frame(renyi_profiles_group[[group]])
  rp$Sample <- rownames(rp)
  rp <- pivot_longer(rp, cols = -Sample, names_to = "scale", values_to = "diversity")
  rp$scale <- as.numeric(rp$scale)
  rp$Group <- group
  return(rp)
})

renyi_group_df <- do.call(rbind, renyi_group_list)
renyi_group_plot <- ggplot(renyi_group_df, aes(x = scale, y = diversity, color = Group, group = interaction(Sample, Group))) +
  geom_line(alpha = 0.5) +
  geom_point() +
  labs(title = "Profili di Rényi per Gruppo di Altitudine", x = "Parametro scala (α)", y = "Diversità") +
  theme_minimal()
png("renyi_group_plot.png", width = 10*300, height = 8*300, res = 300)
print(renyi_group_plot)
dev.off()

#### 6. Analisi dei lepidotteri ####
lep_alt_plot <- ggplot(lep_df, aes(x = Altitude, y = Richness)) +
  geom_point(color = "slateblue1") +
  geom_smooth(method = "gam", se = FALSE, color = "slateblue3") +
  labs(title = "a. Ricch. lepidotteri vs Altitudine",
       x = "Altitudine (m)", y = "Ricchezza di specie") +
  theme_minimal()

leprich_alt <- gam(Richness ~ Altitude, data = lep_df)
print(summary(leprich_alt))

lep_dist_plot <- ggplot(lep_df, aes(x = Distance, y = Richness)) +
  geom_point(color = "orange1") +
  geom_smooth(method = "gam", se = FALSE, color = "orange3") +
  labs(title = "b. Ricch. lepidotteri vs Dist. dalla strada",
       x = "Dist. dalla strada (m)", y = "Ricchezza di specie") +
  theme_minimal()

leprich_dist <- gam(Richness ~ Distance, data = lep_df)
print(summary(leprich_dist))

lep_hs_plot <- ggplot(lep_df, aes(x = H_S_Distance, y = Richness)) +
  geom_point(color = "green1") +
  geom_smooth(method = "gam", se = FALSE, color = "green3") +
  labs(title = "c. Ricch. lepidotteri vs Dist. dagli insediamenti",
       x = "Dist. dagli insediamenti (m)", y = "Ricchezza di specie") +
  theme_minimal()

leprich_hs <- gam(Richness ~ H_S_Distance, data = lep_df)
print(summary(leprich_hs))

# Combinazione dei grafici
combined_lep <- lep_alt_plot + lep_dist_plot + lep_hs_plot
png("combined_lep.png", width = 12*300, height = 5*300, res = 300)
print(combined_lep)
dev.off()


###### 6.1 Relazione tra vegetazione e lepidotteri ########
shannon_diversity_plot <- diversity(veg_matrix_plot_t, index = "shannon")
simpson_diversity_plot <- diversity(veg_matrix_plot_t, index = "invsimpson")
veg_richness_plot      <- colSums(veg_matrix_plot > 0)
total_abundance_plot   <- colSums(veg_matrix_plot_t)
H0_plot <- veg_richness_plot         
H1_plot <- exp(shannon_diversity_plot) 
H2_plot <- simpson_diversity_plot     
aed_plot <- H0_plot + ((H1_plot^2)/(2*H2_plot))
veg_lep_df <- data.frame(
  AED = aed_plot,
  LepRichness = lep_richness_plot
)
veg_lep_lm <- lm(LepRichness ~ AED, data = veg_lep_df)
print(summary(veg_lep_lm))

veg_lep_plot <- ggplot(veg_lep_df, aes(x = AED, y = LepRichness)) +
  geom_point(color = "purple") +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  labs(title = "Relazione tra AED e Ricchezza di Lepidotteri",
       x = "AED", y = "Ricchezza di Lepidotteri") +
  theme_minimal()
print(veg_lep_plot)

veg_lep_gam <- gam(LepRichness ~ s(AED), data = veg_lep_df)
print(summary(veg_lep_gam))

veg_lep_plot_gam <- ggplot(veg_lep_df, aes(x = AED, y = LepRichness)) +
  geom_point(color = "purple") +
  geom_smooth(method = "gam", se = FALSE, color = "purple") +
  labs(title = "Relazione tra AED e Ricchezza di Lepidotteri (GAM)",
       x = "AED", y = "Ricchezza di Lepidotteri") +
  theme_minimal()
print(veg_lep_plot_gam)

cor_aed_lep <- cor.test(aed_plot, lep_richness_plot, method = "spearman", exact = FALSE)
print(cor_aed_lep)
