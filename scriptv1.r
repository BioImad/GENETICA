library(ggplot2)
library(magrittr)
library(dplyr)


# Leggi il file specificando il percorso corretto
roh_indiv <- read.table("D:/plink_win64_20250615nuovo/IBS_MXL_CHB_RoH.hom.indiv", header = TRUE)

# Crea il boxplot delle Run of Homozygosity (RoH) per popolazione
p_roh <- ggplot(roh_indiv, aes(x = FID, y = KB, fill = FID)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Run of Homozygosity per Popolazione",
       x = "Popolazione", y = "Lunghezza totale (KB)") +
  theme(legend.position = "none")

print(p_roh)

data_dir <- "D:/plink_win64_20250615nuovo/"

pca <- read.table(file.path(data_dir, "IBS_MXL_CHB_PCA.eigenvec"), header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:8))

# Converti FID in fattore con ordine specifico
pca$FID <- factor(pca$FID, levels = c("IBS", "MXL", "CHB"))

# Palette colori personalizzata
colors <- c("IBS" = "#1b9e77", "MXL" = "#d95f02", "CHB" = "#7570b3")

# Plot PC1 vs PC2
p1 <- ggplot(pca, aes(x = PC1, y = PC2, color = FID)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC2 - IBS, MXL, CHB",
       x = "PC1", y = "PC2", color = "Popolazione") +
  theme(legend.position = "bottom")

print(p1)

# Plot PC1 vs PC3
p2 <- ggplot(pca, aes(x = PC1, y = PC3, color = FID)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC3 - IBS, MXL, CHB",
       x = "PC1", y = "PC3", color = "Popolazione") +
  theme(legend.position = "bottom")

print(p2)

# Plot PC2 vs PC3
p3 <- ggplot(pca, aes(x = PC2, y = PC3, color = FID)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = "PCA: PC2 vs PC3 - IBS, MXL, CHB",
       x = "PC2", y = "PC3", color = "Popolazione") +
  theme(legend.position = "bottom")

print(p3)

# PCA interattiva con plotly (opzionale)
library(plotly)

pca_plotly <- plot_ly(pca, x = ~PC1, y = ~PC2, color = ~FID, colors = colors,
                      type = 'scatter', mode = 'markers')

pca_plotly <- layout(pca_plotly, title = list(text = "PCA Interattiva: PC1 vs PC2"))

print(pca_plotly)


# ===========================
# 2. ANALISI DEL MISSINGNESS
# ===========================

missingness <- read.table(file.path(data_dir, "1KGPReadyBioInfo2024.imiss"), header = TRUE)

# Filtra per popolazioni
miss_subset <- missingness %>% filter(FID %in% c("IBS", "MXL", "CHB"))

# Statistiche descrittive
miss_stats <- miss_subset %>%
  group_by(FID) %>%
  summarise(
    n_individuals = n(),
    mean_miss = mean(F_MISS),
    median_miss = median(F_MISS),
    sd_miss = sd(F_MISS),
    min_miss = min(F_MISS),
    max_miss = max(F_MISS)
  )
print("Statistiche Missingness per Popolazione:")
print(miss_stats)

# Boxplot missingness
p_miss <- ggplot(miss_subset, aes(x = FID, y = F_MISS, fill = FID)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribuzione del Missingness per Popolazione",
       x = "Popolazione", y = "Proporzione Missing Data") +
  theme(legend.position = "none")

print(p_miss)

# ===========================
# 3. ANALISI DELL'ETEROZIGOSITÀ
# ===========================

het_data <- read.table(file.path(data_dir, "IBS_MXL_CHB_inbreeding.het"), header = TRUE)
het_data$FID <- factor(het_data$FID, levels = c("IBS", "MXL", "CHB"))

het_stats <- het_data %>%
  group_by(FID) %>%
  summarise(
    n_individuals = n(),
    mean_F = mean(F),
    median_F = median(F),
    sd_F = sd(F),
    min_F = min(F),
    max_F = max(F)
  )
print("Statistiche Coefficiente di Inbreeding (F) per Popolazione:")
print(het_stats)

# Boxplot coefficiente F
p_het <- ggplot(het_data, aes(x = FID, y = F, fill = FID)) +
  geom_boxplot() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Coefficiente di Inbreeding per Popolazione",
       x = "Popolazione", y = "Coefficiente F") +
  theme(legend.position = "none")

print(p_het)

# ===========================
# 4. CALCOLO DELLA DIFFERENZA DI FREQUENZA (DAF) E MANHATTAN PLOT
# ===========================

freq <- read.table(file.path(data_dir, "IBS_MXL_CHB_freq.frq.strat"), header = TRUE)

IBS <- freq %>% filter(CLST == "IBS")
MXL <- freq %>% filter(CLST == "MXL")
CHB <- freq %>% filter(CLST == "CHB")

# Estrazione posizione SNP da nome (assumendo formato chr_pos)
positions <- as.integer(sapply(strsplit(as.character(IBS$SNP), "_"), `[`, 2))

# Calcola delta frequenze per coppie
delta_IBS_MXL <- abs(IBS$MAF - MXL$MAF)
delta_IBS_CHB <- abs(IBS$MAF - CHB$MAF)
delta_MXL_CHB <- abs(MXL$MAF - CHB$MAF)

# Crea dataframe per Manhattan plot
deltaFreq_IBS_MXL <- data.frame(chr = IBS$CHR, pos = positions, DAF = delta_IBS_MXL)
deltaFreq_IBS_CHB <- data.frame(chr = IBS$CHR, pos = positions, DAF = delta_IBS_CHB)
deltaFreq_MXL_CHB <- data.frame(chr = IBS$CHR, pos = positions, DAF = delta_MXL_CHB)

# Converti chr in fattore per colore alternato
deltaFreq_IBS_MXL$chr_factor <- factor(deltaFreq_IBS_MXL$chr)
deltaFreq_IBS_CHB$chr_factor <- factor(deltaFreq_IBS_CHB$chr)
deltaFreq_MXL_CHB$chr_factor <- factor(deltaFreq_MXL_CHB$chr)

# Funzione per Manhattan plot migliorato
plot_manhattan <- function(df, title) {
  ggplot(df, aes(x = pos, y = DAF, color = chr_factor)) +
    geom_point(alpha = 0.6, size = 0.7) +
    facet_grid(. ~ chr_factor, scales = "free_x", space = "free_x") +
    scale_color_manual(values = rep(c("gray30", "gray70"), length.out = length(levels(df$chr_factor)))) +
    theme_minimal() +
    labs(title = title, x = "Posizione", y = "Delta Frequency") +
    theme(axis.text.x = element_blank(), legend.position = "none")
}

p4_IBS_MXL <- plot_manhattan(deltaFreq_IBS_MXL, "Manhattan Plot: Differenze Frequenze IBS vs MXL")
p4_IBS_CHB <- plot_manhattan(deltaFreq_IBS_CHB, "Manhattan Plot: Differenze Frequenze IBS vs CHB")
p4_MXL_CHB <- plot_manhattan(deltaFreq_MXL_CHB, "Manhattan Plot: Differenze Frequenze MXL vs CHB")

print(p4_IBS_MXL)
print(p4_IBS_CHB)
print(p4_MXL_CHB)

# Top 10 SNP con maggiori differenze per confronto
top_IBS_MXL <- deltaFreq_IBS_MXL %>% slice_max(DAF, n = 10)
top_IBS_CHB <- deltaFreq_IBS_CHB %>% slice_max(DAF, n = 10)
top_MXL_CHB <- deltaFreq_MXL_CHB %>% slice_max(DAF, n = 10)

print("Top 10 SNP con maggiori differenze IBS vs MXL:")
print(top_IBS_MXL)

print("Top 10 SNP con maggiori differenze IBS vs CHB:")
print(top_IBS_CHB)

print("Top 10 SNP con maggiori differenze MXL vs CHB:")
print(top_MXL_CHB)

# ===========================
# 5. CORRELAZIONE TRA FREQUENZE ALLELICHE
# ===========================

# Scatter plot con base R
par(mfrow = c(2, 2))

smoothScatter(IBS$MAF, MXL$MAF,
              main = "Correlazione Frequenze: IBS vs MXL",
              xlab = "MAF IBS", ylab = "MAF MXL")

smoothScatter(IBS$MAF, CHB$MAF,
              main = "Correlazione Frequenze: IBS vs CHB",
              xlab = "MAF IBS", ylab = "MAF CHB")

smoothScatter(MXL$MAF, CHB$MAF,
              main = "Correlazione Frequenze: MXL vs CHB",
              xlab = "MAF MXL", ylab = "MAF CHB")

par(mfrow = c(1, 1))

# Calcola e stampa correlazioni
cor_IBS_MXL <- cor(IBS$MAF, MXL$MAF)
cor_IBS_CHB <- cor(IBS$MAF, CHB$MAF)
cor_MXL_CHB <- cor(MXL$MAF, CHB$MAF)

cat("Correlazioni tra frequenze alleliche:\n")
cat(sprintf("IBS vs MXL: %.3f\n", cor_IBS_MXL))
cat(sprintf("IBS vs CHB: %.3f\n", cor_IBS_CHB))
cat(sprintf("MXL vs CHB: %.3f\n", cor_MXL_CHB))

# ===========================
# 6. SALVATAGGIO GRAFICI
# ===========================

ggsave(file.path(data_dir, "PCA_PC1_PC2_IBS_MXL_CHB.png"), plot = p1, width = 10, height = 8)
ggsave(file.path(data_dir, "PCA_PC1_PC3_IBS_MXL_CHB.png"), plot = p2, width = 10, height = 8)
ggsave(file.path(data_dir, "Missingness_boxplot.png"), plot = p_miss, width = 8, height = 6)
ggsave(file.path(data_dir, "Inbreeding_boxplot.png"), plot = p_het, width = 8, height = 6)
ggsave(file.path(data_dir, "Manhattan_IBS_MXL.png"), plot = p4_IBS_MXL, width = 15, height = 6)
ggsave(file.path(data_dir, "Manhattan_IBS_CHB.png"), plot = p4_IBS_CHB, width = 15, height = 6)
ggsave(file.path(data_dir, "Manhattan_MXL_CHB.png"), plot = p4_MXL_CHB, width = 15, height = 6)

