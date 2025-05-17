### Author: Annie Squiavinato  
### date: "2025-04-25"  

### Directory for the new library 24-06-2024  
### or directory accessible on all nodes  
mylib <- "/data04/projects04/PatriciaPossik/pdx_la_sanger/lib/Rpackages/4.4.1"  
libdir <- "/data04/tools/R/Rpackages/4.4.1"  

### Command to apply the libPaths selection  
.libPaths(c(mylib, libdir))  
.libPaths()  

# setting directory  
setwd("~/pdx_la_sanger/results/INCA_res")  
### set CRAN Mirror https://cran-r.c3sl.ufpr.br/  
options(repos = c("CRAN" = "https://cran-r.c3sl.ufpr.br/"))  

### Directory for the new library 24-06-2024  
libdir <- "/scr/R/Rpackages/4.4.1"  

##Dependencies for this project are recorded and managed with renv
## To setup your Rstudio environment you might need to run:
renv::restore()
#renv::init() # To initialise revn
#renv::restore()
#renv::init() # To initialise revn
#renv::install("glue") # To add a package from CRAN
#renv::install("BiocManager") # To add a package from Bioconductor
#renv::install("bioC::ggplot2") # To add a package from Bioconductor
#renv::snapshot() # To record any packages you've used and their versions and sources into the renv.lock file
#renv::restore() # To rebuild an environment from the renv.lockfile
# loading packages

library(arsenal)
library(data.table)
library(scales)
library(ggpubr)
library(gridExtra)
library(ggsci)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(data.table) # to put plots side by side
library(dunn.test)
library(here)
renv::snapshot()
# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/02_Figure_S3b_c.R")

# hailstorms_chr11q, 5p, and 22q, one sample per patient, only passed QC sequenza
# Analisys of association FCS and BCS with presence or absence of hailstorm


# table describing samples that have or not hailstorms
Table_hailstorms <- read.csv("./data/allSamples_hailstorm_regions_corrected.csv", header = TRUE, sep = "\t")

### Transposing the dataframe and maitaining colnames and rownames
transpose_data_hail <- setNames(data.frame(t(Table_hailstorms[, -1])), Table_hailstorms[, 1])
transpose_data_hail$Tumor_Sample_Barcode <- NA
transpose_data_hail$Tumor_Sample_Barcode <- rownames(transpose_data_hail)
# Open the table
Metadata <- read.csv("./data/Table_metadata_REVIEWED_18_01_25.csv")

#List of samples to be excluded, cellularity 0.35>x (QC done in Jan/Feb 2025, Sequenza data)
samples <- c(
  "PD53330c", "PD53332d", "PD53333f", "PD53334a", "PD53338a", "PD53342a", "PD53343h", "PD53345a", "PD53347c",
  "PD53348a", "PD53349d", "PD53350d", "PD53351a", "PD53352d", "PD53354a", "PD53355c", "PD53357c", "PD53358a",
  "PD53359c", "PD53361a", "PD53362a", "PD53364c", "PD53365a", "PD53368a", "PD53369a"
)

# Filtrating data
Table_QC_fev_2025 <- Metadata %>% filter(Tumor_Sample_Barcode %in% samples)

# merging the CNApp data
Table_merged_ed <- merge(transpose_data_hail, Table_QC_fev_2025, by = "Tumor_Sample_Barcode")
print(Table_merged_ed)

####### Checking the association with presence of hailstorms and:
# to FCS values
# I gonna analyse all the alterations of chrs, so I will need to applya multiple test correction to a total number of tests performed (each alteration):
colnames(Table_merged_ed)

# Vector of chr regions that have two levels (yes or no), excluding samples that only have one level 11p, 15q, 16q, 10, and samples that have a presence
# of hailstorm in one sample (4p, 6q, 1p, 8p, 17q, 3q, 7q, 5q, 1q, 12p, 20q, 17p)
# Vector of chr
hailstorms_chr <- c("5p", "6p", "11q", "12q", "19q", "21q", "7p", "22q", "4q", "10p")

# Vectors to store raw and adjusted p-values
p_values_fcs <- c()
p_values_bcs <- c()

# Perform Wilcoxon tests for FCS
for (chr in hailstorms_chr) {
  p_value_fcs <- wilcox.test(Table_merged_ed$FCS ~ get(chr, Table_merged_ed), data = Table_merged_ed, exact = FALSE)$p.value
  p_values_fcs <- c(p_values_fcs, p_value_fcs)
}

# Aplicar a correção de Benjamini-Hochberg (FDR) para FCS
adjusted_p_values_fcs <- p.adjust(p_values_fcs, method = "BH")

# Realizar testes de Wilcoxon para BCS
for (chr in hailstorms_chr) {
  p_value_bcs <- wilcox.test(Table_merged_ed$BCS ~ get(chr, Table_merged_ed), data = Table_merged_ed, exact = FALSE)$p.value
  p_values_bcs <- c(p_values_bcs, p_value_bcs)
}

# Apply Benjamini-Hochberg Correction (FDR) to BCS
adjusted_p_values_bcs <- p.adjust(p_values_bcs, method = "BH")

# Create a data.frame with the results
results_table <- data.frame(
  Regiao = hailstorms_chr,
  P_Value_FCS = p_values_fcs,
  Adjusted_P_Value_FCS = adjusted_p_values_fcs,
  P_Value_BCS = p_values_bcs,
  Adjusted_P_Value_BCS = adjusted_p_values_bcs
)

# print the table
print(results_table)

# Save the table to a CSV file 
write.csv(
  results_table,
  here("results", "FCS_BCS_vs_hails_wilcoxon.csv"),
  row.names = TRUE,
  sep = ","
)
# print table
print(results_table)

#### plotting graphics #####
################## FCS
# Find the index of the region "11q" in the results table
indice_11q <- which(results_table$Regiao == "11q")

# Access adjusted p-values for "5p"
adjusted_p_value_fcs_11q <- results_table$Adjusted_P_Value_FCS[indice_11q]

# Plot FCS
p2 <- ggplot(Table_merged_ed, aes(x = `11q`, y = FCS, fill = `11q`)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Hailstorm 11q", y = "FCS") +
  scale_fill_npg() +
  coord_cartesian(ylim = c(1, 1000), expand = TRUE) +
  annotate("text",
    x = 1.5, y = max(Table_merged_ed$FCS) * 1.05, # Adicionar valor-p ajustado manualmente
    label = paste("p.adj =", round(adjusted_p_value_fcs_11q, 4)), color = "black", size = 4
  ) +
  theme_classic()

p2

# saving
ggsave(
  filename = here("results", "FCS_hailstorm_11q_plot.pdf"),
  plot = p2,
  device = "pdf",
  width = 5.52,
  height = 5.7
)

# Find the index of the region "22q" in the results table
indice_22q <- which(results_table$Regiao == "22q")

# Access adjusted p-values for "5p"
adjusted_p_value_fcs_22q <- results_table$Adjusted_P_Value_FCS[indice_22q]

# Plot FCS
p3 <- ggplot(Table_merged_ed, aes(x = `22q`, y = FCS, fill = `22q`)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Hailstorm 22q", y = "FCS") +
  scale_fill_npg() +
  coord_cartesian(ylim = c(1, 1000), expand = TRUE) +
  annotate("text",
    x = 1.5, y = max(Table_merged_ed$FCS) * 1.05, # Add pvalue manual
    label = paste("p.adj =", round(adjusted_p_value_fcs_22q, 4)), color = "black", size = 4
  ) +
  theme_classic()
p3

ggsave(
  filename = here("results", "FCS_hailstorm_22q_plot.pdf"),
  plot = p3,
  device = "pdf",
  width = 5.52,
  height = 5.7
)
# Find the index of the region "5p" in the results table
indice_5p <- which(results_table$Regiao == "5p")

# Access adjusted p-values for "5p"
adjusted_p_value_fcs_5p <- results_table$Adjusted_P_Value_FCS[indice_5p]


p4 <- ggplot(Table_merged_ed, aes(x = `5p`, y = FCS, fill = `5p`)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Hailstorm 5p", y = "FCS") +
  scale_fill_npg() +
  coord_cartesian(ylim = c(1, 1000), expand = TRUE) +
  annotate("text",
    x = 1.5, y = max(Table_merged_ed$FCS) * 1.05, # Add pvalue manual
    label = paste("p.adj =", round(adjusted_p_value_fcs_5p, 4)), color = "black", size = 4
  ) +
  theme_classic()
p4
ggsave(
  filename = here("results", "FCS_hailstorm_5p_plot.pdf"),
  plot = p4,
  device = "pdf",
  width = 5.52,
  height = 5.7
)
combined_plot <- grid.arrange(p2, p3, p4, ncol = 3)


################## BCS

# Encontrar o índice da região "11q" na tabela de resultados
indice_11q <- which(results_table$Regiao == "11q")

# Acessar os valores-p ajustados para 11q
adjusted_p_value_bcs_11q <- results_table$Adjusted_P_Value_BCS[indice_11q]

# Gráfico para FCS
p5 <- ggplot(Table_merged_ed, aes(x = `11q`, y = BCS, fill = `11q`)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Hailstorm 11q", y = "BCS") +
  scale_fill_npg() +
  # coord_cartesian(ylim = c(1, 1000), expand = TRUE) +
  annotate("text",
    x = 1.5, y = max(Table_merged_ed$BCS) * 1.05, # Adicionar valor-p ajustado manualmente
    label = paste("p.adj =", round(adjusted_p_value_bcs_11q, 4)), color = "black", size = 4
  ) +
  theme_classic()

p5
ggsave(
  filename = here("results", "BCS_hailstorm_11q_plot.pdf"),
  plot = p5,
  device = "pdf",
  width = 5.52,
  height = 5.7
)
# Encontrar o índice da região "22q" na tabela de resultados
indice_22q <- which(results_table$Regiao == "22q")

# Acessar os valores-p ajustados para "22q"
adjusted_p_value_bcs_22q <- results_table$Adjusted_P_Value_BCS[indice_22q]

# Gráfico para FCS
p6 <- ggplot(Table_merged_ed, aes(x = `22q`, y = BCS, fill = `22q`)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Hailstorm 22q", y = "BCS") +
  scale_fill_npg() +
  # coord_cartesian(ylim = c(1, 1000), expand = TRUE) +
  annotate("text",
    x = 1.5, y = max(Table_merged_ed$BCS) * 1.05, # Adicionar valor-p ajustado manualmente
    label = paste("p.adj =", round(adjusted_p_value_bcs_22q, 4)), color = "black", size = 4
  ) +
  theme_classic()
p6


p6
ggsave(
  filename = here("results", "BCS_hailstorm_22q_plot.pdf"),
  plot = p6,
  device = "pdf",
  width = 5.52,
  height = 5.7
)
# Encontrar o índice da região "5p" na tabela de resultados
indice_5p <- which(results_table$Regiao == "5p")

# Acessar os valores-p ajustados para "5p"
adjusted_p_value_bcs_5p <- results_table$Adjusted_P_Value_BCS[indice_5p]


p7 <- ggplot(Table_merged_ed, aes(x = `5p`, y = BCS, fill = `5p`)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Hailstorm 5p", y = "BCS") +
  scale_fill_npg() +
  # coord_cartesian(ylim = c(1, 1000), expand = TRUE) +
  annotate("text",
    x = 1.5, y = max(Table_merged_ed$BCS) * 1.05, # Add pvalue manual
    label = paste("p.adj =", round(adjusted_p_value_bcs_5p, 4)), color = "black", size = 4
  ) +
  theme_classic()
p7
ggsave(
  filename = here("results", "BCS_hailstorm_5p_plot.pdf"),
  plot = p7,
  device = "pdf",
  width = 5.52,
  height = 5.7
)

combined_plot <- grid.arrange(p5, p6, p7, ncol = 3)



# Saving all packages used to create these plots
writeLines(
  capture.output(sessionInfo()),
  here("package_version", "versions_02_FigureS3b_c.txt")
)
