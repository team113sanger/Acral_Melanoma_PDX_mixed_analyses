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

## Loading packages
library(dplyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(gridExtra)
library(data.table)
library(ggpubr)
library(ggsci)
library(dunn.test)
library(here)
renv::snapshot()

# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/04_Figure_4a.R")
################### 1) plot comparing the number of SNP_DEL_INS_PDX vs tumor patients

variant_number_table <- read.csv("./data/variant_counts_cna_score_table.csv")
variant_number_table <- as.data.frame(variant_number_table)
write.csv(
  variant_number_table,
  here("results", "Metadata_PDX_tumor_mutations_fev_2025.csv"),
  row.names = TRUE,
  sep = ","
)

# testing normality
shapiro_result <- shapiro.test(variant_number_table$Number_SNP_DEL_INS_tumor)
shapiro_result <- shapiro.test(variant_number_table$Number_SNP_DEL_INS_PDX)

# checking the colnames
glimpse(variant_number_table)

# creating colors personalized
npg_colors <- pal_npg("nrc")(10) # Cores NPG
extra_colors <- c("#FF7F00", "#6A3D9A", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
my_palette <- c(npg_colors, extra_colors) # Total de 16 cores

# correlation using peason since data has normal distribution
cor_result <- cor.test(variant_number_table$Number_SNP_DEL_INS_tumor, variant_number_table$Number_SNP_DEL_INS_PDX, method = "pearson")

# plotting
plot1 <- ggplot(variant_number_table, aes(x = Number_SNP_DEL_INS_tumor, y = Number_SNP_DEL_INS_PDX, color = Patients)) + # PDX: triangle (shape = 17)
  geom_point(aes(shape = "PDX"), alpha = 0.5, size = 4, position = position_jitter(width = 0.2, height = 0.2)) + # Tumor: circle (shape = 16)
  geom_point(aes(x = Number_SNP_DEL_INS_tumor, y = Number_SNP_DEL_INS_tumor, shape = "Tumor"),
    alpha = 0.5, size = 4, position = position_jitter(width = 0.2, height = 0.2)
  ) +
  labs(x = "Number of SNP/DEL/INS Tumor", y = "Number of SNP/DEL/INS PDX") +
  theme_classic() +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +
  # geom_text(x = Inf, y = -Inf, inherit.aes = FALSE,
  # label = paste("R²:", round(cor_result$estimate, 2), "\n",
  # "p-value:", format(cor_result$p.value, scientific = FALSE, 3), sep = ""),
  # hjust = 1.1, vjust = -1, size = 4) +
  scale_color_manual(values = my_palette) + # keep the same color for each patient
  scale_shape_manual(name = "Sample Type", values = c("PDX" = 17, "Tumor" = 16)) + # Different symbols to PDX and Tumor samples
  coord_cartesian(ylim = c(0, 150), xlim = c(3, 150), expand = TRUE)

plot1

# saving
ggsave(
  filename = here("results", "Correlation_number_of_SNP_DEL_INS_pdx_vs_tumor.pdf"),
  plot = plot1,
  device = "pdf",
  width = 5.52,
  height = 6.7
)
################### 2) Plot comparing the FCS and BCS score between PDX-X1 and sample tumors

variant_number_table <- read.delim("./data/metadata_pdx_tumor_compartions_fev_2025.csv", header = TRUE, stringsAsFactors = FALSE)

# testing normality, To say that is a normal distributiom the p value needs to be great than 0,05.
shapiro_result <- shapiro.test(variant_number_table$FCS_pdx)
shapiro_result <- shapiro.test(variant_number_table$FCS_tum)
shapiro_result <- shapiro.test(variant_number_table$BCS_tum)
shapiro_result <- shapiro.test(variant_number_table$BCS_pdx)

# selecting only the columns that I will need
CNApp_FCS <- variant_number_table %>% dplyr::select(Patients, FCS_tum, FCS_pdx)

CNApp_BCS <- variant_number_table %>% dplyr::select(Patients, BCS_tum, BCS_pdx)

# CNApp_GCS <- variant_number_table %>% select(Patients, GCS_tum,  GCS_pdx)

# To calculate the correlation between BCS_tum and BCS_pdx
cor_result_BCS <- cor.test(variant_number_table$BCS_tum, variant_number_table$BCS_pdx, method = "spearman")

# Checking if the object was created normally
print(cor_result_BCS)

# plotting, using the same colors from previus plot
npg_colors <- pal_npg("nrc")(10) # Cores NPG
extra_colors <- c("#FF7F00", "#6A3D9A", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
my_palette <- c(npg_colors, extra_colors) # Total de 16 cores

plot2 <- ggplot(variant_number_table, aes(x = BCS_tum, y = BCS_pdx, color = Patients)) +
  geom_point(aes(shape = "PDX"), alpha = 0.5, size = 4, position = position_jitter(width = 0.2, height = 0.2)) +
  geom_point(aes(x = BCS_tum, y = BCS_tum, shape = "Tumor"), alpha = 0.5, size = 4.4, position = position_jitter(width = 0.2, height = 0.2)) +
  labs(x = "BCS Score Tumor", y = "BCS Score PDX") +
  theme_classic() +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +
  # geom_text(x = Inf, y = -Inf, inherit.aes = FALSE,
  # label = paste("R²:", round(cor_result_BCS$estimate, 2), "\n",
  # "p-value:", format(cor_result_BCS$p.value, scientific = FALSE, 3), sep = ""),
  # hjust = 1.1, vjust = -1, size = 4) +
  scale_color_manual(values = my_palette) + # Mantém a mesma cor para cada paciente
  scale_shape_manual(name = "Sample Type", values = c("PDX" = 17, "Tumor" = 16)) + # Different symbols to PDX and Tumor samples
  coord_cartesian(ylim = c(0, 100), expand = TRUE)

plot2
# saving
ggsave(
  filename = here("results", "Correlation_BCS_pdx_vs_tumor.pdf"),
  plot = plot2,
  device = "pdf",
  width = 5.52,
  height = 5.5
)

# To calculate the correlation between FCS_tum and FCS_pdx
cor_result_FCS <- cor.test(variant_number_table$FCS_tum, variant_number_table$FCS_pdx, method = "spearman")

# Checking if the object was created normally
print(cor_result_FCS)

# plotting
plot3 <- ggplot(variant_number_table, aes(x = FCS_tum, y = FCS_pdx, color = Patients)) +
  geom_point(aes(shape = "PDX"), alpha = 0.5, size = 4, position = position_jitter(width = 0.2, height = 0.2)) +
  geom_point(aes(x = FCS_tum, y = FCS_tum, shape = "Tumor"), alpha = 0.5, size = 4.4, position = position_jitter(width = 0.2, height = 0.2)) +
  labs(x = "FCS Score Tumor", y = "FCS Score PDX") +
  theme_classic() +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +
  # geom_text(x = Inf, y = -Inf, inherit.aes = FALSE,
  # label = paste("R²:", round(cor_result_FCS$estimate, 2), "\n",
  # "p-value:", format(cor_result_FCS$p.value, scientific = FALSE, 3), sep = ""),
  #  hjust = 1.1, vjust = -1, size = 4) +
  scale_color_manual(values = my_palette) + # Mantém a mesma cor para cada paciente
  scale_shape_manual(name = "Sample Type", values = c("PDX" = 17, "Tumor" = 16)) + # Diferentes símbolos para PDX e Tumor
  coord_cartesian(ylim = c(0, 1250), expand = TRUE)

plot3

# saving
library(here)
ggsave(
  filename = here("results", "Correlation_FCS_score_tumor_pdx.pdf"),
  plot = plot3,
  device = "pdf",
  width = 5.52,
  height = 5.5
)

################### 3) Comparing Mutational status vs FCS and BCS score

# Open the table that I previus saved
Metadata <- read.csv("./data/Table_metadata_REVIEWED_18_01_25.csv")

# List of samples to be excluded, cellularity 0.35>x (QC done in Jan/Feb 2025, Sequenza data)
# We used only samples with celularity higher than 0.35m one sample per patient.

# List of samples to be excluded, cellularity 0.35>x (QC done in Jan/Feb 2025, Sequenza data)
samples <- c(
  "PD53330c", "PD53332d", "PD53333f", "PD53334a", "PD53338a", "PD53342a", "PD53343h", "PD53345a", "PD53347c",
  "PD53348a", "PD53349d", "PD53350d", "PD53351a", "PD53352d", "PD53354a", "PD53355c", "PD53357c", "PD53358a",
  "PD53359c", "PD53361a", "PD53362a", "PD53364c", "PD53365a", "PD53368a", "PD53369a"
)

# Filtering data
Table_QC_fev_2025 <- Metadata %>% filter(Tumor_Sample_Barcode %in% samples)


show_col(pal_npg("nrc")(10))

# Changing the data to long format
df_long <- Table_QC_fev_2025 %>%
  pivot_longer(cols = c(Quadruple_negative, RAS_mut, KIT_mut, NF1_mut), names_to = "Mutation_status", values_to = "Status")

# filtering only samples that have yes status
df_filtered <- df_long %>%
  filter(Status == "Yes")

# Reorganazing the level of variable Mutation_status
df_filtered$Mutation_status <- factor(df_filtered$Mutation_status,
  levels = c("RAS_mut", "KIT_mut", "NF1_mut", "Quadruple_negative")
)

######### plot BCS
library(dunn.test)

# plot
plot4 <- ggplot(df_filtered, aes(x = Mutation_status, y = BCS, fill = Mutation_status)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Status Mutacional", y = "BCS") +
  scale_fill_npg() + # Aplicando a paleta NPG
  stat_compare_means(aes(group = Mutation_status), method = "kruskal.test", label.y = max(df_filtered$BCS) * 1.1) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 100), expand = TRUE)
plot4

# saving
ggsave(
  filename = here("results", "BCS_status_mutational.pdf"),
  plot = plot4,
  device = "pdf",
  width = 6.34,
  height = 4.34
)

# doing pos tests: comparing every mutation status vs 4WT, just to check. Kruskal Wallis did not show any significance.
dunn_result <- dunn.test(df_filtered$BCS, df_filtered$Mutation_status, method = "bonferroni")
print(dunn_result)
glimpse(df_filtered)

write.csv(
  dunn_result,
  here("results", "BCS_mutational_status test_dunn_bonferroni"),
  row.names = TRUE,
  sep = ","
)

######### plot FCS

# plot
plot5 <- ggplot(df_filtered, aes(x = Mutation_status, y = FCS, fill = Mutation_status)) +
  geom_boxplot(
    width = 0.6, outlier.shape = 1, outlier.size = 1.6, outlier.alpha = 0,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(alpha = 0.5, color = "black", width = 0.2) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 4, fill = "gray") +
  labs(x = "Status Mutacional", y = "FCS") +
  scale_fill_npg() + # Aplicando a paleta NPG
  stat_compare_means(aes(group = Mutation_status), method = "kruskal.test", label.y = max(df_filtered$FCS) * 0.9) +
  theme_classic() +
  coord_cartesian(ylim = c(1, 1000), expand = TRUE)
plot5
# saving
ggsave(
  filename = here("results", "FCS_status_mutational.pdf"),
  plot = plot5,
  device = "pdf",
  width = 6.34,
  height = 4.34
)

# doing pos tests: comparing every mutation status vs 4WT, just to check. Kruskal Wallis did not show any significance.
dunn_result2 <- dunn.test(df_filtered$FCS, df_filtered$Mutation_status, method = "bonferroni")
print(dunn_result)

write.csv(
  dunn_result2,
  here("results", "FCS_mutational_status test_dunn_bonferroni"),
  row.names = TRUE,
  sep = ","
)

# saving all packages used to create these plots
writeLines(
  capture.output(sessionInfo()),
  here("package_version", "versions_script_04_Figure_4a.txt")
)
