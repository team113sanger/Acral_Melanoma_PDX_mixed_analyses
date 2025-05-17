title: "Mutational and CNA Signatures plots "
author: "Annie"
date: "2025-05-16"

##Loading packages

### diretorio da nova lib 24-06-2024
### ou diretório acessível em todos os nós
mylib <- "/data04/projects04/PatriciaPossik/pdx_la_sanger/lib/Rpackages/4.4.1"
libdir <- "/data04/tools/R/Rpackages/4.4.1"

### comando para efetivar a seleção de libPaths
.libPaths(c(mylib,libdir))
.libPaths()

### set CRAN Mirror https://cran-r.c3sl.ufpr.br/
options(repos = c("CRAN" = "https://cran-r.c3sl.ufpr.br/"))
 
### set CRAN Mirror https://cran-r.c3sl.ufpr.br/  
options(repos = c("CRAN" = "https://cran-r.c3sl.ufpr.br/"))  

##Dependencies for this project are recorded and managed with renv
## To setup your Rstudio environment you might need to run:
renv::restore()

library(dplyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(data.table)
library(here)
library(renv)
library(ggsci)
library(colorspace)

install.packages("devtools")
install.packages("reticulate")

library(reticulate)
library(devtools)
#if (!require("devtools")) install.packages("devtools")
#devtools::install_github("AlexandrovLab/SigProfilerAssignmentR")

library(SigProfilerAssignmentR)

# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/05_Figure_S4c.R")

# First, to CNA signatures

# Load CNA signatures data from CSV file
# Note: Consider using here::here() for more robust path handling
acts_ed <- read.csv("./data/COSMIC_CNV48_activities_fev_25.csv", header = TRUE)

# Convert 'Samples' column to factor with specific order
# This ensures samples appear in consistent order in plots
sample_order <- c("PD53355c", "PD53333f", "PD53330c", "PD53337a", "PD53363a", 
                  "PD53346a", "PD53369a", "PD53352d", "PD53343h", "PD53338a", 
                  "PD53368a", "PD53351a", "PD53350d", "PD53347c", "PD53348a", 
                  "PD53359c", "PD53357c", "PD53362a", "PD53367a", "PD53345a", 
                  "PD53364c", "PD53349d", "PD53361a", "PD53354a", "PD53334a", 
                  "PD53332d", "PD53358a", "PD53366a", "PD53331a", "PD53339a", 
                  "PD53341a", "PD53342a", "PD53353a", "PD53365a")

acts_ed$Samples <- factor(acts_ed$Samples, levels = sample_order)

# Remove any rows with NA values
acts_ed <- na.omit(acts_ed)

# Verify the factor levels are correct
levels(acts_ed$Samples)

# Reshape data from wide to long format for plotting
# This converts signature columns into rows with Signature/Mutations columns
acts_tidy = acts_ed %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations')

# Create custom color palette combining NPG colors and additional colors
npg_colors <- pal_npg("nrc")(10)  # Get 10 colors from NPG palette
extra_colors <- c("#FF7F00", "#6A3D9A", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
my_palette <- c(npg_colors, extra_colors)  # Combine palettes (total 16 colors)

# Create stacked barplot showing absolute contribution of CNA signatures
plot1 <- ggplot(acts_tidy) +
  aes(x = Samples, y = Mutations, fill = Signature) +
  geom_bar(position='stack', stat = 'identity') +  # Stacked bars
  theme_classic() +  # Clean theme
  labs(x = "Samples", y = "CNA Signatures Contribution (Absolute)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels
  scale_fill_manual(values = my_palette)  # Apply custom colors

# Adjust plot limits and save
plot1 <- plot1 + coord_cartesian(xlim = c(0,35), expand = FALSE)

# Save plot as PDF
ggsave(
  filename = here("results", "CNA_signature_absolute.pdf"),
  plot = plot1,
  device = "pdf",
  width = 9.16,
  height = 4.44
)

# Add light gray color for any "empty" category
signature_levels <- unique(acts_tidy$Signature)
custom_palette <- setNames(my_palette[1:length(signature_levels)], signature_levels)
custom_palette["empty"] <- "#D3D3D3"  # Light gray for missing data

# Create relative contribution plot (percentage stacked)
plot2 <- ggplot(acts_tidy, aes(x = Samples, y = Mutations, fill = Signature)) +
  geom_bar(position = 'fill', stat = 'identity') +  # Percentage stacked
  theme_classic() +
  labs(x = "Samples", y = "CNA Signatures Contribution (Relative)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_palette) 

# Adjust plot limits and save
plot2 <- plot2 + coord_cartesian(xlim = c(0, 35), expand = FALSE)

ggsave(
  filename = here("results", "CNA_signature_relative.pdf"),
  plot = plot2,
  device = "pdf",
  width = 9.16,
  height = 4.44
)

# Second to mutational signatures

# Load mutational signatures data
d  <- read.table("./data/COSMIC_SBS96_Activities.txt", header=TRUE)

# Load metadata (though not used in this script)
Table_metadata <- read.csv("./data/Table_metadata_REVIEWED_18_01_25.csv")

mut_sig <- d  # Create working copy

# Apply same sample ordering as before
mut_sig$Samples <- factor(mut_sig$Samples, levels = sample_order)

# Remove NA values
mut_sig <- na.omit(mut_sig)

# Verify factor levels
levels(mut_sig$Samples)

# Reshape mutational signatures data to long format
acts_tidy_mut = mut_sig %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations')

# Create lighter version of NPG palette for relative plot
npg_colors <- pal_npg("nrc")(10)
lighter_npg_colors <- lighten(npg_colors, amount = 0.35)  # Lighten by 35%

# Create palette with light gray for "empty" category
signature_levels <- unique(acts_tidy_mut$Signature)
custom_palette <- setNames(my_palette[1:length(signature_levels)], signature_levels)
custom_palette["empty"] <- "#D3D3D3"

# Create absolute contribution plot for mutational signatures
plot3 <- ggplot(acts_tidy_mut) +
  aes(x = Samples, y = Mutations, fill = Signature) +
  geom_bar(position='stack', stat = 'identity') +
  theme_classic() +
  labs(x = "Samples", y = " Mutational Signatures Contribution")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = npg_colors)  # Use original colors

# Adjust limits and save
plot3 <- plot3 + coord_cartesian(xlim = c(0,35), ylim = c(0,200), expand = FALSE)

ggsave(
  filename = here("results", "Mutational_signatures_absolute.pdf"),
  plot = plot3,
  device = "pdf",
  width = 9.16,
  height = 4.44
)

# Create relative contribution plot using lighter colors
plot4 <- ggplot(acts_tidy_mut) +
  aes(x = Samples, y = Mutations, fill = Signature) +
  geom_bar(position='fill', stat = 'identity') +  # Percentage stacked
  theme_classic() +
  labs(x = "Samples", y = " Relative Mutational Signatures Contribution")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = lighter_npg_colors)  # Use lighter colors

# Adjust limits and save
plot4 <- plot4 + coord_cartesian(xlim = c(0,35), expand = FALSE)

ggsave(
  filename = here("results", "Mutational_signatures_relative.pdf"),
  plot = plot4,
  device = "pdf",
  width = 9.16,
  height = 4.44
)

# Save session info for reproducibility
writeLines(
  capture.output(sessionInfo()),
  here("package_version", "versions_script_05_Figure_S4c.txt")
)