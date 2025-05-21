### Author:Flavia Aguiar, Danielle Carvalho and Annie Squiavinato
### date: "2025-04-25"  

#loading library
library("ComplexHeatmap")
library("circlize")
library("cluster")
library("here")
library("dplyr")

# Load required packages
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Set working directory using here package
here::i_am("heatmap_drugscreen/06_Figure_5a.R")

# Read IC50 data file
data_IC50_edited <- read.csv("./IC50_acral_drugs_annotations.csv", 
                             header = TRUE, 
                             stringsAsFactors = FALSE, 
                             check.names = FALSE, 
                             row.names = 1)

# Convert to data frame
data_IC50_edited <- as.data.frame(data_IC50_edited)

# Exclude first column as it's not required for analysis
data_IC50_edited <- data_IC50_edited[,-1]

# Convert comma decimal separators to points and transform to numeric
data_IC50_edited <- data_IC50_edited %>%
  mutate(across(where(is.character), ~ as.numeric(gsub(",", ".", .)))
         
# Store original row and column names
col_names <- colnames(data_IC50_edited)
row_names <- rownames(data_IC50_edited)
         
# Prepare data for log transformation
 colnames(data_IC50_edited) <- NULL
 rownames(data_IC50_edited) <- NULL
 colnames(data_IC50_edited) <- seq_along(colnames(data_IC50_edited))
         
# Replace NAs with 1000 for log transformation
data_IC50_edited[is.na(data_IC50_edited)] <- 1000
         
# Apply log10 transformation
data_IC50_edited_log10 <- log10(data_IC50_edited)                     
         
# Restore NAs (values that were 1000 before transformation)
data_IC50_edited_log10[data_IC50_edited_log10 == 3.0] <- NA
         
# Restore original names
colnames(data_IC50_edited_log10) <- col_names 
rownames(data_IC50_edited_log10) <- row_names
         
# Define desired row order for heatmap
 row_order <- c("Cobimetinib", "Trametinib", "Ripretinib", "Avapritinib", "XAV939", 
                        "GSK461364", "Volasertib", "MLN8054", "GSK1070916", "BAY-320", 
                        "Abemaciclib", "Palbociclib", "Dinaciclib", "Nutlin-3a", "KU-60019", 
                        "Berzosertib", "Ceralasertib", "CHIR-124", "SCH900776", "Adavosertib", 
                        "Olaparib", "Talazoparib", "BIIB021", "ML-792", "Entinostat", 
                        "Vorinostat", "GSK126", "JQ1", "Molibresib")
         
# Reorder rows
data_IC50_edited_log10_ordered <- data_IC50_edited_log10[row_order, ]
         
# Create matrix with NA values for heatmap
mat_with_na <- as.matrix(data_IC50_edited_log10_ordered)
         
# Define color mapping for heatmap
col_fun <- colorRamp2(c(-3, -2, -1, 1, 2, 3), 
                               c("darkred", "firebrick3", "lightcoral", 
                                 "steelblue1", "steelblue", "royalblue4"), 
                               space = "LAB")
         
# Create heatmap
 Heatmap(mat_with_na, 
         name = "IC50 (log10)", 
         na_col = "white", 
         rect_gp = gpar(col = "gray", lwd = 0.3),  
         col = col_fun, 
         cluster_rows = FALSE, 
         cluster_columns = FALSE,  
          heatmap_legend_param = list(at = seq(-3, 3, length.out = 9)))
          
 
# Save heatmap to PDF
pdf("heatmap2.pdf", width = 8, height = 6)
plot(Heatmap(mat_with_na, 
name = "IC50 (log10)", 
na_col = "white", 
rect_gp = gpar(col = "gray", lwd = 0.3),  
col = col_fun, 
cluster_rows = FALSE, 
cluster_columns = FALSE,  
heatmap_legend_param = list(at = seq(-3, 3, length.out = 9))))
dev.off()
