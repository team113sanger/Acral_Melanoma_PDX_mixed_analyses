### Author: Annie Squiavinato  and Antonio
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
renv::init() # To initialise revn
#renv::restore()
#renv::init() # To initialise revn
#renv::install("glue") # To add a package from CRAN
#renv::install("BiocManager") # To add a package from Bioconductor
#renv::install("bioC::ggplot2") # To add a package from Bioconductor
#renv::snapshot() # To record any packages you've used and their versions and sources into the renv.lock file
#renv::restore() # To rebuild an environment from the renv.lockfile

# Libraries
library(ComplexHeatmap)
library(here)
library(ggplot2)

renv::snapshot()
# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/01_Figure_1b.R")

# Data
clin_mat <- read.table("./data/figure_1_table_23_04_25.tsv", header = TRUE, sep = "\t", row.names = 1)
final_mat <- as.data.frame(t(clin_mat))
rows <- rownames(final_mat)

anno_mat <- read.table("./data/figure_1b_table_23_04_25.tsv", header = TRUE, sep = "\t", row.names = 1)

# Colours
col <- c(Yes = "#00a087")
anato_col <- c(Plantar = "#e64b35", Lymph_Node = "#f1998e", Subungual = "#f9d6d2", In_transit_Metastasis = "#9f2414")
type_col <- c(Primary = "#4dbbd5", Regional_Metastasis = "#278ea5", Recurrence = "#98d8e7")

# Function to draw squares (white background)
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "white")) # white background
  },
  Yes = function(x, y, w, h) {
    grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col["Yes"], col = NA))
  }
)

# Notes (each square visible)
ha <- HeatmapAnnotation(
  df = anno_mat,
  which = "col",
  col = list(
    Anatomic_site = anato_col,
    Sample_type = type_col
  ),
  simple_anno_size = unit(0.6, "cm"),
  annotation_name_side = "left",
  gp = gpar(col = "white") # White borders to avoid overlapping
)

# Heatmap 
p1 <- oncoPrint(
  final_mat,
  alter_fun = alter_fun,
  col = col,
  right_annotation = NULL,
  top_annotation = ha,
  show_pct = FALSE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_order = rows,
  alter_fun_is_vectorized = FALSE, 
  heatmap_legend_param = list(title = "Data and Models")
)

#to save
pdf(file = here("results", "Figure_1B.pdf"), width = 12.38, height = 4.13)
draw(p1)  
dev.off()

# saving all packages used to create these plots
writeLines(
  capture.output(sessionInfo()),
  here("package_version", "versions_script_01_figure_01_plot.txt")
)
