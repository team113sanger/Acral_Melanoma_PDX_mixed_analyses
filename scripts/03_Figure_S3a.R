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

# loading packages  
library(tidyverse)  
library(fs)  
library(GenomicRanges)  
library(biomaRt)  
library(tibble)  
library(ggpubr)  
library(ggsci)  
library(ggplot2)  
library(dplyr)  
library(gridExtra)  
library(here)  
renv::snapshot()

# i_am is a function from the here package that helps to set the working directory relative to th execution script
here::i_am("scripts/03_Figure_S3a.R")
# Choosing a list of genes ----  
# Driver genes for the figure  
genes <- read.csv("./data/genes_to_get_copy_number_hails.csv", header = TRUE, sep = "\t")  

# Gene positions from biomart  
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://feb2021.archive.ensembl.org")  
gene_dat <- getBM(  
  attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),  
  filters = c("hgnc_symbol", "chromosome_name"),  
  values = list("hgnc_symbol" = genes$gene_name, "chromosome_name" = c(1:22, "X", "Y")),  
  mart = mart  
)  
gene_dat$chromosome_name <- gsub("^", "chr", gene_dat$chromosome_name)  
genes_gr <- as(dplyr::rename(gene_dat, "start" = start_position, "end" = end_position), "GRanges")  

# copy number segments  
segment_files <- dir_ls("./data/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/",  
                        glob = "*segments.txt",  
                        recurse = TRUE  
)  
seg_dat <- read_tsv(segment_files, id = "Sample") |>  
  mutate(Sample = gsub("\\_.*", "", basename(Sample))) |>  
  dplyr::rename("start" = start.pos, "end" = end.pos)  

# Sample ploidy data  
ploify_files <- dir_ls("./data/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/",  
                       glob = "*confints_CP.txt",  
                       recurse = TRUE  
)  
ploidy_dat <- read_tsv(ploify_files, id = "Sample") |>  
  mutate(Sample = gsub("\\_.*", "", basename(Sample))) |>  
  group_by(Sample) |>  
  distinct(round(ploidy.mean.cn)) |>  
  dplyr::rename("Ploidy" = "round(ploidy.mean.cn)")  

# Manually setting the ploidy for samples where ploidy_mean doesnt reflect well  
ploidy_dat[34, 2] <- 4  
ploidy_dat[49, 2] <- 3  

# Add ploidy to segments  
final_seg <- left_join(seg_dat, ploidy_dat, by = "Sample")  
seg_gr <- as(final_seg, "GRanges")  

# Searching for overlaps between genes from list and segments  
overlaps <- findOverlaps(seg_gr, genes_gr)  
overlap_list <- list()  

for (i in seq_along(overlaps)) {  
  seg_id <- queryHits(overlaps)[i]  
  gene_id <- subjectHits(overlaps)[i]  
  
  overlap_data <- data.frame(  
    Sample = seg_gr[seg_id]$Sample,  
    Gene = genes_gr[gene_id]$hgnc_symbol,  
    Copy_number = seg_gr[seg_id]$CNt,  
    Ploidy = seg_gr[seg_id]$Ploidy,  
    seg_start = start(seg_gr[seg_id]),  
    seg_end = end(seg_gr[seg_id]),  
    gene_start = start(genes_gr[gene_id]),  
    gene_end = end(genes_gr[gene_id])  
  )  
  
  overlap_list[[i]] <- overlap_data  
}  

copy_number_info <- do.call(rbind, overlap_list)  

processed_cn <- copy_number_info |>  
  mutate(  
    overlap_start = pmax(seg_start, gene_start),  
    overlap_end = pmin(seg_end, gene_end),  
    overlap_length = pmax(0, overlap_end - overlap_start), # Ensure no negative values  
    gene_length = gene_end - gene_start,  
    proportion = overlap_length / gene_length  
  )  

cn_filtered <- processed_cn |>  
  group_by(Sample, Gene) |>  
  slice_max(order_by = proportion, n = 1, with_ties = FALSE) |>  
  ungroup()  

classified <- cn_filtered |>  
  mutate(Alteration = case_when(  
    Copy_number == Ploidy ~ "Neutral",  
    Copy_number < Ploidy & Copy_number != 0 ~ "Loss",  
    Copy_number == 0 ~ "Del",  
    Copy_number > Ploidy & Copy_number < Ploidy * 2 ~ "Gain",  
    Copy_number >= Ploidy * 2 ~ "Amp"  
  ))  

cna_matched_list <- classified |>  
  dplyr::select(Sample, Gene, Alteration, Copy_number)  

# List of samples to exclude, cellularity 0.35>x (QC done in Jan/Feb 2025, Sequenza data)  
ids <- c(  
  "PD53330c", "PD53332d", "PD53333f", "PD53334a", "PD53337a", "PD53338a", "PD53342a", "PD53343h", "PD53345a", "PD53347c",  
  "PD53348a", "PD53349d", "PD53350d", "PD53351a", "PD53352d", "PD53354a", "PD53355c", "PD53357c", "PD53358a",  
  "PD53359c", "PD53361a", "PD53362a", "PD53364c", "PD53365a", "PD53368a", "PD53369a"  
)  

# Filtering the data  
Table_hail_expression <- cna_matched_list %>% filter(Sample %in% ids)  

write_tsv(Table_hail_expression, "Genes_ploydia_cn.tsv")  
# table describing samples that have or not hailstorms  
Table_hailstorms <- read.csv("./data/allSamples_hailstorm_regions_corrected.csv", header = TRUE, sep = "\t")  

# split the chr column  
hail_presence <- split(Table_hailstorms, Table_hailstorms$Arm)  
names(hail_presence)  

# chr11  
dados_11q <- hail_presence[["11q"]]  
dados_11q <- setNames(data.frame(t(dados_11q[, -1])), dados_11q[, 1]) # transposing  
dados_11q <- rownames_to_column(dados_11q, var = "Sample") # placing IDs in the first column  
t_dados_11q <- dados_11q %>% filter(Sample %in% ids)  

# chr12  
dados_12q <- hail_presence[["12q"]]  
dados_12q <- setNames(data.frame(t(dados_12q[, -1])), dados_12q[, 1]) # transposing  
dados_12q <- rownames_to_column(dados_12q, var = "Sample") # placing IDs in the first column  
t_dados_12q <- dados_12q %>% filter(Sample %in% ids)  

# chr12  
dados_22q <- hail_presence[["22q"]]  
dados_22q <- setNames(data.frame(t(dados_22q[, -1])), dados_22q[, 1]) # transposing  
dados_22q <- rownames_to_column(dados_22q, var = "Sample") # placing IDs in the first column  
t_dados_22q <- dados_22q %>% filter(Sample %in% ids)  

# chr5 
dados_5p <- hail_presence[["5p"]]  
dados_5p <- setNames(data.frame(t(dados_5p[, -1])), dados_5p[, 1]) # transposing  
dados_5p <- rownames_to_column(dados_5p, var = "Sample") # placing IDs in the first column  
t_dados_5p <- dados_5p %>% filter(Sample %in% ids)  

# Split the data frame by unique values of chr_arm  
df_split <- split(Table_hail_expression, Table_hail_expression$Gene)  

# View results  
names(df_split) # Shows the names of the list elements (unique values of chr_arm)  
head(df_split[[1]]) # Shows the first subset (data from the first chr_arm)  

# Access data  
dados_CCND1 <- df_split[["CCND1"]]  
dados_CDK4 <- df_split[["CDK4"]]  
dados_CRKL <- df_split[["CRKL"]]  
dados_EP300 <- df_split[["EP300"]]  
dados_GAB2 <- df_split[["GAB2"]]  
dados_LZTR1 <- df_split[["LZTR1"]]  
dados_MDM2 <- df_split[["MDM2"]]  
dados_SOX10 <- df_split[["SOX10"]]  
dados_TERT <- df_split[["TERT"]]  
dados_PAK1 <- df_split[["PAK1"]]  
dados_PRAME <- df_split[["PRAME"]]  

# merging the CCND1 and CRKL and other genes according to chr correspondent  
dados_CCND1_11q <- merge(t_dados_11q, dados_CCND1, by = "Sample")  
dados_CDK4_12q <- merge(t_dados_12q, dados_CDK4, by = "Sample")  
dados_CRKL_22q <- merge(t_dados_22q, dados_CRKL, by = "Sample")  
dados_EP300_22q <- merge(t_dados_22q, dados_EP300, by = "Sample")  
dados_GAB2_11q <- merge(t_dados_11q, dados_GAB2, by = "Sample")  
dados_LZTR1_22q <- merge(t_dados_22q, dados_LZTR1, by = "Sample")  
dados_MDM2_12q <- merge(t_dados_12q, dados_MDM2, by = "Sample")  
dados_SOX10_22q <- merge(t_dados_22q, dados_SOX10, by = "Sample")  
dados_TERT5_5p <- merge(t_dados_5p, dados_TERT, by = "Sample")  
dados_PAK1_11q <- merge(t_dados_11q, dados_PAK1, by = "Sample")  
dados_PRAME_22q <- merge(t_dados_22q, dados_PRAME, by = "Sample")  

# changing the ID to correspond to id_RNA  

# Using gsub to replace all occurrences of "D" with "R" in the column  
dados_CCND1_11q$Sample <- gsub("D", "R", dados_CCND1_11q$Sample)  
dados_CDK4_12q$Sample <- gsub("D", "R", dados_CDK4_12q$Sample)  
dados_CRKL_22q$Sample <- gsub("D", "R", dados_CRKL_22q$Sample)  
dados_EP300_22q$Sample <- gsub("D", "R", dados_EP300_22q$Sample)  
dados_GAB2_11q$Sample <- gsub("D", "R", dados_GAB2_11q$Sample)  
dados_LZTR1_22q$Sample <- gsub("D", "R", dados_LZTR1_22q$Sample)  
dados_MDM2_12q$Sample <- gsub("D", "R", dados_MDM2_12q$Sample)  
dados_SOX10_22q$Sample <- gsub("D", "R", dados_SOX10_22q$Sample)  
dados_TERT5_5p$Sample <- gsub("D", "R", dados_TERT5_5p$Sample)  
dados_PAK1_11q$Sample <- gsub("D", "R", dados_PAK1_11q$Sample)  
dados_PRAME_22q$Sample <- gsub("D", "R", dados_PRAME_22q$Sample)  

TPM_RNAseq <- read.table("./data/TPMS_HTSeq_STAR_ENSv103_ERCC_TPM_M_ENS_gene_name_IDS_QC_PASSED_xenofiltered.txt", header = TRUE)  

genes_exp <- TPM_RNAseq[TPM_RNAseq$external_gene_name %in% c("GAB2", "CCND1", "CDK4", "EP300", "GAB2", "LZTR1", "MDM2", "SOX10", "TERT", "CRKL", "PAK1", "PRAME"), ]  
genes_exp <- as.data.frame(genes_exp)  

genes_exp <- genes_exp[, -2]  
# transposing to each gene get in a individual column  
genes_exp_t <- t(genes_exp)  

# putting the first line as colname  
colnames(genes_exp_t) <- genes_exp_t[1, ]  

# excluding the first line (repeated)  
genes_exp_t <- genes_exp_t[-1, ]  

# rownames to column to do a merge  
genes_exp_t <- as.data.frame(genes_exp_t)  
genes_exp_t <- rownames_to_column(genes_exp_t, var = "Sample")  

# merging with copy number data for CCND1 gene  
CCND1_11q_merge <- merge(dados_CCND1_11q, genes_exp_t, by = "Sample")  
CDK4_12q_merge <- merge(dados_CDK4_12q, genes_exp_t, by = "Sample")  
CRKL_22q_merge <- merge(dados_CRKL_22q, genes_exp_t, by = "Sample")  
EP300_22q_merge <- merge(dados_EP300_22q, genes_exp_t, by = "Sample")  
GAB2_11q_merge <- merge(dados_GAB2_11q, genes_exp_t, by = "Sample")  
LZTR1_22q_merge <- merge(dados_LZTR1_22q, genes_exp_t, by = "Sample")  
MDM2_12q_merge <- merge(dados_MDM2_12q, genes_exp_t, by = "Sample")  
SOX10_22q_merge <- merge(dados_SOX10_22q, genes_exp_t, by = "Sample")  
TERT5_5p_merge <- merge(dados_TERT5_5p, genes_exp_t, by = "Sample")  
PAK1_11q_merge <- merge(dados_PAK1_11q, genes_exp_t, by = "Sample") 
PRAME_22q_merge <- merge (dados_PRAME_22q, genes_exp_t, by = "Sample")

glimpse(CCND1_11q_merge)  
# Transform columns from GAB2 to numeric  
CCND1_11q_merge <- CCND1_11q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

CDK4_12q_merge <- CDK4_12q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

CRKL_22q_merge <- CRKL_22q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

EP300_22q_merge <- EP300_22q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

GAB2_11q_merge <- GAB2_11q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

LZTR1_22q_merge <- LZTR1_22q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

MDM2_12q_merge <- MDM2_12q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

SOX10_22q_merge <- SOX10_22q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

TERT5_5p_merge <- TERT5_5p_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))  

PAK1_11q_merge<- PAK1_11q_merge %>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))

PRAME_22q_merge  <-PRAME_22q_merge%>%  
  mutate(across(GAB2:PRAME, ~ as.numeric(.)))

# testing normality  
shapiro_result <- shapiro.test(CCND1_11q_merge$CCND1) # not normal  
shapiro_result <- shapiro.test(CCND1_11q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(CDK4_12q_merge$CDK4) # not normal  
shapiro_result <- shapiro.test(CDK4_12q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(CRKL_22q_merge$CRKL) # not normal  
shapiro_result <- shapiro.test(CRKL_22q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(EP300_22q_merge$EP300) # normal  
shapiro_result <- shapiro.test(EP300_22q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(GAB2_11q_merge$GAB2) # normal  
shapiro_result <- shapiro.test(GAB2_11q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(LZTR1_22q_merge$LZTR1) # not normal  
shapiro_result <- shapiro.test(LZTR1_22q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(MDM2_12q_merge$MDM2) # not normal  
shapiro_result <- shapiro.test(MDM2_12q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(SOX10_22q_merge$SOX10) # normal  
shapiro_result <- shapiro.test(SOX10_22q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(TERT5_5p_merge$TERT) # not normal  
shapiro_result <- shapiro.test(TERT5_5p_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(PRAME_22q_merge$Copy_number) # not normal  
shapiro_result <- shapiro.test(PRAME_22q_merge$PRAME) # not normal 
shapiro_result <- shapiro.test(PAK1_11q_merge$Copy_number) # not normal 
shapiro_result <- shapiro.test(PAK1_11q_merge$PAK1)# not normal 

# transforming in Log2  
CCND1_11q_merge$CCND1 <- log2(CCND1_11q_merge$CCND1 + 1)  
CDK4_12q_merge$CDK4 <- log2(CDK4_12q_merge$CDK4 + 1)  
CRKL_22q_merge$CRKL <- log2(CRKL_22q_merge$CRKL + 1)  
EP300_22q_merge$EP300 <- log2(EP300_22q_merge$EP300 + 1)  
GAB2_11q_merge$GAB2 <- log2(GAB2_11q_merge$GAB2 + 1)  
LZTR1_22q_merge$LZTR1 <- log2(LZTR1_22q_merge$LZTR1 + 1)  
MDM2_12q_merge$MDM2 <- log2(MDM2_12q_merge$MDM2 + 1)  
SOX10_22q_merge$SOX10 <- log2(SOX10_22q_merge$SOX10 + 1)  
TERT5_5p_merge$TERT <- log2(TERT5_5p_merge$TERT + 1)  
PAK1_11q_merge$PAK1 <- log2(PAK1_11q_merge$PAK1 + 1)  
PRAME_22q_merge$PRAME <- log2(PRAME_22q_merge$PRAME + 1)  

#plotting
npg_colors <- pal_npg("nrc")(10) # NPG colors  
extra_colors <- c("#FF7F00", "#6A3D9A", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")  
my_palette <- c(npg_colors, extra_colors) # Total of 16 colors  

# correlation  
cor_result <- cor.test(CCND1_11q_merge$CCND1, CCND1_11q_merge$Copy_number, method = "spearman")  

# plotting  
# Plot definition  
p1 <- ggplot(CCND1_11q_merge, aes(x = Copy_number, y = CCND1, color = `11q`)) + #  
  geom_point(alpha = 0.5, aes(color = `11q`), size = 2) +  
  labs(x = "Copy Number CCND1", y = "log2 (CCND1 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 12.5), xlim = c(1, 16), expand = TRUE) +
  scale_y_continuous(breaks = seq(0, 12.5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 16, by = 5)) 

p1  

ggsave(
  filename = here("results", "CCND1_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p1,
  device = "pdf",
  width = 3.38,
  height = 3.14
)


cor_result <- cor.test(CDK4_12q_merge$CDK4, CDK4_12q_merge$Copy_number, method = "spearman")  
p2 <- ggplot(CDK4_12q_merge, aes(x = Copy_number, y = CDK4, color = `12q`)) + #  
  geom_point(alpha = 0.5, aes(color = `12q`), size = 2) +  
  labs(x = "Copy Number CDK4", y = "log2 (CDK4 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 13), xlim = c(1, 12.5), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 13, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 12.5, by = 5)) 

p2  

ggsave(
  filename = here("results", "CDK4_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p2,
  device = "pdf",
  width = 3.38,
  height = 3.14
)
cor_result <- cor.test(EP300_22q_merge$EP300, EP300_22q_merge$Copy_number, method = "spearman")  
p3 <- ggplot(EP300_22q_merge, aes(x = Copy_number, y = EP300, color = `22q`)) + #  
  geom_point(alpha = 0.5, aes(color = `22q`), size = 2) +  
  labs(x = "Copy Number EP300", y = "log2 (EP300+ 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 7.5), xlim = c(1, 16), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 7.5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 16, by = 5)) 

p3  


ggsave(
  filename = here("results", "EP300_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p3,
  device = "pdf",
  width = 3.38,
  height = 3.14
)
cor_result <- cor.test(CRKL_22q_merge$CRKL, CRKL_22q_merge$Copy_number, method = "spearman")  
p4 <- ggplot(CRKL_22q_merge, aes(x = Copy_number, y = CRKL, color = `22q`)) + #  
  geom_point(alpha = 0.5, aes(color = `22q`), size = 2) +  
  labs(x = "Copy Number CRKL", y = "log2 (CRKL + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 12.5), xlim = c(1, 21), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 12.5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 21, by = 5)) 

p4  

ggsave(
  filename = here("results", "CRKL_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p4,
  device = "pdf",
  width = 3.38,
  height = 3.14
)
cor_result <- cor.test(GAB2_11q_merge$GAB2, GAB2_11q_merge$Copy_number, method = "spearman")  
p5 <- ggplot(GAB2_11q_merge, aes(x = Copy_number, y = GAB2, color = `11q`)) + #  
  geom_point(alpha = 0.5, aes(color = `11q`), size = 2) +  
  labs(x = "Copy Number GAB2", y = "log2 (GAB2 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 15), expand = TRUE) +
  scale_y_continuous(breaks = seq(0, 10, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 15, by = 5))    

p5  

ggsave(
  filename = here("results", "GAB2_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p5,
  device = "pdf",
  width = 3.38,
  height = 3.14
)

cor_result <- cor.test(LZTR1_22q_merge$LZTR1, LZTR1_22q_merge$Copy_number, method = "spearman")  
p6 <- ggplot(LZTR1_22q_merge, aes(x = Copy_number, y = LZTR1, color = `22q`)) + #  
  geom_point(alpha = 0.5, aes(color = `22q`), size = 2) +  
  labs(x = "Copy Number LZTR1", y = "log2 (LZTR1 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 12.5), xlim = c(1, 21), expand = FALSE)+
  scale_y_continuous(breaks = seq(0, 12.5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 21, by = 5))    
p6 

ggsave(
  filename = here("results", "LZTR1_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p6,
  device = "pdf",
  width = 3.38,
  height = 3.14
)


cor_result <- cor.test(MDM2_12q_merge$MDM2, MDM2_12q_merge$Copy_number, method = "spearman")  
p7 <- ggplot(MDM2_12q_merge, aes(x = Copy_number, y = MDM2, color = `12q`)) + #  
  geom_point(alpha = 0.5, aes(color = `12q`), size = 2) +  
  labs(x = "Copy Number MDM2", y = "log2 (MDM2 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 12.5), xlim = c(1, 21), expand = FALSE)  +
  scale_y_continuous(breaks = seq(0, 12.5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 21, by = 5)) 
p7  

ggsave(
  filename = here("results", "MDM2_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p7,
  device = "pdf",
  width = 3.38,
  height = 3.14
)


cor_result <- cor.test(SOX10_22q_merge$SOX10, SOX10_22q_merge$Copy_number, method = "spearman")  
p8 <- ggplot(SOX10_22q_merge, aes(x = Copy_number, y = SOX10, color = `22q`)) + #  
  geom_point(alpha = 0.5, aes(color = `22q`), size = 2) +  
  labs(x = "Copy Number SOX10", y = "log2 (SOX10 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(1, 12.5), xlim = c(1, 20), expand = TRUE)  +
  scale_y_continuous(breaks = seq(0, 12.5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) 

p8  

ggsave(
  filename = here("results", "SOX10_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p8,
  device = "pdf",
  width = 3.38,
  height = 3.14
)


cor_result <- cor.test(TERT5_5p_merge$TERT, TERT5_5p_merge$Copy_number, method = "spearman")  
p9 <- ggplot(TERT5_5p_merge, aes(x = Copy_number, y = TERT, color = `5p`)) + #  
  geom_point(alpha = 0.5, aes(color = `5p`), size = 2) +  
  labs(x = "Copy Number TERT", y = "log2 (TERT + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(0, 5), xlim = c(0, 16), expand = TRUE) +
  scale_y_continuous(breaks = seq(0, 5, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 16, by = 5)) 

p9  

ggsave(
  filename = here("results", "TERT_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p9,
  device = "pdf",
  width = 3.38,
  height = 3.14
)


cor_result <- cor.test(PAK1_11q_merge$PAK1, PAK1_11q_merge$Copy_number, method = "spearman")  
p10 <- ggplot(PAK1_11q_merge, aes(x = Copy_number, y = PAK1, color = `11q`)) + #  
  geom_point(alpha = 0.5, aes(color = `11q`), size = 2) +  
  labs(x = "Copy Number PAK1", y = "log2 (PAK1 + 1)") +  
  theme_classic() +  
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +  
  geom_text(  
    x = Inf, y = -Inf, inherit.aes = FALSE,  
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",  
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),  
                  sep = ""  
    ),  
    hjust = 1.1, vjust = -1, size = 4  
  ) +  
  scale_color_manual(values = my_palette) + # Maintains the same color for each patient  
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 17), expand = TRUE) +
  scale_y_continuous(breaks = seq(0, 8, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 17, by = 5)) 

p10 

ggsave(
  filename = here("results", "PAK1_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p10,
  device = "pdf",
  width = 3.38,
  height = 3.14
)


p11 <- ggplot(PRAME_22q_merge, aes(x = Copy_number, y = PRAME, color = `22q`)) +
  geom_point(alpha = 0.5, aes(color = `22q`), size = 2) +
  labs(x = "Copy Number PRAME", y = "log2 (PRAME + 1)") +
  theme_classic() +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "darkgrey") +
  geom_text(
    x = Inf, y = -Inf, inherit.aes = FALSE,
    label = paste("Rho:", round(cor_result$estimate, 2), "\n",
                  "p-value:", format(cor_result$p.value, scientific = FALSE, 3),
                  sep = ""
    ),
    hjust = 1.1, vjust = -1, size = 4
  ) +
  scale_color_manual(values = my_palette) +
  coord_cartesian(ylim = c(0, 14), xlim = c(0, 20), expand = TRUE) +
  scale_y_continuous(breaks = seq(0, 14, by = 2.5)) +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) 
p11

ggsave(
  filename = here("results", "PRAME_Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  plot = p11,
  device = "pdf",
  width = 3.38,
  height = 3.14
)

#combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 3)  
#ggsave(
  #filename = here("results", "Spearman_correlation_copynumber_hailstorm_genes_selected.pdf"),
  #plot = combined_plot,
  #device = "pdf",
  #width = 9.11,
  #height = 9
#)


# Saving all packages used to create these plots
writeLines(
  capture.output(sessionInfo()),
  here("package_version", "versions_03_FigureS3a.txt")
)

