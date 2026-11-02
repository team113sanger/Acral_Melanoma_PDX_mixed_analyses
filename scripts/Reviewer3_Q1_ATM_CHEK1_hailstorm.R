
# Author: Annie Squiavinato
# 2026-10
#_______________________________________________________________________________
# ATM and CHEK1 copy-number loss in the BRA cohort, two questions:
#   A  does loss co-occur with a hailstorm on 5p, 11q or 22q?
#  
#
# Both genes are tested on each arm, and correction is applied within the arm,
# across the two genes. The family is therefore one arm by two genes: three
# families of two, not one of six.
#
# COHORT
# The sample list is read ONCE, from CNA_samples_one_ppatient_12_05_26.csv, and
# used by both analyses and by BRA_ATM_CHEK1_expression.R. While each analysis
# carried its own pasted vector the lists drifted apart and reported different
# Loss/Neutral splits for the same cohort.
#
# PLOIDY
# copy_number_BRA_adjPloidy_06_10.tsv uses mean_ploidy_adj from
# Adjusted_ploidy_table.tsv, the manually reviewed ploidy. 
#
# HAILSTORM ABSENCE
# The hailstorm file lists detected events only: a sample with no hailstorm has
# no row. Absence is therefore defined by the cohort list, which is valid only
# if every listed sample was screened.
#
# HAILSTORM IS A PATIENT-LEVEL ATTRIBUTE
# Positive for an arm if detected in any sample from that patient. Sanger ids
# encode the donor in the numeric part and the sample in the trailing letter
# (PD70967a and PD70967c are the same patient). Without this, five patients
# whose hailstorm was recorded on a sibling sample become false negatives.
#
# INPUTS  (under IN)
#   CNA_samples_one_ppatient_12_05_26.csv
#   copy_number_BRA_adjPloidy_06_10.tsv
#   Supplementary_Detected_hailstorms_in_all_samples_up_23_09_26.tsv
#
# OUTPUTS (under OUT)
#   ATM_CHEK1_vs_hailstorm_fisher_<tag>.tsv
#   ATM_CHEK1_vs_hailstorm_burden_<tag>.tsv
#   
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from INCA_res/, the directory that holds renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(tidyr)
library(readr)
library(tibble)

IN  <- "data"
OUT <- "results"
if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DATA_TAG      <- "adjPloidy_06_10"
PLOIDY_SOURCE <- "mean_ploidy_adj"

F_IDS  <- file.path(IN, "CNA_samples_one_ppatient_12_05_26.csv")
F_CNA  <- file.path(IN, paste0("copy_number_BRA_", DATA_TAG, ".tsv"))
F_HAIL <- file.path(IN,
                    "Supplementary_Detected_hailstorms_in_all_samples_up_23_09_26.tsv")
for (f in c(F_IDS, F_CNA, F_HAIL))
  if (!file.exists(f)) stop("missing input: ", f)

# The three arms most frequently affected by hailstorms in this cohort. 5p is
# included as a negative control of sorts: neither ATM nor CHEK1 is on it, so
# an association there would point to overall structural burden rather than to
# anything local.
ARMS    <- c("5p", "11q", "22q")
DRIVERS <- c("ATM", "CHEK1")
TARGETS <- c("GAB2","CCND1","FADD","PAK1","NUMA1","CDK4","PTPRB","MDM2",
             "LRIG3","HMGA2","EP300","SOX10","MRTFA","CRKL","BCR","TERT",
             "RICTOR","LIFR","GOLPH3")

pid <- function(x) sub("[a-z]$", "", sub("_hum$", "", x))

## ----------------------------------------------------------------- cohort
# This block is identical to the one in BRA_ATM_CHEK1_expression.R, down to
# the variable name, so the two analyses cannot diverge.
#
# AM96 and AM98 belong to another project. They are already absent from the
# file; EXCLUDE is a guard in case they are ever added back.
EXCLUDE <- c("PD70974a", "PD70975a")

# The CSV has shipped with a trailing ",," line, which read.csv turns into a
# row with ID = "". An empty id matches nothing and surfaces as an error with
# no name in it, so blanks are dropped explicitly.
ids    <- trimws(as.character(read.csv(F_IDS)$ID))
blank  <- is.na(ids) | !nzchar(ids)
if (any(blank)) message("blank rows dropped from the cohort file: ", sum(blank))

ids_ok <- setdiff(ids[!blank], EXCLUDE)
stopifnot(length(ids_ok) == 39, !anyDuplicated(ids_ok))
message("cohort: ", length(ids_ok), " samples, one per patient")

## ------------------------------------------------------------------ input
cna  <- read_tsv(F_CNA,  show_col_types = FALSE)
hail_raw <- read_tsv(F_HAIL, show_col_types = FALSE)
stopifnot(all(c("Sample", "AM_id", "Gene", "Alteration") %in% names(cna)))

missing <- setdiff(ids_ok, cna$Sample)
if (length(missing)) stop("cohort sample(s) with no copy-number data: ",
                          paste(missing, collapse = ", "))
missing <- setdiff(c(DRIVERS, TARGETS), unique(cna$Gene))
if (length(missing)) message("gene(s) absent from the copy-number table: ",
                             paste(missing, collapse = ", "))

# In the hailstorm file `Sample` holds the AM identifier and `Sanger_RNA_ID`
# the PD identifier; the PD one matches the copy-number and cohort tables.
unknown <- setdiff(pid(hail_raw$Sanger_RNA_ID), pid(ids_ok))
message("patients with hailstorms but outside the cohort: ", length(unknown),
        if (length(unknown)) paste0(" -> ", paste(unknown, collapse = ", ")) else "")

## ============================== A. hailstorm ==============================
hs_pat <- hail_raw |> mutate(patient = pid(Sanger_RNA_ID)) |> distinct(patient, Arm)

hail <- expand_grid(Sample = ids_ok, Arm = ARMS) |>
  mutate(patient = pid(Sample)) |>
  left_join(mutate(hs_pat, hit = TRUE), by = c("patient", "Arm")) |>
  mutate(hailstorm = factor(if_else(is.na(hit), "Absent", "Present"),
                            levels = c("Present", "Absent"))) |>
  select(Sample, Arm, hailstorm)

stopifnot(nrow(hail) == length(ids_ok) * length(ARMS))
cat("\nhailstorm status per arm:\n"); print(table(hail$Arm, hail$hailstorm))

burden <- tibble(Sample = ids_ok) |>
  mutate(patient = pid(Sample)) |>
  left_join(dplyr::count(hs_pat, patient, name = "n_hailstorms"), by = "patient") |>
  mutate(n_hailstorms = replace_na(n_hailstorms, 0L)) |>
  select(Sample, n_hailstorms)

cat("\nhailstorm burden per patient:\n"); print(table(burden$n_hailstorms))

# Loss vs Neutral only. A gained or amplified locus is neither, so those
# samples are dropped from that gene's test rather than folded into a
# "not lost" group. With 39 samples, losing two matters: it is printed.
gene_alt <- filter(cna, Sample %in% ids_ok, Gene %in% DRIVERS)

dropped <- filter(gene_alt, !Alteration %in% c("Loss", "Neutral"))
if (nrow(dropped)) {
  cat("\nexcluded, neither Loss nor Neutral:\n")
  print(as.data.frame(dropped[, c("Sample", "Gene", "Alteration", "Copy_number")]),
        row.names = FALSE)
}

df <- gene_alt |>
  filter(Alteration %in% c("Loss", "Neutral")) |>
  mutate(alt_bin = factor(Alteration, levels = c("Loss", "Neutral"))) |>
  inner_join(hail, by = "Sample", relationship = "many-to-many") |>
  left_join(burden, by = "Sample")

cat("\nLoss / Neutral per gene (cohort level):\n")
print(as.data.frame(distinct(df, Gene, Sample, alt_bin) |>
                      dplyr::count(Gene, alt_bin)), row.names = FALSE)

# MULTIPLE TESTING, AND WHY THE FAMILY MATTERS
# Correction is applied WITHIN each arm, across the two genes tested on it, so
# the family is "one arm, the genes we looked at": three families of two. The
# alternative, correcting over all six tests at once, is reported alongside.
#
# The two are written out because the choice changes which associations clear
# 0.05, and it was made after the raw p values were known. Burying it would be
# the problem; showing both is not.
#
# Bonferroni is kept alongside BH but is conservative here: ATM and CHEK1 sit
# ~17 Mb apart on the same arm and their loss is correlated, so the tests are
# not independent.
results <- df |>
  group_by(Gene, Arm) |>
  group_modify(~{
    ct <- table(.x$hailstorm, .x$alt_bin)
    ft <- fisher.test(ct)
    tibble(n_samples = nrow(.x),
           hail_Loss      = ct["Present", "Loss"],
           hail_Neutral   = ct["Present", "Neutral"],
           nohail_Loss    = ct["Absent",  "Loss"],
           nohail_Neutral = ct["Absent",  "Neutral"],
           zero_cell  = any(ct == 0),
           odds_ratio = unname(ft$estimate),
           CI_low     = ft$conf.int[1],
           CI_high    = ft$conf.int[2],
           p_value    = ft$p.value)
  }) |>
  ungroup() |>
  # Correction within each gene, across the arms tested for it. ATM and CHEK1
  # are two independent analyses and are never pooled: no p value from one
  # gene enters the other's adjustment.
  group_by(Gene) |>
  mutate(padj_per_gene = p.adjust(p_value, method = "BH"),
         bonf_per_gene = p.adjust(p_value, method = "bonferroni")) |>
  ungroup() |>
  mutate(bh_family     = paste0("within gene (", length(ARMS), " arms)"),
         ploidy_source = PLOIDY_SOURCE,
         cn_table      = basename(F_CNA)) |>
  relocate(Gene, Arm, n_samples, p_value, padj_per_gene, bonf_per_gene) |>
  arrange(Gene, p_value)

cat("\n===== Fisher, ", nrow(results), " pre-specified tests =====\n", sep = "")
print(as.data.frame(results), row.names = FALSE, digits = 4)
write_tsv(results, file.path(OUT,
                             paste0("ATM_CHEK1_vs_hailstorm_fisher_", DATA_TAG, ".tsv")))

# A patient with hailstorms elsewhere already carries a higher structural load,
# so an arm-specific association can be confounded by total burden.
burden_test <- df |>
  distinct(Gene, Sample, alt_bin, n_hailstorms) |>
  group_by(Gene) |>
  group_modify(~ tibble(
    n_Loss      = sum(.x$alt_bin == "Loss"),
    n_Neutral   = sum(.x$alt_bin == "Neutral"),
    med_Loss    = median(.x$n_hailstorms[.x$alt_bin == "Loss"]),
    med_Neutral = median(.x$n_hailstorms[.x$alt_bin == "Neutral"]),
    p_value     = if (min(table(.x$alt_bin)) >= 3)
      wilcox.test(n_hailstorms ~ alt_bin, data = .x, exact = FALSE)$p.value
    else NA_real_)) |>
  ungroup() |>
  mutate(ploidy_source = PLOIDY_SOURCE, cn_table = basename(F_CNA))

cat("\n===== sensitivity: hailstorm burden vs gene loss =====\n")
print(as.data.frame(burden_test), row.names = FALSE, digits = 4)
write_tsv(burden_test, file.path(OUT,
                                 paste0("ATM_CHEK1_vs_hailstorm_burden_", DATA_TAG, ".tsv")))

