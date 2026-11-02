# Author: Annie Squiavinato
# 2026-10
#_______________________________________________________________________________
# Reviewer response, Figure R4: does copy-number loss of ATM or CHEK1 translate
# into lower expression? One tumour per patient, Loss vs Neutral, BRA cohort.
#
# Same design and the same metrics as MX_ATM_CHEK1_expression.R, so the two
# cohorts can be put side by side:
#   effect size   DESeq2 Wald on raw counts, design ~ Status, giving
#                 log2FoldChange and its standard error directly
#   interval      log2FC +/- 1.96 * lfcSE, back-transformed to a fold change
#   companion     Wilcoxon on the VST values, the same test in both cohorts
#   display       DESeq2 VST (blind), fitted per gene on the samples tested
#
# The effect size is NOT taken as a difference in VST values: VST is close to
# log2 at high counts but compressed at low ones, so 2^(VST difference) is not
# a fold change.
#
# Status is taken per gene. A sample that is Gain or Amp for one gene and
# Neutral for the other therefore enters one test and not the other, which is
# why the fit is run separately for each gene.
#
# The DNA -> RNA mapping is READ from the metadata, not derived by substituting
# PD for PR: the letter suffix is not the same in both assays.
#
# ERCC spike-in rows are removed before library sizes are computed; spike-ins
# are added per unit of input and would distort the normalisation.
#
# INPUTS  (under IN)
#   copy_number_BRA_adjPloidy_06_10.tsv
#   2775_3112_3743_Counts_HTSeq_STAR_ENSv103_ERCC_ENS_gene_name_IDS_QC_PASSED_KEEPONLY.tsv.gz
#   6591_2775_3112_3743_METADATA_PDX_from_Latin_America_RNAseq_for_release.tsv
#
# OUTPUTS (under OUT)
#   BRA_ATM_CHEK1_cn_vs_expression_data.tsv
#   BRA_ATM_CHEK1_effect_size.tsv
#   BRA_ATM_expression_by_CN_status.pdf
#   BRA_CHEK1_expression_by_CN_status.pdf
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from INCA_res/, the directory that holds renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggpubr)

IN  <- "data"
OUT <- "results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

F_CNA  <- file.path(IN, "copy_number_BRA_adjPloidy_06_10.tsv")
F_EXPR <- file.path(IN,
                    "2775_3112_3743_Counts_HTSeq_STAR_ENSv103_ERCC_ENS_gene_name_IDS_QC_PASSED_KEEPONLY.tsv.gz")
F_META <- file.path(IN,
                    "6591_2775_3112_3743_METADATA_PDX_from_Latin_America_RNAseq_for_release.tsv")

for (f in c(F_CNA, F_EXPR, F_META))
  if (!file.exists(f)) stop("missing input: ", f)

GENES      <- c(ATM = "ENSG00000149311", CHEK1 = "ENSG00000149554")
EXPR_LABEL <- "VST-normalised counts"

# One sample per patient, cellularity > 0.35. DNA identifiers.
#
# Read from the cohort file, not from a vector typed into the script. The
# hailstorm analysis reads the same file; while each script carried its own
# pasted list the two drifted apart and reported different Loss/Neutral
# splits for the same cohort.
#
# AM96 and AM98 belong to another project. They are already absent from the
# file; EXCLUDE is a guard in case they are ever added back.
EXCLUDE <- c("PD70974a", "PD70975a")

F_IDS  <- file.path(IN, "CNA_samples_one_ppatient_12_05_26.csv")
if (!file.exists(F_IDS)) stop("missing input: ", F_IDS)

# The CSV has shipped with a trailing ",," line, which read.csv turns into a
# row with ID = "". An empty id matches nothing and surfaces as an error with
# no name in it, so blanks are dropped explicitly.
ids    <- trimws(as.character(read.csv(F_IDS)$ID))
blank  <- is.na(ids) | !nzchar(ids)
if (any(blank)) message("blank rows dropped from the cohort file: ", sum(blank))

ids_ok <- setdiff(ids[!blank], EXCLUDE)
stopifnot(length(ids_ok) == 39, !anyDuplicated(ids_ok))
message("cohort: ", length(ids_ok), " samples, one per patient")

## --------------------------------------------------------- copy number
cna_all <- read_tsv(F_CNA, show_col_types = FALSE) |>
  filter(Sample %in% ids_ok, Gene %in% names(GENES))

cna <- cna_all |>
  filter(Alteration %in% c("Loss", "Neutral")) |>
  transmute(DNA = Sample, AM_id, Gene,
            Status = factor(Alteration, levels = c("Neutral", "Loss")))

message("copy-number status (DNA level, one sample per patient):")
print(table(cna$Gene, cna$Status))

# Gain, Amp and Del are a different question; mixing them in would change what
# is being tested, so they are excluded and reported.
dropped <- filter(cna_all, !Alteration %in% c("Loss", "Neutral"))
if (nrow(dropped)) {
  message("\nATM/CHEK1 rows not Loss or Neutral, excluded:")
  print(table(dropped$Gene, dropped$Alteration))
}

## ------------------------------------------------------- DNA -> RNA map
meta <- read_tsv(F_META, show_col_types = FALSE, guess_max = 10000)

need <- c("Sanger_DNA_ID", "Sanger_RNA_ID", "final_sample_name_used",
          "samp_type", "OK_to_analyse_RNA?", "SAMPLE_USE", "Case ID")
missing_col <- setdiff(need, names(meta))
if (length(missing_col))
  stop("columns missing from the metadata: ", paste(missing_col, collapse = ", "))

# Tumour rows only. PDX rows carry final_sample_name_used ending in _hum and
# belong to the fidelity analyses.
map <- meta |>
  transmute(DNA    = str_squish(Sanger_DNA_ID),
            RNA    = str_squish(final_sample_name_used),
            RNA_id = str_squish(Sanger_RNA_ID),
            LabID  = str_squish(`Case ID`),
            type   = str_squish(samp_type),
            ok_rna = str_squish(`OK_to_analyse_RNA?`),
            use    = str_squish(SAMPLE_USE)) |>
  filter(type == "Tumour", toupper(ok_rna) == "YES",
         !is.na(RNA), RNA != "NA", nzchar(RNA))

message("\nmetadata rows: ", nrow(meta), " | tumour rows, RNA OK: ", nrow(map))

# Technical duplicates (PR53330c2, PR53350e2, PR53357d2) would break
# one-tumour-per-patient. Keep the row whose Case ID does not say "duplicate";
# if that does not resolve it, stop rather than pick arbitrarily.
dup <- map |> group_by(DNA) |> filter(n() > 1) |> ungroup()
if (nrow(dup)) {
  message("\nDNA samples with more than one RNA library:")
  print(as.data.frame(dup[, c("DNA", "RNA", "LabID")]), row.names = FALSE)
  map <- map |>
    group_by(DNA) |>
    filter(n() == 1 | !str_detect(tolower(LabID), "duplicate")) |>
    ungroup()
  # dplyr:: qualified: S4Vectors (loaded with DESeq2) exports a count()
  # generic that masks dplyr::count and fails with "Argument 'x' is not a
  # vector: list".
  unresolved <- map |> dplyr::count(DNA) |> dplyr::filter(n > 1)
  if (nrow(unresolved))
    stop("duplicates not resolved for: ", paste(unresolved$DNA, collapse = ", "))
  message("duplicates dropped; one RNA library per DNA sample kept")
}
stopifnot(!anyDuplicated(map$DNA))

no_rna <- setdiff(ids_ok, map$DNA)
if (length(no_rna))
  message("\ncohort samples with no usable RNA (dropped): ", length(no_rna),
          " -> ", paste(no_rna, collapse = ", "))

map <- filter(map, DNA %in% ids_ok)
message("cohort samples with RNA: ", nrow(map))
stopifnot(nrow(map) > 0)

## ---------------------------------------------------------- expression
expr <- read_tsv(F_EXPR, show_col_types = FALSE, guess_max = 10000)
message("\ncounts matrix: ", nrow(expr), " rows x ", ncol(expr), " columns")

# The gene id column is found by CONTENT: the header varies between exports and
# an assumed name silently returns nothing.
i_id <- which(vapply(expr, \(x)
                     mean(grepl("^ENSG\\d+", as.character(x))) > 0.5, logical(1)))[1]
if (is.na(i_id)) stop("no column of ENSG identifiers in ", basename(F_EXPR))
message("gene id column: '", names(expr)[i_id], "'")

ens <- sub("\\..*$", "", as.character(expr[[i_id]]))   # drops version suffix

no_gene <- names(GENES)[!GENES %in% ens]
if (length(no_gene))
  stop("genes not found in the matrix: ", paste(no_gene, collapse = ", "))

not_in_matrix <- setdiff(map$RNA, names(expr))
if (length(not_in_matrix)) {
  message("\nRNA ids in the metadata but not in the matrix (dropped): ",
          length(not_in_matrix), " -> ", paste(not_in_matrix, collapse = ", "))
  map <- filter(map, RNA %in% names(expr))
}
stopifnot(nrow(map) > 0)
message("samples carried into the analysis: ", nrow(map))

# ERCC rows identified by the id column, so a spike-in cannot leak into the
# library size.
is_ercc <- grepl("^ERCC-", as.character(expr[[i_id]]))
message("ERCC spike-in rows removed: ", sum(is_ercc))

M <- as.matrix(expr[!is_ercc, map$RNA, drop = FALSE])
rownames(M) <- ens[!is_ercc]
storage.mode(M) <- "double"
stopifnot(!any(is.na(M)))

# A tsv round-trip can deliver integral counts as double; rounding is a no-op
# then, and the check catches a matrix that is genuinely not counts.
if (max(abs(M - round(M))) > 1e-6)
  stop("the matrix is not integer counts -- the Wald fit and VST do not apply")
storage.mode(M) <- "integer"

if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is required")

## ------------------------------------------------------- per gene: test
run_gene <- function(g) {
  st <- cna |> filter(Gene == g) |>
    inner_join(select(map, DNA, RNA), by = "DNA") |>
    select(RNA, Status)
  st   <- st[match(colnames(M), st$RNA), ]
  keep <- !is.na(st$Status)
  Mg   <- M[, keep, drop = FALSE]
  sg   <- droplevels(st$Status[keep])
  
  n <- table(sg)
  message("\n", g, ": ", paste(names(n), n, sep = " = ", collapse = ", "))
  if (length(n) < 2 || min(n) < 3)
    return(tibble(gene = g, n_neutral = NA_integer_, n_loss = NA_integer_,
                  testable = FALSE))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = Mg,
    colData   = data.frame(row.names = colnames(Mg), Status = sg),
    design    = ~ Status)
  
  # All-zero genes destabilise the dispersion fit. The target gene is asserted
  # to survive: dropping a pre-specified gene for low expression would remove
  # the effect under test.
  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  stopifnot(GENES[[g]] %in% rownames(dds))
  
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  r   <- DESeq2::results(dds, contrast = c("Status", "Loss", "Neutral"))
  id  <- GENES[[g]]
  
  # vst() needs enough genes for its parametric fit; with few samples it can
  # fail, and varianceStabilizingTransformation is the exact, slower form.
  vsd <- tryCatch(DESeq2::vst(dds, blind = TRUE),
                  error = function(e) {
                    message("vst() failed (", conditionMessage(e),
                            ") -- using varianceStabilizingTransformation()")
                    DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
                  })
  e <- SummarizedExperiment::assay(vsd)[id, ]
  w <- stats::wilcox.test(e ~ sg, conf.int = TRUE, exact = FALSE)
  
  assign("VST_LONG", bind_rows(
    if (exists("VST_LONG", envir = .GlobalEnv)) get("VST_LONG", envir = .GlobalEnv),
    tibble(Gene = g, RNA = names(e), expr = unname(e), Status = sg)),
    envir = .GlobalEnv)
  
  tibble(gene = g, n_neutral = n[["Neutral"]], n_loss = n[["Loss"]],
         testable = TRUE,
         log2FC = r[id, "log2FoldChange"], lfcSE = r[id, "lfcSE"],
         log2FC_lo = r[id, "log2FoldChange"] - 1.96 * r[id, "lfcSE"],
         log2FC_hi = r[id, "log2FoldChange"] + 1.96 * r[id, "lfcSE"],
         p_wald = r[id, "pvalue"], p_wilcox = w$p.value)
}

if (exists("VST_LONG", envir = .GlobalEnv)) rm("VST_LONG", envir = .GlobalEnv)
res <- bind_rows(lapply(names(GENES), run_gene)) |>
  mutate(cohort = "BRA", fold_change = 2^log2FC,
         fc_lo = 2^log2FC_lo, fc_hi = 2^log2FC_hi,
         estimator = "DESeq2 Wald on raw counts",
         normalisation = EXPR_LABEL,
         cn_table      = basename(F_CNA),
         ploidy_source = "mean_ploidy_adj") |>
  relocate(cohort, gene, n_neutral, n_loss)

print(as.data.frame(res), row.names = FALSE, digits = 3)
write_tsv(res, file.path(OUT, "BRA_ATM_CHEK1_effect_size.tsv"))
write_tsv(VST_LONG, file.path(OUT, "BRA_ATM_CHEK1_cn_vs_expression_data.tsv"))

## ----------------------------------------------------------------- plot
# The caption carries the effect size and its interval, not only the p value:
# with these group sizes a non-significant p says little on its own, whereas
# the interval shows how large an effect the data still allow.
for (g in names(GENES)) {
  d <- filter(VST_LONG, Gene == g)
  if (!nrow(d)) next
  s <- filter(res, gene == g)
  
  cap <- sprintf("n = %d (%d neutral, %d loss)", nrow(d),
                 sum(d$Status == "Neutral"), sum(d$Status == "Loss"))
  if (nrow(s) == 1 && isTRUE(s$testable))
    cap <- paste0(cap, sprintf(
      "\nfold change %.2f (95%% CI %.2f to %.2f)\nP = %.3f (Wilcoxon), P = %.3f (Wald)",
      s$fold_change, s$fc_lo, s$fc_hi, s$p_wilcox, s$p_wald))
  
  p <- ggplot(d, aes(Status, expr, fill = Status)) +
    geom_boxplot(alpha = .7, width = .6, outlier.shape = NA) +
    geom_jitter(width = .15, size = 3, shape = 21, colour = "black", alpha = .8) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 4.5) +
    scale_fill_manual(values = c(Neutral = "#4E79A7", Loss = "#E15759")) +
    labs(title = paste(g, "expression, BRA cohort"),
         subtitle = "one tumour per patient, cellularity > 0.35",
         x = "Copy-number status", y = EXPR_LABEL,
         caption = cap) +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(hjust = .5, face = "bold"),
          plot.subtitle = element_text(hjust = .5),
          plot.caption  = element_text(hjust = 0, size = 8),
          legend.position = "none")
  
  ggsave(file.path(OUT, paste0("BRA_", g, "_expression_by_CN_status.pdf")),
         p, width = 4.2, height = 5.2)
}

cat("\nfor the text:\n")
for (i in seq_len(nrow(res)))
  if (isTRUE(res$testable[i]))
    cat(sprintf("  BRA, %s: fold change %.2f (95%% CI %.2f to %.2f), P = %.3f; n = %d vs %d\n",
                res$gene[i], res$fold_change[i], res$fc_lo[i], res$fc_hi[i],
                res$p_wilcox[i], res$n_neutral[i], res$n_loss[i]))

message("\ndone")

