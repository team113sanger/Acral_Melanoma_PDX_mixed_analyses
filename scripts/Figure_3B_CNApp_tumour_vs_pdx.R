### Author: Annie Squiavinato
### 2025-04-25, revised 2026-10
#_______________________________________________________________________________
# Figure 3B: CNApp copy-number scores (FCS, BCS, GCS) of patient tumours
# against their matched PDX-X1, by Spearman correlation.
#
# INPUTS  (under IN)
#   Cnapp_Scores_all_Samples_25_09_26.xlsx
#       the first row is a title; the real header is on the second, hence
#       skip = 1. Without it every score column is read as text.
#   metadata_copy_number_data_filtered_by_cellularity_all_samples.csv
#
# OUTPUTS (under OUT)
#   Correlation_FCS_tumour_vs_pdx_<DATA_TAG>.pdf   (also BCS, GCS)
#   metadata_pdx_tumor_pairs_<DATA_TAG>.csv
#   CNApp_tumour_vs_pdx_correlation_<DATA_TAG>.tsv
#
# Post-processing in Affinity for the submitted figure
# Package versions are managed with renv: run renv::restore() once.
# Run this script from INCA_res/, the directory that holds renv.lock.
#_______________________________________________________________________________

library(tidyverse)
library(readxl)
library(ggsci)

IN  <- "data"
OUT <- "results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DATA_TAG        <- "25_09_26"
MIN_CELL        <- 0.35
ONE_PER_PATIENT <- FALSE     # FALSE = 30 pairs / 23 patients; TRUE = 23 / 23
COLOUR_BY       <- "lesion"  # lesion, tissue, patient or none

SCORES_FILE <- file.path(IN, "Cnapp_Scores_all_Samples_25_09_26.xlsx")
META_FILE   <- file.path(IN,
  "metadata_copy_number_data_filtered_by_cellularity_all_samples.csv")

for (f in c(SCORES_FILE, META_FILE))
  if (!file.exists(f)) stop("missing input: ", f)

## ================================================================ input
meta       <- read_csv(META_FILE, show_col_types = FALSE)
scores_raw <- read_excel(SCORES_FILE, skip = 1)
message("columns read: ", paste(names(scores_raw), collapse = " | "))

# The PD identifier is located by CONTENT, not by column name: the name varies
# between exports (Sanger_id, Sanger ID, Sanger_RNA_ID) but the PDxxxxx format
# does not.
find_col <- function(df, pattern) {
  ok <- vapply(df, function(x)
    mean(grepl(pattern, trimws(as.character(x)))) > 0.8, logical(1))
  if (!any(ok)) NA_character_ else names(df)[which(ok)[1]]
}

pd_col <- find_col(scores_raw, "^PD\\d+")
am_col <- find_col(scores_raw, "^AM\\d+")
if (is.na(pd_col))
  stop("no PD identifier column in the score file. Columns: ",
       paste(names(scores_raw), collapse = ", "))
message("PD identifier: '", pd_col, "' | AM: '", am_col, "'")

missing_score <- setdiff(c("FCS", "BCS", "GCS"), names(scores_raw))
if (length(missing_score))
  stop("score column(s) absent after skip = 1: ",
       paste(missing_score, collapse = ", "), ". Check the title row.")

scores <- scores_raw |>
  transmute(Sample = trimws(as.character(.data[[pd_col]])),
            FCS = suppressWarnings(as.numeric(FCS)),
            BCS = suppressWarnings(as.numeric(BCS)),
            GCS = suppressWarnings(as.numeric(GCS))) |>
  filter(!is.na(Sample), nzchar(Sample))

# as.numeric() on text returns NA silently.
n_na <- colSums(is.na(select(scores, FCS, BCS, GCS)))
if (any(n_na > 0)) {
  message("non-numeric values became NA -> ",
          paste(names(n_na), n_na, sep = "=", collapse = ", "))
  print(head(filter(scores, is.na(FCS) | is.na(BCS) | is.na(GCS))))
}
message("scores read: ", nrow(scores), " samples")
stopifnot(!anyDuplicated(scores$Sample))

## =============================================================== pairing
m <- meta |>
  mutate(Cellularity = as.numeric(gsub(",", ".", Cellularity)),
         Ploidy      = as.numeric(gsub(",", ".", Ploidy)),
         AM_label    = Lab_id |>
           str_remove("_(gDNA|RNA|gRNA)_tumour$") |>
           str_replace_all("_X", "-X"),
         AM_patient  = str_extract(Lab_id, "^AM\\d+"))

n0 <- nrow(m)
m  <- filter(m, Cellularity > MIN_CELL)
message("cellularity > ", MIN_CELL, ": ", nrow(m), "/", n0,
        " samples (lowest kept: ", min(m$Cellularity), ")")

tum <- m |>
  filter(Sample_type == "Patient") |>
  select(AM_patient, AM_tum = AM_label, PD_tum = Sample,
         Tissue_tum = Tissue_histology, Cell_tum = Cellularity,
         Ploidy_tum = Ploidy)

pdx <- m |>
  filter(Sample_type == "PDX-X1") |>
  mutate(AM_tum = str_remove(AM_label, "-X1$")) |>
  select(AM_tum, AM_pdx = AM_label, PD_pdx = Sample,
         Cell_pdx = Cellularity, Ploidy_pdx = Ploidy)

pairs_df   <- inner_join(tum, pdx, by = "AM_tum")
pdx_orphan <- anti_join(pdx, tum, by = "AM_tum")
tum_no_pdx <- anti_join(tum, pdx, by = "AM_tum")

message("\npairs formed          : ", nrow(pairs_df))
message("distinct patients     : ", n_distinct(pairs_df$AM_patient))
message("PDX-X1 with no tumour : ", nrow(pdx_orphan),
        if (nrow(pdx_orphan))
          paste0(" (", paste(pdx_orphan$AM_pdx, collapse = ", "), ")") else "")
message("tumours with no PDX-X1: ", nrow(tum_no_pdx))

dup <- pairs_df |> dplyr::count(AM_patient) |> filter(n > 1)
if (nrow(dup))
  message("patients with >1 paired lesion: ",
          paste0(dup$AM_patient, " (", dup$n, ")", collapse = ", "))

paired <- pairs_df |>
  left_join(rename_with(scores, ~ paste0(.x, "_tum"), -Sample),
            by = c("PD_tum" = "Sample")) |>
  left_join(rename_with(scores, ~ paste0(.x, "_pdx"), -Sample),
            by = c("PD_pdx" = "Sample")) |>
  rename(Patients = AM_patient) |>
  relocate(Patients, AM_tum, AM_pdx, PD_tum, PD_pdx)

no_score <- paired |> filter(if_any(matches("^(FCS|BCS|GCS)_(tum|pdx)$"), is.na))
if (nrow(no_score)) {
  message("\n", nrow(no_score), " pair(s) with a missing score:")
  print(as.data.frame(select(no_score, Patients, PD_tum, PD_pdx,
                             matches("^(FCS|BCS|GCS)_"))), row.names = FALSE)
}

if (ONE_PER_PATIENT) {
  # Explicit criterion: the lesion with the highest tumour cellularity, not
  # "the first", which would depend on row order.
  paired <- paired |>
    group_by(Patients) |>
    slice_max(Cell_tum, n = 1, with_ties = FALSE) |>
    ungroup()
  message("\nONE_PER_PATIENT = TRUE -> ", nrow(paired), " pairs")
}

stopifnot(nrow(paired) > 0, !anyDuplicated(paired$PD_pdx))

write_csv(paired, file.path(OUT, paste0("metadata_pdx_tumor_pairs_", DATA_TAG, ".csv")))

cat("\n== pairs ==\n")
print(as.data.frame(select(paired, Patients, AM_tum, AM_pdx, PD_tum, PD_pdx,
                           matches("^(FCS|BCS|GCS)_"))),
      row.names = FALSE, digits = 4)

## =============================================================== colours
# No colour repeats. The npg + extras palette held 16 colours for 23 patients,
# so ggplot would recycle silently and two patients would look identical.
paired <- paired |>
  mutate(Lesion = AM_tum,
         Tissue = Tissue_tum |>
           str_remove("^Patient") |>
           str_remove("^Tumour_") |>
           str_replace_all("_", " ") |>
           str_to_sentence())

colour_var <- switch(COLOUR_BY,
                     lesion = "Lesion", tissue = "Tissue",
                     patient = "Patients", none = NA_character_,
                     stop("COLOUR_BY must be lesion, tissue, patient or none"))

if (!is.na(colour_var)) {
  levs       <- sort(unique(paired[[colour_var]]))
  my_palette <- pal_igv("default")(max(length(levs), 3))[seq_along(levs)]
  names(my_palette) <- levs
  message("\ncolouring by '", colour_var, "': ", length(levs), " levels")
  print(table(paired[[colour_var]]))
  stopifnot(length(my_palette) == length(levs),
            !anyDuplicated(my_palette), !any(is.na(my_palette)))
} else {
  my_palette <- NULL
}

colour_lab <- switch(COLOUR_BY, lesion = "Lesion", tissue = "Tissue",
                     patient = "Patient", none = NULL)

## ================================================================= plots
plot_score <- function(df, score, ylim = NULL, xlim = NULL) {
  xv <- paste0(score, "_tum"); yv <- paste0(score, "_pdx")
  if (!all(c(xv, yv) %in% names(df))) return(NULL)

  d <- df[!is.na(df[[xv]]) & !is.na(df[[yv]]), ]
  if (nrow(d) < nrow(df))
    message(score, ": ", nrow(df) - nrow(d), " pair(s) dropped for NA")
  if (nrow(d) < 3) { message(score, ": too few pairs"); return(NULL) }

  # Spearman assumes monotonicity only. exact = FALSE because with ties no
  # exact p exists.
  ct  <- cor.test(d[[xv]], d[[yv]], method = "spearman", exact = FALSE)
  lab <- paste0("Rho = ", sprintf("%.2f", ct$estimate),
                "\np = ",  format.pval(ct$p.value, digits = 2, eps = 1e-4),
                "\nn = ",  nrow(d))
  message(score, ": rho = ", sprintf("%.3f", ct$estimate),
          "  p = ", signif(ct$p.value, 3), "  n = ", nrow(d))

  p <- ggplot(d, aes(x = .data[[xv]], y = .data[[yv]])) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                colour = "grey35", fill = "grey80", alpha = 0.4)

  if (is.na(colour_var)) {
    p <- p + geom_point(shape = 17, size = 4, alpha = 0.75, colour = "grey25")
  } else {
    # my_palette is NAMED, so a level keeps its colour whatever the order and
    # the three panels stay consistent with each other.
    p <- p +
      geom_point(aes(colour = .data[[colour_var]]),
                 shape = 17, size = 4, alpha = 0.8) +
      scale_colour_manual(values = my_palette, drop = FALSE)
  }

  p <- p +
    annotate("text", x = Inf, y = -Inf, label = lab,
             hjust = 1.05, vjust = -0.25, size = 3.6, lineheight = 1.05) +
    labs(x = paste(score, "patient tumours"), y = paste(score, "AM-PDXs"),
         colour = colour_lab) +
    theme_classic(base_size = 13)

  if (!is.na(colour_var) && length(my_palette) > 12)
    p <- p + guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3)))

  # coord_cartesian only crops the view; scale_*_continuous(limits=) would drop
  # points and change the fitted line silently.
  if (!is.null(ylim) || !is.null(xlim))
    p <- p + coord_cartesian(ylim = ylim, xlim = xlim, expand = TRUE)
  if (!is.null(ylim) && any(d[[yv]] < ylim[1] | d[[yv]] > ylim[2]))
    warning(score, ": point(s) outside ylim -- hidden, but still in rho")

  list(plot = p, test = ct, n = nrow(d))
}

res <- list(BCS = plot_score(paired, "BCS", ylim = c(0, 100)),
            FCS = plot_score(paired, "FCS", ylim = c(0, 1250)),
            GCS = plot_score(paired, "GCS"))
res <- Filter(Negate(is.null), res)

for (s in names(res))
  ggsave(file.path(OUT, paste0("Correlation_", s, "_tumour_vs_pdx_", DATA_TAG, ".pdf")),
         res[[s]]$plot, device = "pdf", width = 5.52, height = 5.5)

## ------------------------------------------- statistics for the legend
summary_tbl <- bind_rows(lapply(res, function(r) tibble(
  n = r$n, rho = unname(r$test$estimate), p_value = r$test$p.value,
  method = "Spearman")), .id = "score") |>
  mutate(one_per_patient = ONE_PER_PATIENT,
         n_patients      = n_distinct(paired$Patients))

print(as.data.frame(summary_tbl), row.names = FALSE, digits = 4)
write_tsv(summary_tbl,
          file.path(OUT, paste0("CNApp_tumour_vs_pdx_correlation_", DATA_TAG,
                                if (ONE_PER_PATIENT) "_1perpatient" else "", ".tsv")))
