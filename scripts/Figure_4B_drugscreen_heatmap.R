# Author: Flavia Aguiar, Danielle Carvalho and Annie Squiavinato
# 2025-04-25, revised 2026-11-02
#_______________________________________________________________________________
# Drug-screen heatmap: log10 IC50 of 30 small-molecule inhibitors across four
# acral melanoma cell lines, grouped by functional class.
#
# WHITE CELLS
# White means one thing only: no IC50 was determined, i.e. viability did not
# fall by 50% at any concentration assayed. These arrive as NA in the input
# file; the script creates none of its own.
#
# A previous version also blanked out IC50 values that extrapolated above the
# highest dose, but only for avapritinib. That rule was dropped, for two
# reasons. It was applied to one compound out of thirty, so XAV939 in AM007
# (15.78 uM against a 10 uM ceiling, the same situation) stayed on the plot
# while avapritinib did not. And a high extrapolated IC50 is informative: on
# this scale blue means a large IC50, which is the correct reading for a
# compound that did not work. Blanking it discarded that.
#
# EXTRAPOLATED VALUES
# Nothing is censored now, so some plotted values lie outside the range
# assayed for their compound. The audit block below measures how far outside,
# in dilution steps of 3.162 (the spacing actually used between doses), and
# writes the result next to the figure. Seven of the eleven sit less than half
# a step beyond the last dose tested, which is interpolation to the edge rather
# than extrapolation. Entinostat is the exception at three to four steps in all
# four lines, and is named in the figure legend for that reason.
#
# The audit does not change the figure. It exists so that a future update to
# the input file cannot introduce a new extreme case unnoticed.
#
# INPUTS
#   data/IC50_acral_drugs_annotations.csv
#       rows = compounds, col 1 = functional class, remaining cols = cell lines
#       decimal separator is a comma
#
# OUTPUTS
#   results/drugscreen_heatmap.pdf
#   results/drugscreen_IC50_log10.tsv          matrix as plotted, white as NA
#   results/drugscreen_outside_dose_range.tsv  audit of extrapolated values
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from the INCA_res/ directory, the one holding renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

IN  <- "data"
OUT <- "results"
if (!dir.exists(IN))
  stop("'data' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CELL_LINES <- c("AM007", "AM016", "AM015a", "AM015b")

ROW_GROUPS <- list(
  "RTK /\nMAPK"            = c("Cobimetinib", "Trametinib", "Ripretinib",
                               "Avapritinib"),
  "b-catenin"              = c("XAV939"),
  "Cell cycle"             = c("GSK461364", "Volasertib", "MLN8054",
                               "GSK1070916", "BAY-320", "Abemaciclib",
                               "Palbociclib", "Dinaciclib", "Nutlin-3a"),
  "DDR &\nRepair"          = c("KU-60019", "Berzosertib", "Ceralasertib",
                               "CHIR-124", "SCH900776", "Adavosertib",
                               "Olaparib", "Talazoparib"),
  "Protein\nstability"     = c("BIIB021", "ML-792"),
  "Epigenetic\nregulation" = c("Entinostat", "Vorinostat", "GSK126",
                               "Tazemetostat", "JQ1", "Molibresib"))

row_order <- unlist(ROW_GROUPS, use.names = FALSE)

# Concentration range assayed per compound, from Extended Data Table S7.
# Used only to audit the plotted values; it censors nothing.
DOSE_RANGE <- tibble::tribble(
  ~Compound,       ~top,  ~bottom,
  "Abemaciclib",   10,    0.03,
  "Adavosertib",   10,    0.03,
  "Avapritinib",   10,    0.03,
  "BAY-320",       31.6,  0.1,
  "BIIB021",       1,     0.003,
  "Berzosertib",   10,    0.03,
  "CHIR-124",      3.16,  0.01,
  "Ceralasertib",  10,    0.03,
  "Cobimetinib",   1,     0.003,
  "Dinaciclib",    0.316, 0.001,
  "Entinostat",    31.6,  0.1,
  "GSK1070916",    3.16,  0.01,
  "GSK126",        31.6,  0.1,
  "GSK461364",     1,     0.003,
  "JQ1",           10,    0.03,
  "KU-60019",      31.6,  0.1,
  "ML-792",        3.16,  0.01,
  "MLN8054",       10,    0.03,
  "Molibresib",    3.16,  0.01,
  "Nutlin-3a",     10,    0.03,
  "Olaparib",      3.16,  0.01,
  "Palbociclib",   3.16,  0.01,
  "Ripretinib",    3.16,  0.01,
  "SCH900776",     10,    0.03,
  "Talazoparib",   0.316, 0.001,
  "Tazemetostat",  10,    0.03,
  "Trametinib",    3.16,  0.01,
  "Volasertib",    3.16,  0.01,
  "Vorinostat",    3.16,  0.01,
  "XAV939",        10,    0.03)

DILUTION <- 3.162   # half-log steps, as used in the assay

## ------------------------------------------------------------------ input
# read.csv2, not read.csv: this file is semicolon-separated with a comma as the
# decimal mark. read.csv would see one column per row and stop with
# "more columns than column names".
ic50 <- read.csv2(file.path(IN, "IC50_acral_drugs_annotations_FA.csv"),
                  header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE, row.names = 1)

# Seven compounds carry a trailing space. Indexing a data frame by a name that
# does not match returns a row of NA rather than an error.
rownames(ic50) <- trimws(rownames(ic50))
colnames(ic50) <- trimws(colnames(ic50))

ic50 <- ic50[, -1]                       # drop the functional-class column
ic50 <- ic50 |>
  mutate(across(where(is.character), ~ as.numeric(gsub(",", ".", .x))))

names(ic50) <- sub("\\s*\\(IC50\\)\\s*$", "", names(ic50))

missing <- setdiff(CELL_LINES, colnames(ic50))
if (length(missing))
  stop("cell line(s) not found: ", paste(missing, collapse = ", "),
       "\ncolumns present: ", paste(colnames(ic50), collapse = ", "))
ic50 <- ic50[, CELL_LINES]

# Checked both ways: a compound absent from ROW_GROUPS would be dropped from the
# figure with no warning, which is how Tazemetostat was lost previously.
missing <- setdiff(rownames(ic50), row_order)
if (length(missing)) stop("compound in the data but not in ROW_GROUPS: ",
                          paste(missing, collapse = ", "))
missing <- setdiff(row_order, rownames(ic50))
if (length(missing)) stop("compound in ROW_GROUPS but not in the data: ",
                          paste(missing, collapse = ", "))

raw <- as.matrix(ic50)[row_order, ]
cat(sprintf("compounds %d | cell lines %d\n", nrow(raw), ncol(raw)))

## ----------------------------------------------------------------- transform
# No censoring: every fitted IC50 is plotted, and NA comes only from the input.
mat <- log10(raw)

cat(sprintf("log10 range plotted: %.2f to %.2f | white cells: %d of %d\n",
            min(mat, na.rm = TRUE), max(mat, na.rm = TRUE),
            sum(is.na(mat)), length(mat)))

all_na <- rownames(mat)[rowSums(!is.na(mat)) == 0]
if (length(all_na))
  message("compound(s) plotted as an all-white row: ",
          paste(all_na, collapse = ", "))

## --------------------------------------------------------------------- audit
# A compound missing from DOSE_RANGE would silently drop out of the join and
# never be audited, so check the two lists match before joining.
stopifnot(setequal(DOSE_RANGE$Compound, rownames(raw)))

outside <- as.data.frame(as.table(raw)) |>
  setNames(c("Compound", "Cell_line", "IC50")) |>
  filter(!is.na(IC50)) |>
  left_join(DOSE_RANGE, by = "Compound") |>
  filter(IC50 > top | IC50 < bottom) |>
  mutate(side           = if_else(IC50 > top, "above top dose", "below bottom dose"),
         fold_outside   = if_else(IC50 > top, IC50 / top, bottom / IC50),
         dilution_steps = log(fold_outside) / log(DILUTION)) |>
  arrange(desc(dilution_steps))

write_tsv(outside, file.path(OUT, "drugscreen_outside_dose_range.tsv"))

cat(sprintf("\nplotted values outside the range assayed: %d of %d\n",
            nrow(outside), sum(!is.na(raw))))
print(as.data.frame(outside), row.names = FALSE, digits = 3)

far <- filter(outside, dilution_steps >= 1)
if (nrow(far))
  message("\nmore than one dilution step beyond the last dose tested, ",
          "name these in the figure legend:\n  ",
          paste(unique(far$Compound), collapse = ", "))

## ------------------------------------------------------------------- heatmap
# Scale fixed at +/-4, as in the published panel, not to the range of the data.
LIM     <- 4
CELL_MM <- 7      # side of each cell in mm

# Two constraints at once.
# White is reserved for cells with no IC50 (na_col), so it must not appear in
# the ramp: there is no stop at zero, and the midpoint interpolates straight
# from lightcoral to lightsteelblue.
# The legend keeps the published +/-4, but the data only span -3 to 1.5, so the
# stops saturate inside that interval and the ends are flat. Anchoring the
# strongest colours at +/-4 instead would leave every cell washed out.
col_fun <- colorRamp2(c(-LIM, -3, -2, -1, -0.3, 0.3, 1, 1.5, LIM),
                      c("darkred", "darkred", "firebrick3", "indianred",
                        "lightcoral", "lightsteelblue", "steelblue",
                        "royalblue4", "royalblue4"), space = "LAB")

near_white <- vapply(seq(-LIM, LIM, 0.05),
                     \(v) all(col2rgb(col_fun(v)) > 247), logical(1))
if (any(near_white))
  stop("the colour ramp passes through white, which is the 'no IC50' category")

row_split <- factor(rep(names(ROW_GROUPS), lengths(ROW_GROUPS)),
                    levels = names(ROW_GROUPS))

ht <- Heatmap(mat,
              name             = "IC50\n(log10, uM)",
              na_col           = "white",
              col              = col_fun,
              # Grey, not white: a white border would make the "no IC50" cells
              # disappear into the background instead of reading as empty.
              rect_gp          = gpar(col = "gray70", lwd = 0.4),
              cluster_rows     = FALSE,
              cluster_columns  = FALSE,
              row_split        = row_split,
              row_title_rot    = 0,
              row_gap          = unit(1.5, "mm"),
              row_names_side   = "left",
              column_names_rot = 45,
              # Fixed cell size, so the panel keeps its proportions whatever the
              # page size. Change CELL_MM, not the pdf width and height.
              width            = unit(ncol(mat) * CELL_MM, "mm"),
              height           = unit(nrow(mat) * CELL_MM, "mm"),
              heatmap_legend_param = list(at = seq(-LIM, LIM, LIM/2)))

pdf(file.path(OUT, "drugscreen_heatmap.pdf"), width = 6, height = 10)
# Extra bottom padding, otherwise the rotated cell-line labels are clipped.
draw(ht, padding = unit(c(12, 2, 2, 2), "mm"))
dev.off()

as.data.frame(mat) |> rownames_to_column("Compound") |>
  write_tsv(file.path(OUT, "drugscreen_IC50_log10.tsv"), na = "")

message("saved: ", file.path(OUT, "drugscreen_heatmap.pdf"))