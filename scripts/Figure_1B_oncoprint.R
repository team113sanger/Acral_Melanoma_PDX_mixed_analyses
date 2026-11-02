### Author: Annie Squiavinato and Antonio
### 2025-04-25, revised 2026-10
#_______________________________________________________________________________
# Figure 1a: data and models available per acral melanoma case.
# Body = presence/absence matrix; top = anatomic site and sample type.
#
# INPUTS  (under IN)
#   figure_1a_table_23_04_25.tsv   clinical annotation, one row per case,
#                                  columns Anatomic_site and Sample_type
#   data_fig_1a_sept_26.xlsx       Yes/NA matrix, cases in rows, data types
#                                  in columns
#
# OUTPUTS (under OUT)
#   figure_1b_oncoprint.pdf
#
# Post-processing in Affinity for the submitted figure
# Package versions are managed with renv: run renv::restore() once to rebuild
# the environment recorded in renv.lock. Run this script from INCA_res/, the
# directory that holds renv.lock, data/, scripts/ and results/.
#_______________________________________________________________________________

library(ComplexHeatmap)
library(grid)
library(readxl)

IN  <- "data"
OUT <- "results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

F_CLIN <- file.path(IN, "figure_1a_table_23_04_25.tsv")   # clinical annotation
F_DATA <- file.path(IN, "data_fig_1a_241026.xlsx")       # Yes/NA matrix
for (f in c(F_CLIN, F_DATA)) if (!file.exists(f)) stop("missing input: ", f)

DROP_EMPTY <- FALSE   # TRUE hides cases with no data at all

# ------------------------------------------------------------------ data
# The two files play opposite roles and are easy to swap: the xlsx is the
# body of the oncoPrint, the tsv is the annotation on top.
clin_mat <- read.table(F_CLIN, header = TRUE, sep = "\t", row.names = 1,
                       stringsAsFactors = FALSE)

x        <- read_excel(F_DATA)           # read.table cannot open .xlsx
anno_mat <- as.data.frame(x[, -1])
rownames(anno_mat) <- x[[1]]

final_mat <- t(as.matrix(anno_mat))      # variables in rows, samples in columns

if (DROP_EMPTY) {
  n0 <- ncol(final_mat)
  final_mat <- final_mat[, colSums(!is.na(final_mat)) > 0, drop = FALSE]
  message("cases with no data dropped: ", n0 - ncol(final_mat))
}

# --------------------------------------------------------- display names
# Row labels as they should read in the figure. Matching on the file's own
# names, so a renamed column in the source fails here instead of silently
# dropping a row.
ROW_LABEL <- c(
  Patient_WES               = "Patient WES",
  Patient_RNAseq            = "Patient RNAseq",
  PDX_WES                   = "AM-PDX WES",
  PDX_RNAseq                = "AM-PDX RNAseq",
  Ancestry                  = "Ancestry",
  PDX_model_available       = "AM-PDX model",
  `Mestastasis_in_mice(Xn)` = "Metastasis in mice",
  Cell_line_available       = "AM cell line"
)
no_label <- setdiff(rownames(final_mat), names(ROW_LABEL))
if (length(no_label)) stop("rows with no display label: ",
                           paste(no_label, collapse = ", "))

final_mat <- final_mat[names(ROW_LABEL), , drop = FALSE]   # fixes row order
rownames(final_mat) <- unname(ROW_LABEL)

# ------------------------------------------------------- align annotation
shared <- intersect(colnames(final_mat), rownames(clin_mat))
message("matrix: ", ncol(final_mat), " | annotation: ", nrow(clin_mat),
        " | in both: ", length(shared))
no_annotation <- setdiff(colnames(final_mat), rownames(clin_mat))
if (length(no_annotation))
  message("no clinical annotation: ", paste(no_annotation, collapse = ", "))
stopifnot(length(shared) > 0)

final_mat <- final_mat[, shared, drop = FALSE]

# HeatmapAnnotation matches by POSITION, not by name: the i-th row of the
# data frame is drawn over the i-th column. Reordering here is what stops
# the annotation from being silently shifted, which would look fine and be
# wrong.
ha_df <- clin_mat[shared, , drop = FALSE]
stopifnot(identical(rownames(ha_df), colnames(final_mat)))
message("annotation columns: ", paste(names(ha_df), collapse = ", "))

# Lower-case level labels for the legend, recoded in the data so the
# palette, the annotation and the legend cannot drift apart.
LVL <- c(Plantar               = "plantar",
         Subungual             = "subungual",
         In_transit_Metastasis = "in transit",
         Lymph_Node            = "lymph node",
         Primary               = "primary",
         Recurrence            = "recurrence",
         Regional_Metastasis   = "regional metastasis")

unmapped <- setdiff(unlist(lapply(ha_df, \(v) na.omit(unique(v)))), names(LVL))
if (length(unmapped)) stop("levels with no label: ",
                           paste(unmapped, collapse = ", "))

ha_df[] <- lapply(ha_df, \(v) unname(LVL[v]))
names(ha_df) <- c(Anatomic_site = "Anatomic site",
                  Sample_type   = "Sample type")[names(ha_df)]

# --------------------------------------------------------------- colours
col       <- c(Yes = "#00a087")
anato_col <- c(plantar      = "#e64b35", subungual    = "#f9d6d2",
               `in transit` = "#9f2414", `lymph node` = "#f1998e")
type_col  <- c(primary = "#4dbbd5", recurrence = "#98d8e7",
               `regional metastasis` = "#278ea5")

ANN_COL <- list(`Anatomic site` = anato_col, `Sample type` = type_col)

# A palette whose name does not match a column is ignored without warning
# and ComplexHeatmap invents colours instead.
if (!all(names(ANN_COL) %in% names(ha_df)))
  stop("palette names do not match the annotation columns: ",
       paste(names(ha_df), collapse = ", "))

# --------------------------------------------------------------- drawing
alter_fun <- list(
  background = function(x, y, w, h)
    grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "white")),
  Yes = function(x, y, w, h)
    grid.rect(x, y, w * 0.92, h * 0.92, gp = gpar(fill = col["Yes"], col = NA))
)

# Every legend is built by hand from the same SW. Letting oncoPrint and
# HeatmapAnnotation build their own produced different swatch sizes, and
# grid_height/grid_width did not fix it because 4 mm is already the
# default -- the two legends were not taking the same code path. Building
# all three here makes the size identical by construction, and also fixes
# the order of the labels.
SW <- unit(4, "mm")

# Cell height, shared by the annotation tracks and the body rows. The two
# are controlled separately -- simple_anno_size for the annotation, total
# heatmap height divided by the number of rows for the body -- so without
# a common constant they never match. Column width is already shared,
# since both use the same columns.
CELL <- unit(6, "mm")

ha <- HeatmapAnnotation(
  df     = ha_df,
  which  = "col",
  col    = ANN_COL,
  na_col = "white",
  simple_anno_size     = CELL,
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 11),
  gp          = gpar(col = "white"),
  show_legend = FALSE            # legends built below instead
)

p1 <- oncoPrint(
  final_mat,
  alter_fun               = alter_fun,
  col                     = col,
  top_annotation          = ha,
  right_annotation        = NULL,
  show_pct                = FALSE,
  show_column_names       = TRUE,
  column_names_gp         = gpar(fontsize = 7),
  row_names_side          = "left",
  row_names_gp            = gpar(fontsize = 11),
  row_order               = rownames(final_mat),   # keep the order above
  alter_fun_is_vectorized = FALSE,
  show_heatmap_legend     = FALSE,                 # replaced by lgd_data
  height                  = CELL * nrow(final_mat) # same cell height as ha
)

# All three legends from the same call and the same SW.
mk <- function(pal, title, nrow)
  Legend(labels = names(pal), title = title,
         legend_gp      = gpar(fill = unname(pal)),
         nrow           = nrow,
         grid_height    = SW,
         grid_width     = SW,
         title_position = "topleft",
         title_gp       = gpar(fontsize = 11, fontface = "bold"),
         labels_gp      = gpar(fontsize = 10))

# Label order set here, not by the palette definition: two columns of two
# for the site, one column of three for the type.
lgd_site <- mk(anato_col[c("in transit", "lymph node", "plantar", "subungual")],
               "Anatomic site", 2)
lgd_type <- mk(type_col[c("primary", "recurrence", "regional metastasis")],
               "Sample type", 3)
lgd_data <- mk(setNames(col["Yes"], "yes"), "Data and models", 1)

FIG <- file.path(OUT, "figure_1a_oncoprint.pdf")
pdf(FIG, width = 14, height = 6)
draw(p1,
     annotation_legend_side = "bottom",
     annotation_legend_list = list(lgd_site, lgd_type, lgd_data))
dev.off()

message("figure: ", FIG)
