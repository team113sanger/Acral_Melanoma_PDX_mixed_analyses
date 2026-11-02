# Author: Annie Squiavinato
# 2026-10
#_______________________________________________________________________________
# Figure S5C, clinical annotation tracks: tumour type, clinical stage and
# anatomic site, one tile per sample.
#
# The sample list and order are the same as Figure_S5C_signatures.R, so the
# three tracks stack under the substitution-type, SBS and CNV signature panels.
# They are repeated here rather than imported so the script runs on its own.
#
# INPUTS  (under IN)
#   clinical_info_one_sp_april_26.csv
#       columns ID_Sanger, ID_lab, Mutational_status, Anatomic_site,
#       Tumor_type, Clinical_Stage. ID_lab is already in AM0xx form.
#
# OUTPUTS (under OUT)
#   clinical_annotation_tracks.pdf
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from INCA_res/, the directory that holds renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(patchwork)

IN  <- "data"
OUT <- "results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

F_CLIN <- file.path(IN, "clinical_info_one_sp_april_26.csv")
if (!file.exists(F_CLIN)) stop("missing input: ", F_CLIN)

## ------------------------------------------------------------- ID mapping
sanger_to_am <- c(
  PD53330c = "AM001", PD53331a = "AM002", PD53332d = "AM003a", PD53333f = "AM004",
  PD53334a = "AM005", PD53337a = "AM008a", PD53338a = "AM009", PD53339a = "AM010",
  PD53341a = "AM013", PD53342a = "AM014", PD53343h = "AM015a", PD53345a = "AM017",
  PD53346a = "AM018", PD53347c = "AM019", PD53348a = "AM020", PD53349d = "AM021a",
  PD53350d = "AM022a", PD53351a = "AM023", PD53352d = "AM025a", PD53353a = "AM026",
  PD53354a = "AM027", PD53355c = "AM029", PD53357c = "AM032a", PD53358a = "AM033",
  PD53359c = "AM034", PD53361a = "AM042", PD53362a = "AM043", PD53363a = "AM044",
  PD53364c = "AM045", PD53365a = "AM046", PD53366a = "AM047", PD53367a = "AM048",
  PD53368a = "AM049", PD53369a = "AM050", PD70963a = "AM051", PD70964a = "AM052a",
  PD70965a = "AM053", PD70967a = "AM058", PD70968a = "AM060a", PD70969a = "AM062",
  PD70970a = "AM065", PD70971a = "AM067", PD70972a = "AM068a", PD70973a = "AM069",
  PD70976a = "AM057a")

sample_order <- c(
  "PD53346a", "PD53369a", "PD70971a", "PD53343h", "PD53352d", "PD70970a",
  "PD53355c", "PD53333f", "PD53330c", "PD53337a", "PD53363a", "PD70965a",
  "PD70967a", "PD70972a", "PD70973a", "PD53350d", "PD53347c", "PD70976a",
  "PD53338a", "PD53368a", "PD53351a", "PD53357c", "PD70968a", "PD53348a",
  "PD70963a", "PD53359c", "PD70969a", "PD53362a", "PD53367a", "PD53364c",
  "PD53345a", "PD53349d", "PD53361a", "PD53334a", "PD70964a", "PD53332d",
  "PD53331a", "PD53339a", "PD53341a", "PD53342a", "PD53353a", "PD53354a",
  "PD53358a", "PD53365a", "PD53366a")

sample_order_am <- unname(sanger_to_am[sample_order])
stopifnot(!anyNA(sample_order_am), !anyDuplicated(sample_order))

## ------------------------------------------------------------------ input
meta <- read.csv(F_CLIN) |>
  mutate(Samples_AM = factor(ID_lab, levels = sample_order_am)) |>
  filter(!is.na(Samples_AM))

missing <- setdiff(sample_order_am, meta$Samples_AM)
if (length(missing))
  message("samples with no clinical row (left blank in the track): ",
          paste(missing, collapse = ", "))

## --------------------------------------------------------------- palettes
# Blue for tumour type, darkest for the commonest category.
PAL_TUMOUR <- c("Primary"               = "#2C3E70",
                "Recurrence"            = "#7B8FC9",
                "Lymph_node_metastasis" = "#D6DCEF")

# Clinical stage is ORDINAL (Ib < IIa < IIb < IIc < III < IV), so the scale
# runs light to dark; NA falls through to na.value.
PAL_STAGE <- c("Ib"  = "#D3F2E5", "IIa" = "#A8E6CF", "IIb" = "#6FCF9E",
               "IIc" = "#2E9E74", "III" = "#1B6B52", "IV"  = "#0D3D2E")

PAL_SITE <- c("Plantar"               = "#B23A2E",
              "Subungual"             = "#D98B76",
              "Lymph_node"            = "#F2D0C4",
              "in-transit metastasis" = "#FBEDE7")

# A level with no colour would be drawn as NA grey and read as missing data.
for (v in list(list("Tumor_type", PAL_TUMOUR), list("Clinical_Stage", PAL_STAGE),
               list("Anatomic_site", PAL_SITE))) {
  lev <- setdiff(unique(na.omit(meta[[v[[1]]]])), names(v[[2]]))
  lev <- lev[nzchar(lev) & lev != "NA"]
  if (length(lev)) stop(v[[1]], ": level(s) with no colour: ",
                        paste(lev, collapse = ", "))
}

## ---------------------------------------------------------------- tracks
# coord_fixed, not coord_cartesian: it forces one x unit to occupy the same
# physical size as one y unit, so each tile is square rather than stretched
# across the full figure width.
make_track <- function(df, var, palette, var_label, show_x_axis = FALSE) {
  ggplot(df, aes(x = Samples_AM, y = 1, fill = .data[[var]])) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = palette, name = var_label, na.value = "grey70") +
    coord_fixed(ratio = 1,
                xlim = c(0.5, length(levels(df$Samples_AM)) + 0.5), expand = FALSE) +
    labs(x = NULL, y = var_label) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y  = element_blank(),
          axis.ticks   = element_blank(),
          axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 9),
          axis.text.x  = if (show_x_axis)
            element_text(angle = 45, hjust = 1, size = 7) else element_blank(),
          panel.grid   = element_blank(),
          legend.key.size = unit(0.35, "cm"),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8))
}

p_tumour <- make_track(meta, "Tumor_type",     PAL_TUMOUR, "Tumour type")
p_stage  <- make_track(meta, "Clinical_Stage", PAL_STAGE,  "Clinical stage")
p_site   <- make_track(meta, "Anatomic_site",  PAL_SITE,   "Anatomic site",
                       show_x_axis = TRUE)

p_tracks <- p_tumour / p_stage / p_site + plot_layout(guides = "collect")

# Width and height are only the paper available; coord_fixed keeps the tiles
# square regardless. Reduce height if white space opens between the tracks.
ggsave(file.path(OUT, "clinical_annotation_tracks.pdf"), p_tracks,
       width = 9.16, height = 3)

message("samples drawn: ", nlevels(meta$Samples_AM),
        " | saved: ", file.path(OUT, "clinical_annotation_tracks.pdf"))

# To stack these under the signature panels, run Figure_S5C_signatures.R first
# in the same session and then:
#   plot1 / p_tumour / p_stage / p_site +
#     plot_layout(heights = c(6, 1, 1, 1), guides = "collect")
