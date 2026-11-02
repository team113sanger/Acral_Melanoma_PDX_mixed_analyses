# title: "Mutational and CNA signature plots"
# author: Annie Squiavinato
# 2025-05-16, revised 2026-10
#_______________________________________________________________________________
# Figure S5C: stacked bars of COSMIC SBS and copy-number signature activities
# per tumour, in absolute and relative contribution. Samples are labelled with
# the AM identifier.
#
# Low-cellularity samples are kept on the x axis as a grey "Low cellularity"
# bar instead of being dropped, so the sample order stays comparable across
# panels. This applies to the CNA panels only: cellularity was not an
# exclusion criterion for SBS.
#
# INPUTS  (under IN)
#   COSMIC_CNV48_Activities_sept_26.txt
#   COSMIC_SBS96_Activities_april_2026.txt
#
# OUTPUTS (under OUT)
#   CNA_signature_absolute_sept26.pdf
#   CNA_signature_relative_sept26.pdf
#   Mutational_signatures_absolute_ago26.pdf
#   Mutational_signatures_relative_ago26.pdf
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from INCA_res/, the directory that holds renv.lock.
#_______________________________________________________________________________

library(tidyverse)
library(ggsci)

IN  <- "data"
OUT <- "results"

if (!dir.exists(IN))
  stop("'data/' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

F_CNA <- file.path(IN, "COSMIC_CNV48_Activities_sept_26.txt")
F_SBS <- file.path(IN, "COSMIC_SBS96_Activities_april_2026.txt")
for (f in c(F_CNA, F_SBS)) if (!file.exists(f)) stop("missing input: ", f)

# Height of the grey "Low cellularity" bar on the absolute panels. The y axis
# is the absolute number of alterations (~100-370 in real samples), so a value
# of 1 would be invisible and the sample would appear to be missing.
LOW_CELL_BAR_HEIGHT <- 100

## ------------------------------------------------------------- ID mapping
# Sanger ID, as it appears in the data files -> AM ID, as shown on the x axis.
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

## ---------------------------------------------------------------- helpers
# Adds one row per sample listed in sample_order but absent from the file.
# Those rows get 0 for every real signature and LOW_CELL_BAR_HEIGHT for the
# new "Low cellularity" column, so the sample keeps its place on the x axis.
add_missing_samples <- function(df, sample_order) {
  df$Samples <- as.character(df$Samples)
  missing <- setdiff(sample_order, df$Samples)

  # Only create the column if absent, so a file that already carries real
  # values in it is not overwritten with 0.
  if (!"Low cellularity" %in% colnames(df)) df[["Low cellularity"]] <- 0

  if (length(missing)) {
    cat("samples with no data (shown as 'Low cellularity'):\n  ",
        paste(missing, collapse = ", "), "\n")
    new_rows <- as.data.frame(matrix(0, nrow = length(missing), ncol = ncol(df)))
    colnames(new_rows) <- colnames(df)
    new_rows$Samples <- missing
    new_rows[["Low cellularity"]] <- LOW_CELL_BAR_HEIGHT
    df <- rbind(df, new_rows)
  } else {
    cat("no sample missing from this file -- all", length(sample_order),
        "have data\n")
  }

  extras <- setdiff(df$Samples, sample_order)
  if (length(extras))
    cat("rows in the file but not in the sample list (left out of the plot):\n  ",
        paste(extras, collapse = ", "), "\n")
  df
}

# For samples that DO have data but are known to be low-cellularity: zero every
# signature and mark "Low cellularity", so the sample stays on the axis with the
# real values discarded rather than disappearing.
force_low_cellularity <- function(df, samples_to_force) {
  if (!"Low cellularity" %in% colnames(df)) df[["Low cellularity"]] <- 0
  hit <- df$Samples %in% samples_to_force
  if (any(hit)) {
    signal_cols <- setdiff(colnames(df), c("Samples", "Samples_AM", "Low cellularity"))
    df[hit, signal_cols] <- 0
    df[hit, "Low cellularity"] <- LOW_CELL_BAR_HEIGHT
    cat("CNA: forced to 'Low cellularity' (real values discarded, kept on the plot):\n  ",
        paste(df$Samples[hit], collapse = ", "), "\n")
  }
  df
}

## ============================================================= CNA panels
# Low cellularity -> kept, but the whole bar becomes grey.
# Not in this cohort -> removed from the figure entirely.
cna_low_cellularity <- c("PD53331a", "PD53337a", "PD53339a", "PD53346a",
                         "PD53353a", "PD53363a", "PD53366a", "PD70976a")
cna_not_in_cohort   <- c("PD70974a", "PD70975a")

sample_order_cna    <- setdiff(sample_order, cna_not_in_cohort)
sample_order_cna_am <- unname(sanger_to_am[sample_order_cna])

acts_ed <- read.delim(F_CNA, header = TRUE)
acts_ed <- na.omit(acts_ed)

n_before <- nrow(acts_ed)
acts_ed  <- filter(acts_ed, !Samples %in% cna_not_in_cohort)
cat("CNA: removed", n_before - nrow(acts_ed), "row(s) not in the cohort\n")

acts_ed <- add_missing_samples(acts_ed, sample_order_cna)
acts_ed <- force_low_cellularity(acts_ed, cna_low_cellularity)

acts_ed$Samples_AM <- factor(unname(sanger_to_am[acts_ed$Samples]),
                             levels = sample_order_cna_am)

# Rows outside the sample list became NA above; without this the plot gains a
# phantom "NA" bar on the x axis.
n_before <- nrow(acts_ed)
acts_ed  <- filter(acts_ed, !is.na(Samples_AM))
cat("CNA: removed", n_before - nrow(acts_ed), "extra row(s) (Samples_AM = NA)\n")

acts_tidy <- acts_ed |>
  pivot_longer(cols = !c(Samples, Samples_AM),
               names_to = "Signature", values_to = "Alterations")

npg_colors   <- pal_npg("nrc")(10)
extra_colors <- c("#FF7F00", "#6A3D9A", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
my_palette   <- c(npg_colors, extra_colors)

cna_pastel_colors <- c(
  "#C7A6D3", "#D9C2A6", "#8FBF9F", "#D9A48F", "#9FB8D9", "#7FA6A6",
  "#C9A68F", "#E8C98A", "#E0C9A6", "#7FA6A6", "#D99FA9", "#A6C2D9",
  "#B8D98F", "#A6D9C2", "#D9A6A3", "#B0A6D9")

# Named palette, so "empty" and "Low cellularity" always get grey and never
# take a colour belonging to a real signature.
signature_levels <- unique(acts_tidy$Signature)
cna_palette <- setNames(cna_pastel_colors[seq_along(signature_levels)],
                        signature_levels)
cna_palette["empty"]            <- "grey50"
cna_palette["Low cellularity"]  <- "#D3D3D3"

n_samples <- length(sample_order_cna_am)

cna_plot <- function(position, ylab) {
  ggplot(acts_tidy, aes(x = Samples_AM, y = Alterations, fill = Signature)) +
    geom_bar(position = position, stat = "identity") +
    theme_classic() +
    labs(x = "Samples", y = ylab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = cna_palette) +
    coord_cartesian(xlim = c(0.5, n_samples + 0.5), expand = FALSE)
}

ggsave(file.path(OUT, "CNA_signature_absolute_sept26.pdf"),
       cna_plot("stack", "CNA Signatures Contribution (Absolute)"),
       device = "pdf", width = 9.89, height = 4.03)

ggsave(file.path(OUT, "CNA_signature_relative_sept26.pdf"),
       cna_plot("fill", "CNA Signatures Contribution (Relative)"),
       device = "pdf", width = 9.89, height = 4.03)

## ===================================================== mutational panels
# Full sample_order, without the CNA exclusions: cellularity was not an
# exclusion criterion for SBS.
mut_sig <- read.table(F_SBS, header = TRUE)
mut_sig <- na.omit(mut_sig)
mut_sig <- add_missing_samples(mut_sig, sample_order)

mut_sig$Samples_AM <- factor(unname(sanger_to_am[mut_sig$Samples]),
                             levels = sample_order_am)
mut_sig <- filter(mut_sig, !is.na(Samples_AM))

acts_tidy_mut <- mut_sig |>
  pivot_longer(cols = !c(Samples, Samples_AM),
               names_to = "Signature", values_to = "Mutations")

# Fixed colour per signature name, so a given SBS keeps its colour whatever
# set of samples or signatures is present.
mut_sig_colors <- c(
  SBS1   = "#E8A0A8",  # dusty rose
  SBS33  = "#B4998A",  # taupe
  SBS38  = "#D9EDE6",  # pale mint
  SBS40b = "#AAB4DD",  # periwinkle
  SBS5   = "#4FB8A8",  # teal
  SBS7a  = "#5A8FC2",  # steel blue
  SBS7b  = "#F6C3AE",  # pale peach
  SBS7d  = "#F2894E")  # orange

signature_levels_mut <- unique(acts_tidy_mut$Signature)
unlisted <- setdiff(signature_levels_mut,
                    c(names(mut_sig_colors), "empty", "Low cellularity"))
mut_palette <- c(mut_sig_colors,
                 setNames(my_palette[seq_along(unlisted)], unlisted))
mut_palette["empty"]           <- "grey50"
mut_palette["Low cellularity"] <- "#D3D3D3"

n_samples_mut <- length(sample_order_am)

mut_plot <- function(position, ylab, ylim = NULL) {
  ggplot(acts_tidy_mut, aes(x = Samples_AM, y = Mutations, fill = Signature)) +
    geom_bar(position = position, stat = "identity") +
    theme_classic() +
    labs(x = "Samples", y = ylab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = mut_palette) +
    coord_cartesian(xlim = c(0.5, n_samples_mut + 0.5), ylim = ylim,
                    expand = FALSE)
}

ggsave(file.path(OUT, "Mutational_signatures_absolute_ago26.pdf"),
       mut_plot("stack", " Mutational Signatures Contribution", c(0, 200)),
       device = "pdf", width = 9.89, height = 4.03)

ggsave(file.path(OUT, "Mutational_signatures_relative_ago26.pdf"),
       mut_plot("fill", " Relative Mutational Signatures Contribution"),
       device = "pdf", width = 9.89, height = 4.03)

message("done")
