# Author: Annie Squiavinato
# 2026-09-04, revised 2026-10
#_______________________________________________________________________________
# Extended Data Figure 21: expression change of the CRISPR-essential genes that
# are differentially spliced in both AM cell lines, coloured by the predicted
# protein consequence of the splicing event.
#
# Genes with more than one event in a line are drawn as a bar cut corner to
# corner into two triangles, one per consequence category. This replaces
# colouring the gene by its largest-|dPSI| event, which was unstable: for SF1
# the winning event led by 0.1 dPSI and for UFD1L by 0.8, enough to flip the
# colour.
#
# Two figures are written, differing only in which bars are split:
#   _by_event     a bar is split when the gene has > 1 event   <- Fig. 21
#   _by_category  a bar is split when the events fall into > 1 category
#
# INPUTS  (under IN_DIR)
#   Functional_impact_diff_exons_introns_UTR_NonCoding_filtered.xlsx
#       the UTR_NonCoding variant, as used for Extended Data Table 11
#   Extended_Data_Table_9_shared_vulnerabilities.tsv   823 genes, one column
#   AM202607_sgSF3B4_vs_sgCtrol_edgeR.xlsx
#
# OUTPUTS (under OUT_DIR)
#   essential_genes_split_bars_by_event.svg
#   essential_genes_split_bars_by_category.svg
#   essential_genes_split_bars.csv
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from the project directory that holds renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggsci)
library(grid)

IN_DIR  <- "data"
OUT_DIR <- "results"

if (!dir.exists(IN_DIR))
  stop("'", IN_DIR, "' not found. Run this script from the project directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

IMPACT_PATH <- file.path(IN_DIR,
  "Functional_impact_diff_exons_introns_UTR_NonCoding_filtered.xlsx")
CRISPR_PATH <- file.path(IN_DIR, "Extended_Data_Table_9_shared_vulnerabilities.tsv")
EDGER_PATH  <- file.path(IN_DIR, "AM202607_sgSF3B4_vs_sgCtrol_edgeR.xlsx")

for (f in c(IMPACT_PATH, CRISPR_PATH, EDGER_PATH))
  if (!file.exists(f)) stop("missing input: ", f)

LINES <- c("AM007", "AM016")

## ------------------------------------------------------------------ input
# Each sheet carries BOTH dPSI columns, so the one matching the sheet's own
# line must be chosen explicitly.
load_sheet <- function(line) {
  wanted <- paste0("dPSI_", line)
  df <- read_excel(IMPACT_PATH, sheet = paste0("dPSI_", line))
  if (!wanted %in% names(df))
    stop("sheet dPSI_", line, " lacks '", wanted, "'. Columns: ",
         paste(names(df), collapse = ", "))
  df |> mutate(dPSI = .data[[wanted]], cell_line = line)
}
imp <- bind_rows(lapply(LINES, load_sheet))

# Both label sets are mapped so the script runs with either impact file: the
# plain "_filtered" one, which merges UTR and non-coding under a single label,
# or the "_UTR_NonCoding_filtered" variant used for Extended Data Table 11.
# They collapse to the same four categories here, so the figure is identical.
impact_recode <- c("ORF disruption" = "ORF disruption",
                   "Alt. isoforms"  = "Alternative protein isoforms",
                   "UTR, NonCoding" = "Noncoding",
                   "NonCoding"      = "Noncoding",
                   "5 UTR"          = "Noncoding",
                   "3 UTR"          = "Noncoding",
                   "NA, Uncertain"  = "NA / Uncertain")
imp <- imp |> mutate(category = unname(impact_recode[IMPACT]))
if (any(is.na(imp$category)))
  stop("unmapped IMPACT label(s): ",
       paste(unique(imp$IMPACT[is.na(imp$category)]), collapse = ", "))

crispr <- read_tsv(CRISPR_PATH, show_col_types = FALSE)
genes_crispr <- crispr[[1]] |> as.character() |> trimws() |> unique()
genes_crispr <- genes_crispr[nzchar(genes_crispr) & !is.na(genes_crispr)]
if (length(genes_crispr) < 800) warning("expected 823 genes, got ",
                                        length(genes_crispr))

# Gene set: CRISPR-essential and differentially spliced in both lines. This is
# gene-level concordance -- the same event is not required.
g007 <- unique(na.omit(imp$GENE[imp$cell_line == "AM007"]))
g016 <- unique(na.omit(imp$GENE[imp$cell_line == "AM016"]))
genes_overlap <- intersect(genes_crispr, intersect(g007, g016))
cat(sprintf("genes in figure: %d\n", length(genes_overlap)))

imp <- filter(imp, GENE %in% genes_overlap)

edger <- bind_rows(lapply(LINES, function(l)
  read_excel(EDGER_PATH, sheet = paste0(l, "_sgSF3B4_vs_", l, "_sgCtrol")) |>
    transmute(GENE = NAME, logFC, FDR, cell_line = l))) |>
  filter(GENE %in% genes_overlap) |>
  distinct(GENE, cell_line, .keep_all = TRUE)

## ------------------------------------------------------- categories, palette
SEV <- c("ORF disruption", "Alternative protein isoforms", "Noncoding",
         "NA / Uncertain")              # most to least severe
npg <- pal_npg("nrc")(10)
PAL <- c("ORF disruption"               = npg[8],
         "Alternative protein isoforms" = npg[3],
         "Noncoding"                    = npg[9],
         "NA / Uncertain"               = "grey70")

cat_tbl <- imp |>
  distinct(GENE, cell_line, category) |>
  mutate(rank = match(category, SEV)) |>
  arrange(GENE, cell_line, rank) |>
  group_by(GENE, cell_line) |>
  summarise(n_cat = n(), cat_1 = category[1],
            cat_2 = if (n() >= 2) category[2] else NA_character_,
            all_cats = paste(category, collapse = " + "), .groups = "drop")

# cat_2 falls back to cat_1 so a bar split by the event rule but holding a
# single category is drawn in one colour instead of producing an NA fill.
plot_dat <- cat_tbl |>
  left_join(count(imp, GENE, cell_line, name = "n_events"),
            by = c("GENE", "cell_line")) |>
  inner_join(edger, by = c("GENE", "cell_line")) |>
  mutate(cat_2 = coalesce(cat_2, cat_1),
         split_by_category = n_cat >= 2,
         split_by_event    = n_events >= 2)

# A diagonal shows exactly two categories, so a third would be dropped without
# any sign of it in the figure.
over_two <- filter(plot_dat, n_cat > 2)
if (nrow(over_two))
  stop("the diagonal holds two categories, but these have more:\n  ",
       paste(sprintf("%s (%s): %s", over_two$GENE, over_two$cell_line,
                     over_two$all_cats), collapse = "\n  "))

write_csv(plot_dat, file.path(OUT_DIR, "essential_genes_split_bars.csv"))
cat("\nsplit by category:", sum(plot_dat$split_by_category),
    "| split by event:", sum(plot_dat$split_by_event), "\n")

## ------------------------------------------------------------------ drawing
# Half-height of a bar on the discrete y axis: genes are 1 unit apart, so
# 0.375 fills 75% of the slot and matches geom_col(width = 0.75).
BAR_H   <- 0.375
DIAG_LW <- 0.5   # raise to 0.8-1.2 if the diagonal must read at print size

# ggplot2 has no glyph for a diagonally split key, so one is supplied here.
# Coordinates are npc units and follow the same layout as mk_tri() below.
draw_key_split <- function(data, params, size) {
  grobTree(
    polygonGrob(x = c(0, 1, 0), y = c(0, 0, 1),
                gp = gpar(fill = PAL[["ORF disruption"]], col = NA)),
    polygonGrob(x = c(0, 1, 1), y = c(1, 0, 1),
                gp = gpar(fill = PAL[["Alternative protein isoforms"]], col = NA)),
    segmentsGrob(x0 = 0, y0 = 1, x1 = 1, y1 = 0, gp = gpar(col = "white", lwd = 1.4)))
}

make_plot <- function(dat, gene_levels, split_col, out_svg, key_label, height_in) {
  dat <- dat |>
    filter(GENE %in% gene_levels) |>
    mutate(GENE = factor(GENE, levels = gene_levels),
           y = as.numeric(GENE), .split = .data[[split_col]])

  singles <- filter(dat, !.split)
  splits  <- filter(dat, .split)

  # Diagonal from the top-left corner (0, y+h) to the bottom-right (logFC, y-h):
  #   lo = (0, y-h), (logFC, y-h), (0, y+h)      -> cat_1, the more severe
  #   up = (0, y+h), (logFC, y-h), (logFC, y+h)  -> cat_2
  # The shared edge is the same pair of points as the geom_segment below, so
  # the drawn separation and the drawn line coincide.
  mk_tri <- function(df, part) {
    if (!nrow(df)) return(tibble(poly_id = character(), cell_line = character(),
                                 category = character(), x = numeric(), yy = numeric()))
    idx  <- rep(seq_len(nrow(df)), each = 3)
    vert <- rep(1:3, times = nrow(df))
    d <- df[idx, ]
    tibble(
      poly_id   = paste0(d$GENE, "_", d$cell_line, "_", part),
      cell_line = d$cell_line,
      category  = if (part == "lo") d$cat_1 else d$cat_2,
      x  = if (part == "lo") ifelse(vert == 2, d$logFC, 0)
           else              ifelse(vert == 1, 0, d$logFC),
      yy = if (part == "lo") ifelse(vert == 3, d$y + BAR_H, d$y - BAR_H)
           else              ifelse(vert == 2, d$y - BAR_H, d$y + BAR_H))
  }
  tris <- bind_rows(mk_tri(splits, "lo"), mk_tri(splits, "up"))

  p <- ggplot() +
    geom_rect(data = singles,
              aes(xmin = pmin(0, logFC), xmax = pmax(0, logFC),
                  ymin = y - BAR_H, ymax = y + BAR_H, fill = cat_1),
              colour = "white", linewidth = 0.1) +
    geom_polygon(data = tris,
                 aes(x = x, y = yy, group = poly_id, fill = category),
                 colour = "white", linewidth = 0.1) +
    # The diagonal is its own segment rather than a thicker polygon stroke,
    # which would also thicken the outer edges of the bar.
    geom_segment(data = splits,
                 aes(x = 0, xend = logFC, y = y + BAR_H, yend = y - BAR_H),
                 colour = "white", linewidth = DIAG_LW, lineend = "butt") +
    # Zero-area invisible layer whose only purpose is the split-bar legend
    # entry. show.legend must name the aesthetics: a bare TRUE would stamp
    # this glyph over the four Consequence swatches as well.
    geom_polygon(data = tibble(x = rep(0, 3), yy = rep(1, 3), key = key_label),
                 aes(x = x, y = yy, colour = key), fill = NA, alpha = 0,
                 key_glyph = draw_key_split, inherit.aes = FALSE,
                 show.legend = c(colour = TRUE, fill = FALSE)) +
    scale_colour_manual(values = setNames("white", key_label), name = NULL) +
    guides(fill = guide_legend(order = 1),
           colour = guide_legend(order = 2, keyheight = unit(0.55, "cm"))) +
    geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.4) +
    facet_wrap(~cell_line, ncol = 2) +
    scale_fill_manual(values = PAL, name = "Consequence", drop = FALSE) +
    scale_y_continuous(breaks = seq_along(gene_levels), labels = gene_levels,
                       expand = expansion(add = 0.7)) +
    scale_x_continuous(expand = expansion(mult = 0.04)) +
    labs(x = "logFC (SF3B4-KD vs control)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 6.5, colour = "grey15"),
          axis.line.x = element_line(colour = "grey40", linewidth = 0.4),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          legend.position = "right",
          legend.key.size = unit(0.45, "cm"),
          legend.title = element_text(face = "bold", size = 9.5),
          plot.margin = margin(10, 14, 10, 10))

  ggsave(out_svg, p, width = 10, height = height_in, dpi = 200, limitsize = FALSE)
  cat("saved:", out_svg, sprintf("(%d genes)\n", length(gene_levels)))
  p
}

# Shared gene order: mean logFC across the two lines, lowest at the bottom.
gene_order <- plot_dat |>
  group_by(GENE) |> summarise(ord = mean(logFC), .groups = "drop") |>
  arrange(ord) |> pull(GENE)
h <- max(6, 1.4 + length(gene_order) * 0.13)

make_plot(plot_dat, gene_order, "split_by_event",
          file.path(OUT_DIR, "essential_genes_split_bars_by_event.svg"),
          "More than one splicing event", h)

make_plot(plot_dat, gene_order, "split_by_category",
          file.path(OUT_DIR, "essential_genes_split_bars_by_category.svg"),
          "Two consequence categories", h)

# For the figure legend: in a split bar the lower-left triangle is the more
# severe category and the upper-right the other, in the order ORF disruption >
# alternative protein isoforms > noncoding > NA/uncertain. The equal division
# does not reflect how many events fall in each category.
cat("\ngenes with more than one category:\n")
print(plot_dat |> filter(split_by_category) |>
        select(cell_line, GENE, n_events, n_cat, all_cats) |> arrange(cell_line, GENE),
      n = 30)
