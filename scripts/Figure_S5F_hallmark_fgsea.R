# Author: Annie Squiavinato
# 2026-08-19, revised 2026-10
#_______________________________________________________________________________
# Figure S5F: whole-transcriptome fgsea of the 50 MSigDB Hallmark gene sets after SF3B4
# knockout, in AM007 and AM016. The figure shows only the pathways significant
# (padj < 0.05) in BOTH lines AND agreeing on the sign of the NES.
#
# A pathway can reach padj < 0.05 in both lines with OPPOSITE NES, which a
# plain "significant in both" filter would keep; the direction check removes it.
#
# INPUTS  (under IN_DIR)
#   AM202607_sgSF3B4_vs_sgCtrol_edgeR.xlsx
#       sheets AM007_sgSF3B4_vs_AM007_sgCtrol and AM016_..., columns
#       NAME, logFC, PValue, FDR
#
# OUTPUTS (under OUT_DIR)
#   hallmark_fgsea_all_50.csv              50 gene sets x 2 cell lines
#   hallmark_fgsea_significant_both.svg
#
# fgseaMultilevel is stochastic, so set.seed() below fixes the p-values.
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from INCA_res/, the directory that holds renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggsci)
library(scales)
library(grid)
library(fgsea)
library(msigdbr)
library(svglite)


IN_DIR  <- "data"
OUT_DIR <- "results"

if (!dir.exists(IN_DIR))
  stop("'", IN_DIR, "' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EDGER_PATH <- file.path(IN_DIR, "AM202607_sgSF3B4_vs_sgCtrol_edgeR.xlsx")
if (!file.exists(EDGER_PATH)) stop("missing input: ", EDGER_PATH)

set.seed(1)

# npg positions 4 (navy) and 5 (salmon): a light/dark pair, so the two lines
# stay distinguishable in greyscale.
col_main <- setNames(pal_npg("nrc")(10)[c(4, 5)], c("AM007", "AM016"))

## ----------------------------------------------------- 1. edgeR, both lines
load_edger <- function(line)
  read_excel(EDGER_PATH, sheet = paste0(line, "_sgSF3B4_vs_", line, "_sgCtrol")) |>
    select(GENE = NAME, logFC, PValue, FDR)

edger_007 <- load_edger("AM007")
edger_016 <- load_edger("AM016")

## ------------------------------------------------- 2. the 50 Hallmark sets
# msigdbr renamed `category` to `collection` in v10; both are tried so the
# script runs with either version.
hallmark_df <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H"))

col_symbol <- intersect(c("gene_symbol", "human_gene_symbol"), names(hallmark_df))[1]
hallmark_list <- split(hallmark_df[[col_symbol]], hallmark_df$gs_name)
stopifnot(length(hallmark_list) == 50)

readable_names <- names(hallmark_list) |>
  gsub("^HALLMARK_", "", x = _) |>
  gsub("_", " ", x = _) |>
  tolower()

## ------------------------------------------------------------- 3. ranking
ranking_fgsea <- function(edger_tab) {
  stat <- sign(edger_tab$logFC) * -log10(edger_tab$PValue)
  names(stat) <- edger_tab$GENE
  stat <- stat[!is.na(stat) & !duplicated(names(stat)) & names(stat) != ""]
  sort(stat, decreasing = TRUE)
}

## ---------------------------------------------------------------- 4. fgsea
res_hallmark <- bind_rows(
  as_tibble(fgseaMultilevel(hallmark_list, ranking_fgsea(edger_007), eps = 0)) |>
    mutate(cell_line = "AM007"),
  as_tibble(fgseaMultilevel(hallmark_list, ranking_fgsea(edger_016), eps = 0)) |>
    mutate(cell_line = "AM016")) |>
  mutate(gs_name_official = pathway,
         pathway_label = readable_names[match(pathway, names(hallmark_list))]) |>
  select(pathway_label, gs_name_official, cell_line, NES, pval, padj, size)

write_csv(res_hallmark, file.path(OUT_DIR, "hallmark_fgsea_all_50.csv"))
cat(sprintf("\nsaved: hallmark_fgsea_all_50.csv (%d rows = 50 gene sets x 2 lines)\n",
            nrow(res_hallmark)))

## ----------------------------------------- 5. significance and direction
sig_007 <- res_hallmark |> filter(cell_line == "AM007", padj < 0.05) |> pull(pathway_label)
sig_016 <- res_hallmark |> filter(cell_line == "AM016", padj < 0.05) |> pull(pathway_label)

nes_wide <- res_hallmark |>
  select(pathway_label, cell_line, NES) |>
  pivot_wider(names_from = cell_line, values_from = NES)
same_direction <- nes_wide |> filter(sign(AM007) == sign(AM016)) |> pull(pathway_label)

sig_both_raw <- intersect(sig_007, sig_016)
sig_both     <- intersect(sig_both_raw, same_direction)

if (length(sig_both_raw) > length(sig_both))
  cat(sprintf("\n%d pathway(s) padj<0.05 in both lines but with opposite NES, excluded: %s\n",
              length(sig_both_raw) - length(sig_both),
              paste(setdiff(sig_both_raw, sig_both), collapse = ", ")))

cat(sprintf("\npadj<0.05: AM007 = %d | AM016 = %d | both, same direction = %d of 50\n",
            length(sig_007), length(sig_016), length(sig_both)))

## ----------------------------------------------------------------- 6. plot
if (length(sig_both) == 0) {
  cat("\nno Hallmark significant in both lines in the same direction\n")
} else {
  order_nes <- res_hallmark |>
    filter(pathway_label %in% sig_both) |>
    group_by(pathway_label) |>
    summarise(mean_nes = mean(NES), .groups = "drop") |>
    arrange(mean_nes) |> pull(pathway_label)

  plot_df <- res_hallmark |>
    filter(pathway_label %in% sig_both) |>
    mutate(pathway_label = factor(pathway_label, levels = order_nes),
           cell_line = factor(cell_line, levels = c("AM007", "AM016")),
           sig_label = case_when(padj < 0.001 ~ "***",
                                 padj < 0.01  ~ "**",
                                 padj < 0.05  ~ "*",
                                 TRUE         ~ ""))

  # Title Case, then fix the acronyms it mangles ("Tnfa" -> "TNFA"). Applied
  # to the factor LEVELS, so the NES ordering above is preserved. Entries not
  # present in the current set are no-ops.
  acronyms <- c("Tnfa" = "TNFA", "Nfkb" = "NF-kB", "Tgf" = "TGF",
                "E2f" = "E2F", "G2m" = "G2M", "P53" = "p53", "Dna" = "DNA",
                "Mtorc1" = "mTORC1", "Il6" = "IL6", "Il2" = "IL2",
                "Jak" = "JAK", "Stat3" = "STAT3", "Stat5" = "STAT5",
                "Kras" = "KRAS", "Myc" = "MYC", "Uv" = "UV", "Wnt" = "WNT",
                "Pi3k" = "PI3K", "Akt" = "AKT", "Mtor" = "mTOR",
                "Reactive Oxygen Species" = "Reactive oxygen species")

  pretty_labels <- gsub("\\b(\\w)", "\\U\\1", levels(plot_df$pathway_label), perl = TRUE)
  for (i in seq_along(acronyms))
    pretty_labels <- gsub(paste0("\\b", names(acronyms)[i], "\\b"),
                          acronyms[[i]], pretty_labels)
  levels(plot_df$pathway_label) <- pretty_labels

  height_hallmark <- max(4, 1.4 + length(sig_both) * 0.32)

  # Symmetric x limits put zero in the centre; the 12% headroom keeps the
  # significance asterisks from being clipped at the panel edge.
  nes_max <- max(abs(plot_df$NES)) * 1.12

  p <- ggplot(plot_df, aes(x = NES, y = pathway_label, fill = cell_line)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68,
             colour = "white", linewidth = 0.15) +
    geom_text(aes(label = sig_label, hjust = ifelse(NES >= 0, -0.15, 1.15)),
              position = position_dodge(width = 0.75), size = 2.8, colour = "grey25") +
    geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.4) +
    scale_fill_manual(values = col_main) +
    scale_x_continuous(limits = c(-nes_max, nes_max),
                       breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = 0)) +
    labs(x = "Normalised enrichment score (NES)", y = NULL,
         title = "MSigDB Hallmark pathways altered by SF3B4-KD",
         subtitle = sprintf("%d of 50 hallmarks significant (padj < 0.05) in both cell lines, in the same direction",
                            length(sig_both)),
         caption = "fgsea, ranked by sign(logFC) x -log10(P).   * padj<0.05   ** padj<0.01   *** padj<0.001") +
    theme_classic(base_size = 10) +
    theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.line.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y  = element_text(size = 9, colour = "grey15"),
          axis.line.x  = element_line(colour = "grey40", linewidth = 0.4),
          axis.title.x = element_text(size = 9, colour = "grey30",
                                      margin = margin(t = 6)),
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey35",
                                       margin = margin(b = 10)),
          plot.caption  = element_text(size = 7.5, colour = "grey45", hjust = 0),
          legend.title     = element_blank(),
          legend.position  = "top",
          legend.justification = "left",
          legend.key.size  = unit(0.45, "cm"),
          legend.text      = element_text(size = 9),
          plot.margin      = margin(10, 16, 10, 10))

  ggsave(file.path(OUT_DIR, "hallmark_fgsea_significant_both.svg"), p,
         width = 8.5, height = height_hallmark, dpi = 200)
  cat("\nsaved: hallmark_fgsea_significant_both.svg\n")
}
