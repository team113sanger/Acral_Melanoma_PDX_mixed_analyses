# Author: Annie Squiavinato
# 2026-10
#_______________________________________________________________________________
# Figure 5, panels e, g and h, all built from Extended Data Table 11 so the
# figure and the table cannot disagree.
#
#   e) differentially spliced events per cell line and shared, by event type,
#      plus the shared genes
#   g) overlap between CRISPR-essential genes and genes spliced in both lines
#   h) the 71 ES/IR events in 62 genes, by predicted protein impact
#
# Event type is taken from the EVENT ID prefix, not from the COMPLEX column.
# COMPLEX splits exon skipping into S, C1, C2, C3, ANN and MIC; filtering on
# COMPLEX == "S" alone drops 1,030 of the 2,678 ES events.
#
# INPUTS  (under IN_DIR)
#   Extended_Data _Table_11.xlsx                   3,795 events
#       must be the rebuilt copy. A version of this file was opened in Excel
#       under a pt-BR locale, which turned AGO1 and AGO4 into dates and 403
#       numeric cells into text; the checks below stop on that file.
#   Extended_Data_Table_9_shared_vulnerabilities.tsv   823 genes
#
# OUTPUTS (under OUT_DIR)
#   Figure_5e_events_and_genes.pdf
#   Figure_5g_venn.pdf
#   Figure_5h_predicted_impact.pdf
#   Figure_5_panel_counts.tsv     every number drawn, for the legend
#
# Package versions are managed with renv: run renv::restore() once.
# Run this script from the project directory that holds renv.lock.
#_______________________________________________________________________________

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(patchwork)

IN_DIR  <- "data"
OUT_DIR <- "results"

if (!dir.exists(IN_DIR))
  stop("'", IN_DIR, "' not found. Run this script from the project directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

S11_PATH <- file.path(IN_DIR, "SupplementaryTable_S11_events_REBUILT.xlsx")

s11 <- read_excel(S11_PATH)   
CRISPR_PATH <- file.path(IN_DIR, "Extended_Data_Table_9_shared_vulnerabilities.tsv")
for (f in c(S11_PATH, CRISPR_PATH)) if (!file.exists(f)) stop("missing input: ", f)

CLASSES <- c("Exon skipping (ES)", "Intron retention (IR)",
             "Alt 3' splice site", "Alt 5' splice site")
IMPACTS <- c("ORF disruption", "Alternative protein isoforms", "Noncoding")

COL_CLASS  <- c("Exon skipping (ES)"    = "#E8755C",
                "Intron retention (IR)" = "#4FB8A0",
                "Alt 3' splice site"    = "#7FC4DE",
                "Alt 5' splice site"    = "#2E4D7B")
COL_IMPACT <- c("ORF disruption"               = "#7B7FB5",
                "Alternative protein isoforms" = "#9FD9C4",
                "Noncoding"                    = "#E89279")

## ------------------------------------------------------------------ input
# skip = 1: the first row is a title, the real header is on the second.
s11 <- read_excel(S11_PATH, skip = 1)

# Guards against the Excel-corrupted copy. A gene symbol read as a date or a
# dPSI column read as text means the wrong file; the counts would still come
# out, just wrong.
dates <- names(s11)[vapply(s11, function(x)
  inherits(x, "Date") || inherits(x, "POSIXct") ||
    any(grepl("^20\\d\\d-\\d\\d-\\d\\d", as.character(x))), logical(1))]
if (length(dates))
  stop("date-like values in: ", paste(dates, collapse = ", "),
       "\nThis is the copy that was opened in Excel. Use the rebuilt table.")
if (!all(c("AGO1", "AGO4") %in% s11$GENE))
  stop("AGO1 and/or AGO4 missing: in a pt-BR locale Excel converts them to ",
       "dates. Use the rebuilt table.")

s11 <- s11 |>
  mutate(across(c(dPSI_AM007, dPSI_AM016), as.numeric),
         `EVENT CLASS` = factor(`EVENT CLASS`, levels = CLASSES))

stopifnot(nrow(s11) == 3795, !anyDuplicated(s11$EVENT), !anyNA(s11$`EVENT CLASS`))

# "Not annotated" is one label covering several loci and must never be counted
# as a gene.
genes_of <- function(df) setdiff(unique(df$GENE), "Not annotated")

essential <- read_tsv(CRISPR_PATH, show_col_types = FALSE)[[1]] |>
  as.character() |> trimws() |> unique()
essential <- essential[nzchar(essential) & !is.na(essential)]
if (length(essential) < 800) warning("expected 823 genes, got ", length(essential))

in_007 <- s11$`DETECTED IN` %in% c("AM007 only", "Both cell lines")
in_016 <- s11$`DETECTED IN` %in% c("AM016 only", "Both cell lines")
shared <- s11$`DETECTED IN` == "Both cell lines"

## ================================================================ panel e
events_e <- bind_rows(
  count(s11[in_007, ], `EVENT CLASS`, .drop = FALSE) |> mutate(set = "AM007"),
  count(s11[in_016, ], `EVENT CLASS`, .drop = FALSE) |> mutate(set = "AM016"),
  count(s11[shared, ], `EVENT CLASS`, .drop = FALSE) |> mutate(set = "shared")) |>
  mutate(set = factor(set, levels = c("AM007", "AM016", "shared")))

# Genes are counted per event class, so a gene carrying both an ES and an IR
# event appears in both segments. The distinct total is reported separately.
genes_e <- s11[shared, ] |>
  filter(GENE != "Not annotated") |>
  distinct(GENE, `EVENT CLASS`) |>
  count(`EVENT CLASS`, .drop = FALSE) |>
  mutate(set = "shared")

n_genes_distinct <- length(genes_of(s11[shared, ]))

p_ev <- ggplot(events_e, aes(set, n, fill = `EVENT CLASS`)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = COL_CLASS, name = "event type") +
  scale_y_continuous(limits = c(0, 3000), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Differentially spliced events") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        legend.position = "none")

p_gn <- ggplot(genes_e, aes(set, n, fill = `EVENT CLASS`)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = COL_CLASS, name = "event type") +
  scale_y_continuous(limits = c(0, 1500), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Differentially spliced genes") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"))

ggsave(file.path(OUT_DIR, "Figure_5e_events_and_genes.pdf"),
       p_ev + p_gn + plot_layout(widths = c(1.6, 1)),
       width = 6.2, height = 3.4)

## ================================================================ panel g
# Gene-level concordance: significant in AM007 AND in AM016, not necessarily
# through the same event.
spliced_both <- intersect(genes_of(s11[in_007, ]), genes_of(s11[in_016, ]))
overlap      <- intersect(essential, spliced_both)

venn <- tibble(
  x = c(-0.55, 0.55), y = 0, r = 1,
  lab = c("Essential\ngenes", "Differentially\nspliced genes"))
counts <- tibble(
  x = c(-1.00, 0, 1.00), y = 0,
  n = c(length(setdiff(essential, spliced_both)), length(overlap),
        length(setdiff(spliced_both, essential))))

circle <- function(x0, r, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = x0 + r * cos(t), y = r * sin(t))
}
circles <- bind_rows(
  mutate(circle(-0.55, 1), set = "Essential genes"),
  mutate(circle( 0.55, 1), set = "Differentially spliced genes"))

p_venn <- ggplot() +
  geom_polygon(data = circles, aes(x, y, group = set, fill = set),
               alpha = 0.55, colour = NA) +
  geom_text(data = counts, aes(x, y, label = n), size = 4.2) +
  geom_text(data = venn, aes(x, 1.25, label = lab), size = 3.6, lineheight = 0.95) +
  scale_fill_manual(values = c("Essential genes" = "#E8A79B",
                               "Differentially spliced genes" = "#A8D4E6")) +
  coord_equal(xlim = c(-2.1, 2.1), ylim = c(-1.3, 1.8)) +
  theme_void() + theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "Figure_5g_venn.pdf"), p_venn, width = 4, height = 2.6)

## ================================================================ panel h
# The cascade: essential and spliced -> same event in both lines -> ES or IR
# with an impact annotation available. "Impact available" is not the same as
# "functional consequence predicted": 5 of the 71 are non-coding.
both_ev <- filter(s11, CRISPR_essential == "Yes", `DETECTED IN` == "Both cell lines")
esir    <- filter(both_ev, `EVENT CLASS` %in% CLASSES[1:2], `IMPACT AVAILABLE` == "Yes")

impact_h <- esir |>
  mutate(`EVENT CLASS` = droplevels(`EVENT CLASS`),
         `PREDICTED IMPACT` = factor(`PREDICTED IMPACT`, levels = IMPACTS)) |>
  count(`EVENT CLASS`, `PREDICTED IMPACT`, .drop = FALSE)

p_imp <- ggplot(impact_h, aes(`EVENT CLASS`, n, fill = `PREDICTED IMPACT`)) +
  geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = COL_IMPACT, name = NULL) +
  scale_x_discrete(labels = c("Exon skipping (ES)" = "ES",
                              "Intron retention (IR)" = "IR")) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Differentially spliced events",
       subtitle = sprintf("(n = %d genes)", n_distinct(esir$GENE))) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "right", legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8.5))

ggsave(file.path(OUT_DIR, "Figure_5h_predicted_impact.pdf"), p_imp,
       width = 4.2, height = 3.2)

## ---------------------------------------------------- counts and assertions
out <- bind_rows(
  transmute(events_e, panel = "e", unit = "events", group = set,
            class = as.character(`EVENT CLASS`), n),
  transmute(genes_e,  panel = "e", unit = "genes",  group = set,
            class = as.character(`EVENT CLASS`), n),
  tibble(panel = "g", unit = "genes",
         group = c("essential only", "overlap", "spliced only"),
         class = NA_character_, n = counts$n),
  transmute(impact_h, panel = "h", unit = "events", group = as.character(`EVENT CLASS`),
            class = as.character(`PREDICTED IMPACT`), n))
write_tsv(out, file.path(OUT_DIR, "Figure_5_panel_counts.tsv"))

cat("\npanel e, events:\n"); print(pivot_wider(events_e, names_from = set, values_from = n))
cat("\npanel e, genes per class (distinct genes overall:", n_genes_distinct, ")\n")
print(genes_e)
cat(sprintf("\npanel g: essential %d | spliced in both %d | overlap %d\n",
            length(essential), length(spliced_both), length(overlap)))
cat(sprintf("panel h: %d genes, %d events (from %d genes / %d events with the same event)\n",
            n_distinct(esir$GENE), nrow(esir),
            n_distinct(both_ev$GENE), nrow(both_ev)))
cat("\n"); print(pivot_wider(impact_h, names_from = `PREDICTED IMPACT`, values_from = n))

# The numbers quoted in the figure and in the text.
stopifnot(
  sum(events_e$n[events_e$set == "AM007"])  == 2824,
  sum(events_e$n[events_e$set == "AM016"])  == 2397,
  sum(events_e$n[events_e$set == "shared"]) == 1426,
  length(essential) == 823, length(spliced_both) == 1158, length(overlap) == 81,
  n_distinct(both_ev$GENE) == 76, nrow(both_ev) == 98,
  n_distinct(esir$GENE)    == 62, nrow(esir)    == 71)

message("\nsaved three panels and Figure_5_panel_counts.tsv in ", OUT_DIR)
