# Author: Annie Squiavinato
# 2026-10-25
#_______________________________________________________________________________
# Overlap between CRISPR-essential genes and genes differentially spliced after
# SF3B4 knockout, by hypergeometric test and by matched permutation.
#
# The two tests use different gene counts because the permutation needs a
# measurable event plus expression and length for every gene:
#
#                                           essential   spliced   overlap
#   full lists                                   823     1,158        81
#   + event measurable in both cell lines        780     1,158        81
#   + expression and length also available       746     1,117        78
#
# INPUTS  (under IN_DIR)
#   Extended_Data_Table_9_shared_vulnerabilities.tsv   823 CRISPR-essential genes
#        (not all_common_essential_genes_ranking.tsv, which holds only 361)
#   dPSI_DiffEvents_AM_sgSF3B4_dPSI10_percondition.xlsx   sheets dPSI_AM007/AM016
#   IncTable_dPSI_values_AMcelllines_sgSF3B4_n356613.txt  all measurable events
#   AM202607_sgSF3B4_vs_sgCtrol_edgeR.xlsx                logCPM
#   ENSEMBL_v103_human_gene_biotypes_full_gene_information.txt   gene coordinates
#
# OUTPUT  enrichment_two_tests.tsv
#
# Package versions are managed with renv: run renv::restore() once.
# Paths are relative to the repository root; run the script from there.
#_______________________________________________________________________________

library(dplyr)
library(readr)
library(readxl)
library(tibble)

IN_DIR  <- "data"
OUT_DIR <- "results"
if (!dir.exists(IN_DIR))
  stop("'data/crispr-spicling' not found. Run this script from the INCA_res/ directory ",
       "(the one containing renv.lock). Working directory is: ", getwd())

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_EXPRESSED <- 16000
N_PERM      <- 10000
set.seed(1)

## ------------------------------------------------------------------ input
essential <- read_tsv(file.path(IN_DIR, "Extended_Data_Table_9_shared_vulnerabilities.tsv"),
                      show_col_types = FALSE)[[1]] |>
  as.character() |> trimws() |> unique()
essential <- essential[nzchar(essential) & !is.na(essential)]
if (length(essential) < 800)
  stop("expected ~823 genes, got ", length(essential),
       ": this looks like the 361-gene ranking file, not the shared-vulnerability list")

gene_of <- function(sheet) {
  read_excel(file.path(IN_DIR, "dPSI_DiffEvents_AM_sgSF3B4_dPSI10_percondition.xlsx"),
             sheet = sheet)$GENE |> as.character() |> trimws() |> unique()
}
spliced <- intersect(gene_of("dPSI_AM007"), gene_of("dPSI_AM016"))
spliced <- spliced[nzchar(spliced) & !is.na(spliced)]
overlap <- intersect(essential, spliced)

cat(sprintf("essential %d | spliced in both lines %d | overlap %d\n",
            length(essential), length(spliced), length(overlap)))

## ------------------------------------------- test 1: hypergeometric, N = 16,000
K1 <- length(essential); n1 <- length(spliced); k1 <- length(overlap)
e1 <- K1 * n1 / N_EXPRESSED

res1 <- tibble(test = "hypergeometric, all expressed genes",
               N = N_EXPRESSED, K = K1, n = n1, observed = k1,
               expected = round(e1, 1), fold = round(k1 / e1, 2),
               p_value = phyper(k1 - 1, K1, N_EXPRESSED - K1, n1, lower.tail = FALSE))

## ------------------------ test 2: permutation matched on expression and length
inc <- read_tsv(file.path(IN_DIR, "IncTable_dPSI_values_AMcelllines_sgSF3B4_n356613.txt"),
                show_col_types = FALSE,
                col_select = c(GENE, dPSI_AM007, dPSI_AM016))

# The IncTable is expected to hold only events with a dPSI in both lines;
# asserted rather than assumed, since the background depends on it.
if (anyNA(inc$dPSI_AM007) || anyNA(inc$dPSI_AM016)) {
  message("IncTable carries NA dPSI -- keeping only events measurable in both lines")
  inc <- filter(inc, !is.na(dPSI_AM007), !is.na(dPSI_AM016))
}
testable <- inc$GENE |> as.character() |> trimws() |> unique()
testable <- testable[nzchar(testable) & !is.na(testable)]

# A spliced gene missing from the IncTable means the gene symbols disagree
# between the two files, which would corrupt every count below.
fora <- setdiff(spliced, testable)
if (length(fora)) stop(length(fora), " spliced gene(s) absent from the IncTable: ",
                       paste(head(fora, 10), collapse = ", "))

edger <- file.path(IN_DIR, "AM202607_sgSF3B4_vs_sgCtrol_edgeR.xlsx")
cpm <- read_excel(edger, sheet = "AM007_sgSF3B4_vs_AM007_sgCtrol") |>
  select(GENE = NAME, a = logCPM) |> distinct(GENE, .keep_all = TRUE) |>
  inner_join(read_excel(edger, sheet = "AM016_sgSF3B4_vs_AM016_sgCtrol") |>
               select(GENE = NAME, b = logCPM) |> distinct(GENE, .keep_all = TRUE),
             by = "GENE") |>
  transmute(GENE, logCPM = (a + b) / 2)

# Genomic span. Where a symbol maps to several entries the longest is kept, so
# the covariate is defined the same way for every gene.
len <- read_tsv(file.path(IN_DIR,
                          "ENSEMBL_v103_human_gene_biotypes_full_gene_information.txt"),
                show_col_types = FALSE) |>
  filter(!is.na(hgnc_symbol)) |>
  transmute(GENE = trimws(hgnc_symbol), length = end_position - start_position) |>
  arrange(desc(length)) |> distinct(GENE, .keep_all = TRUE)

universe <- tibble(GENE = testable) |>
  inner_join(cpm, by = "GENE") |>
  inner_join(len, by = "GENE") |>
  mutate(is_ess = GENE %in% essential, is_ds = GENE %in% spliced)

N2 <- nrow(universe); K2 <- sum(universe$is_ess)
n2 <- sum(universe$is_ds); k2 <- sum(universe$is_ess & universe$is_ds)
cat(sprintf("\nmatched universe %d | essential %d | spliced %d | overlap %d\n",
            N2, K2, n2, k2))

cat("\ncovariate medians:\n")
print(universe |> group_by(is_ess) |>
        summarise(genes = n(), logCPM = round(median(logCPM), 2),
                  length_kb = round(median(length) / 1000, 1), .groups = "drop"))

# Deciles of each covariate give up to 100 strata. Each draw reproduces the
# observed set's composition across strata and samples without replacement
# within them, so expression and length cannot contribute to the null.
dec <- function(x) {
  r <- rank(x, ties.method = "first")
  cut(r, quantile(r, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
}
key  <- interaction(dec(universe$logCPM), dec(universe$length), drop = TRUE)
idx  <- split(seq_len(N2), key)
need <- table(key[universe$is_ds])
ess  <- universe$is_ess

short <- names(need)[need > lengths(idx)[names(need)]]
if (length(short)) stop("stratum too small to match: ", paste(short, collapse = ", "))

null <- replicate(N_PERM, sum(unlist(lapply(names(need), function(s)
  ess[sample(idx[[s]], need[[s]])]))))

res2 <- tibble(test = "permutation, matched on expression and length",
               N = N2, K = K2, n = n2, observed = k2,
               expected = round(mean(null), 1), fold = round(k2 / mean(null), 2),
               p_value = (sum(null >= k2) + 1) / (N_PERM + 1))

## ------------------------------------------------------------------ output
res <- bind_rows(res1, res2)
print(as.data.frame(res), row.names = FALSE, digits = 3)
write_tsv(res, file.path(OUT_DIR, "enrichment_two_tests.tsv"))

# Values quoted in the manuscript. Counts are exact; the permutation mean
# shifts by about 0.2 between seeds, so it is checked within tolerance.
stopifnot(res1$observed == 81, res1$expected == 59.6, res1$fold == 1.36,
          round(res1$p_value, 4) == 0.0027,
          res2$N == 10648, res2$K == 746, res2$n == 1117, res2$observed == 78,
          abs(res2$expected - 80.6) < 0.5, res2$p_value > 0.5)

cat("\nsaved:", file.path(OUT_DIR, "enrichment_two_tests.tsv\n"))
library(readr); library(dplyr)
inc <- read_tsv("data/IncTable_dPSI_values_AMcelllines_sgSF3B4_n356613.txt",
                show_col_types = FALSE)
names(inc)[1:6]                 # confira o nome da coluna de gene
universe <- n_distinct(inc$GENE)

spliced   <- 1158               # genes com evento significativo nas duas linhagens
essential <- 823
overlap   <- 81                 # essenciais entre os spliced

c(bg_rate        = spliced / universe,
  essential_rate = overlap / essential)

fisher.test(matrix(c(overlap,
                     essential - overlap,
                     spliced   - overlap,
                     universe - essential - spliced + overlap), 2))
library(readr); library(dplyr); library(readxl)

d <- read_excel("data/dPSI_DiffEvents_AM_sgSF3B4_dPSI10_percondition.xlsx")
names(d)                       # localize a coluna de dPSI e a de EVENT

# pegue um evento com |dPSI| grande
e <- d |> slice_max(abs(dPSI), n = 1)
e

inc <- read_tsv("data/IncTable_dPSI_values_AMcelllines_sgSF3B4_n356613.txt",
                show_col_types = FALSE)
inc |> filter(EVENT == e$EVENT) |> glimpse()   # veja as colunas de PSI por amostra

