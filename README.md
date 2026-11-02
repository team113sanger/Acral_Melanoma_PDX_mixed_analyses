# Acral Melanoma PDX — Genomic and Functional Analyses

Analysis scripts, metadata and results for the acral melanoma (AM) patient-derived
xenograft study, including the code that generates the figures of the manuscript
and the additional analyses prepared in response to peer review.

## Repository structure

```
INCA_res/
├── renv.lock          R package versions
├── data/              inputs (see Data sources)
├── results/           everything the scripts write
└── scripts/           one script per figure or analysis
```

Every script follows the same format: a header listing its INPUTS and OUTPUTS,
relative paths only, and assertions on the numbers quoted in the manuscript.

## Running

```bash
cd INCA_res
R -e 'renv::restore()'        # once
Rscript scripts/<script>.R
```

A script that stops with `'data/' not found` was run from the wrong directory.

### Paths

Scripts use relative paths (`IN <- "data"`, `OUT <- "results"`) and are intended
to be run from the project root. The `here` package is not used, and no absolute
paths to institutional file systems are recorded anywhere in the repository.

Do not set `.libPaths()` or `RENV_PROJECT` inside a script: either overrides the
environment recorded in `renv.lock`. `renv::snapshot()` and `renv::init()` do not
belong in a script either, since they rewrite the lockfile on every run.

## Figures

| script | figure |
|---|---|
| `Figure_1B_oncoprint.R` | 1B |
| `Figure_3B_CNApp_tumour_vs_pdx.R` | 3B |
| `Figure_4B_drugscreen_heatmap.R` | 4B |
| `Figure_S5C_signatures.R` | S5C, mutational signature panels |
| `Figure_S5C_annotation_tracks.R` | S5C, clinical annotation tracks |
| `Figure_S5E_G_H_splicing_panels.R` | S5E, S5G, S5H |
| `Figure_S5F_hallmark_fgsea.R` | S5F |
| `Figure_S21_essential_genes_split_bars.R` | S21 |

## Reviewer response

Analyses added during revision. These do not reproduce published panels; they
generate the figures and statistics included in the response to reviewers.

| script | analysis |
|---|---|
| `Reviewer3_Q1_expression_atm_chek1_bra_cohort.R` | ATM/CHEK1 copy-number loss versus expression, Brazilian cohort |
| `Reviewer3_Q1_ATM_CHEK1_hailstorm.R` | ATM/CHEK1 loss versus hailstorms
| `Reviewer_enrichment_two_tests.R` | essential genes among differentially spliced genes |

## Data sources

**Sequenza segmentation output** —
[Figshare](https://figshare.com/articles/dataset/Sequenza_segmentation_output_files/34745433)

**Copy-number scores (FCS/BCS/GCS)** — generated with CNApp, run from a local installation, not the web interface, using the SEQUENZA segmentation as input, default parameters, and the ploidy that SEQUENZA estimated for that same sample. The per-sample ploidy values supplied to CNApp are in data/copy_number_BRA_adjPloidy_06_10.tsv.

**Drug screen** — IC50 values and the concentration range assayed per compound
are in Extended Data Table S7.

**CRISPR dependency data** — DepMap Public 24Q4.

**Gene sets** — MSigDB Hallmark collection, via `msigdbr`.

## Environment

- R version: 4.4.1
- Package versions: `renv.lock`

```R
if (!require("renv")) install.packages("renv")
renv::restore()
```

## Authorship note

The drug-screen analysis (`Figure_4B_drugscreen_heatmap.R`) was developed by
Aguiar F.C., Carvalho D.G. and Sousa-Squiavinato A.C.M.

## Contact

**Annie Cristhine Moraes Sousa-Squiavinato**
annie.sousa@ensino.inca.gov.br · as81@sanger.ac.uk