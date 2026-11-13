# GLP2R expression in colorectal cancer

Code, processed results, tables and figures for:

> Lu X, Tao T, Li J, Cao H. *Reduced GLP2R expression in colorectal cancer: a confounder-adjusted multi-cohort analysis of its prognostic, immune and epigenetic correlates.* (manuscript submitted)

All analyses use public, de-identified data: TCGA-COAD/READ (GDC, UCSC Xena, cBioPortal, PanCanAtlas) and the GEO series GSE39582 and GSE17536.

## Repository layout

| Path | Contents |
|---|---|
| `05_Analysis/scripts/` | R pipeline, run in numerical order (see below) |
| `05_Analysis/results/` | Processed result tables (`.csv`) and analysis objects (`.rds`) written by the pipeline; `rm_*.csv` are per-module result ledgers |
| `05_Analysis/results_master.csv` | Every number quoted in the manuscript, with its key, value and description |
| `05_Analysis/data_manifest.csv` | Input files: source, URL, size and MD5; column `in_repository` marks what is included here |
| `05_Analysis/data_cache/` | Small inputs that are included (gene sets, marker lists, GDC file listings); see `data_cache/README.md` for the rest |
| `05_Analysis/sessionInfo.txt` | R session of the main analysis (R 4.4.1) |
| `input/Table_S1_GLP_genes.csv` | The 50-gene GLP-signalling candidate set |
| `results_presubmission/` | Output of `S01_presubmission_sensitivity.R` (tissue-source-site analysis, MSI-stratum heterogeneity, top-100 turquoise hub genes) |
| `02_Figures/` | Main Figures 1–7 and Supplementary Figures S1–S6 (PNG, 300 dpi, and PDF) |
| `03_Tables/` | Tables 1–2 and Supplementary Tables S1–S8 (CSV and one Excel workbook) |

## Requirements

- R 4.4.1 for the full pipeline (package versions in `05_Analysis/sessionInfo.txt`). `00_setup.R` to `00e_check_pkgs.R` install the CRAN and Bioconductor dependencies into `_work/Rlib`.
- Figures 5, 7 and S1 and `S01_presubmission_sensitivity.R` were run with R 4.6.1 and ggplot2 4.0.3. Under R ≥ 4.5, `S00_install_r46.R` installs the plotting packages into the user library.
- Two download steps call `curl` (`22a_gdc_paired_download.R`; `62_epi_extra.R` uses `curl.exe` on Windows).

## Setting up the paths

The scripts read their locations from environment variables:

| Variable | Meaning | Value for this repository |
|---|---|---|
| `GLP2R_RB` | root containing `05_Analysis/` | the repository folder |
| `GLP2R_SCRIPTS` | script folder | `<repo>/05_Analysis/scripts` |
| `GLP2R_OUT` | where figures, tables and logs are written | the repository folder (or any other folder) |
| `GLP2R_INPUT` | folder with `Table_S1_GLP_genes.csv` | `<repo>/input` |
| `GLP2R_ORIG` | folder with `TCGA_CRC_combined.rds` | wherever you place it |

Example (R):

```r
repo <- "C:/path/to/GLP2R"
Sys.setenv(GLP2R_RB = repo, GLP2R_OUT = repo,
           GLP2R_SCRIPTS = file.path(repo, "05_Analysis/scripts"),
           GLP2R_INPUT   = file.path(repo, "input"),
           GLP2R_ORIG    = "C:/path/to/tcga_input")
source(file.path(Sys.getenv("GLP2R_SCRIPTS"), "10_cohort.R"))
```

or from a shell, `Rscript --vanilla 05_Analysis/scripts/10_cohort.R` after setting the same variables.

## Pipeline

| Step | Scripts | Main outputs |
|---|---|---|
| Setup | `00_setup.R`, `00b`–`00e`, `S00_install_r46.R` | R package library |
| Public inputs | `02_extract_external.R`, `02b_extract_cbio_thorsson.R`, `03_geo_prepare.R` | files in `data_cache/` |
| Cohort | `10_cohort.R` | `cohort.rds`, `cohort_clinical_harmonised.csv` (591 primary tumours, 33 normal samples) |
| Expression | `20_deseq.R`, `21_expr_tests.R`, `22a_gdc_paired_download.R`, `22b_paired_deseq.R` | DESeq2 (unpaired and 32 matched pairs), VST matrices, expression tests |
| Survival | `30_survival.R`, `32_geo_glp2r.R` | nested Cox models M0–M2, sensitivity analyses, GEO cohorts, cross-cohort synthesis |
| Immune | `40_immune.R` | ssGSEA (28 cell types), MCP-counter, ESTIMATE, checkpoints, immune subtypes |
| Pathways | `50_pathway.R` | partial-correlation ranking, fgsea, compartment markers |
| Methylation, CNV, mutation | `60_epi_cnv.R`, `62_epi_extra.R`, `63_cimp.R` | promoter methylation and decision rule, GISTIC2, MC3, TMB, CIMP/BRAF/MMR |
| Signature | `70_lasso.R` | LASSO-Cox signature, bootstrap stability, permutation reference, external evaluation |
| Network | `80_wgcna.R` | signed WGCNA and module annotation |
| Figures | `90_fig_helpers.R` (sourced), `91_fig1.R`–`98_figS.R` | `02_Figures/` |
| Ledgers and tables | `99a_results_master.R`, `99b_tables.R` | `results_master.csv`, `03_Tables/` |
| Sensitivity analyses | `S01_presubmission_sensitivity.R` | `results_presubmission/` |

The random seed is 20240617 throughout.

Because `05_Analysis/results/` is included, the figure scripts can be re-run without repeating the analyses, provided the derived matrices they read from `data_cache/` (`cohort.rds`, `vst_tumor_591.rds`, `vst_all_624.rds`, `geo_cohorts.rds`) have been rebuilt or downloaded.

## Notes

- `TCGA_CRC_combined.rds` (gene-level STAR counts and clinical data for the 625 TCGA-COAD/READ samples from the GDC; 210 MB) is not included because it exceeds GitHub's file-size limit. It can be rebuilt from the GDC or obtained from the corresponding author.
- `30_survival.R` and `70_lasso.R` contain two reconciliation steps that read objects from an earlier internal analysis (`cox_full.rds`, `lasso_combined.rds`). These objects are not distributed; the steps are skipped when the files are absent and do not affect any reported result.
- Supplementary Tables S2b, S3b–S3d, S4c, S6a–S6b, S7g and S8h in `03_Tables/` contain columns or rows added after `99b_tables.R` was run (Fisher z comparison of MSI strata, Hartung–Knapp interval, covariate descriptions and the top-100 hub table from `S01_presubmission_sensitivity.R`). The versions in `03_Tables/` are the ones reported in the manuscript.

## Licence

Code: MIT License (`LICENSE`). Results, tables and figures: CC BY 4.0 (`LICENSE-data.md`). Public input data remain under the terms of their original providers.

## Contact

Xiangchen Lu, Department of Critical Care Medicine, Chengdu Xinhua Hospital Affiliated to North Sichuan Medical College, Chengdu 610066, China. E-mail: luxc27@live.hk
