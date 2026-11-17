# Reduced GLP2R expression in colorectal cancer: a confounder-adjusted multi-cohort analysis — reproducibility package

[![DOI](https://zenodo.org/badge/DOI/ZENODO_DOI_PLACEHOLDER.svg)](https://doi.org/ZENODO_DOI_PLACEHOLDER)

- **Repository:** https://github.com/LXC27/GLP2R
- **Version:** 2.0.0 (2026-11-17)
- **Companion manuscript:** *Reduced GLP2R expression in colorectal cancer: a confounder-adjusted multi-cohort analysis of its prognostic, immune and epigenetic correlates*
- **Code licence:** MIT — see `LICENSE`
- **Licence for derived results and tables:** CC BY 4.0 — see `LICENSE-DATA.md`

This repository contains the analysis code, the numeric result ledger, the manuscript table files and the pre-submission sensitivity outputs for the analysis reported in the manuscript. It is intended to make the reported numbers traceable to code and to allow the analyses to be re-run where the required inputs can be supplied.

**The package is about 6 MB (133 files) and can be cloned directly. It deliberately excludes the analysis data cache (≈1.1 GB), the R library (≈1.1 GB), intermediate `.rds` objects and all third-party raw data.** Read [What can and cannot be reproduced from this archive](#what-can-and-cannot-be-reproduced-from-this-archive) before attempting a re-run.

---

## Authors

Xiangchen Lu¹, Xiaoru Sun², Jinming Li³, Xin Yang², Tianyi Tao⁴,*, Haiquan Cao⁵,*

1. Chengdu Xinhua Hospital Affiliated to North Sichuan Medical College, Chengdu 610066, China
2. Department of Otolaryngology Head & Neck Surgery, West China Hospital, Sichuan University, Chengdu 610041, China
3. Department of Colorectal Surgery, Fudan University Shanghai Cancer Center, Shanghai 200032, China
4. West China Institute for Smart Healthcare, West China Second University Hospital, Sichuan University, Chengdu 610041, China
5. Department of Critical Care Medicine, Beijing Anzhen Nanchong Hospital of Capital Medical University & Nanchong Central Hospital, Nanchong 637000, China

\* Correspondence: Tianyi Tao (tytao@scu.edu.cn) and Haiquan Cao (caohaiquan1978@163.com).

ORCID iDs: Xiangchen Lu 0009-0001-8824-8961 · Xiaoru Sun 0000-0001-8623-6935 · Jinming Li 0000-0002-4870-9815 · Xin Yang 0009-0003-3319-5094 · Tianyi Tao 0009-0007-9481-7909 · Haiquan Cao 0009-0000-3434-7932.

**Author contributions (CRediT).** Xiangchen Lu — conceptualization, methodology, formal analysis, data curation, software, visualization, writing–original draft; Xiaoru Sun — methodology, formal analysis, validation; Jinming Li — data curation, validation, resources, software, visualization; Xin Yang — data curation, investigation, resources; Tianyi Tao — conceptualization, project administration, writing–review and editing; Haiquan Cao — supervision, writing–review and editing.

**Declarations.** Funding — this study received no specific grant from any funding agency. Competing interests — the authors declare none. Ethics — the analyses used only publicly available, de-identified data from TCGA and GEO, collected under the approvals and consent procedures of the original studies, so no additional ethical approval was required. Data and code availability — the analysis scripts, software versions (sessionInfo) and data manifest are provided in this repository (https://github.com/LXC27/GLP2R).

---

## Study summary

The study analyses GLP2R (glucagon-like peptide-2 receptor) expression in colorectal cancer (CRC) across one TCGA cohort and two microarray cohorts, and asks which of the previously reported prognostic and immune associations **persist after adjustment** for microsatellite instability (MSI), tumour purity, primary site and clinical covariates, and whether the GLP2R promoter is hypermethylated.

**Cohorts.** TCGA-COADREAD: 591 primary tumours (431 colon, 160 rectum), 33 solid-tissue normals and 32 matched normal/tumour pairs (retrieved from the GDC). Overall survival comes from the TCGA-CDR (566 patients with OS, 114 deaths). External validation: **GSE39582** (556 with OS) and **GSE17536** (177 with OS). Only these two GEO series are analysed; no other GEO series is used by any script.

**Main findings.** GLP2R mRNA is strongly lower in tumours than normal mucosa (DESeq2 log2FC −4.07, FDR 2.2×10⁻⁶⁶; lower in 32/32 matched pairs) and lower in MSI-high (P = 3.5×10⁻¹⁷) and CIMP-high tumours. Higher expression is associated with longer overall survival **before** adjustment (HR per SD 0.82, 95% CI 0.68–0.98); adjusted for age, sex, stage, site and MSI the estimate is 0.86 (0.69–1.07), with little added discrimination (C-index 0.762 → 0.767). Age/sex/stage-adjusted estimates are 0.83 (TCGA), 1.00 (GSE39582) and 0.71 (GSE17536). Adjusted ssGSEA immune correlations are small (|ρ| ≤ 0.24); the strongest correlate is the MCP-counter endothelial score (ρ = 0.42), and GLP2R belongs to a stromal co-expression module. Eight promoter CpGs are not hypermethylated (Δβ = −0.03) and copy number is unrelated to expression after adjustment. An **eight-gene** LASSO-Cox signature does not generalise (held-out C-index 0.518).

**Pre-submission sensitivity analyses** (`results/S01_presubmission_sensitivity.txt` and the `results/A_*`, `D_*`, `E_*` files):
- **A — tissue-source-site.** Restricting to the 312 tumours / 33 normals from tissue-source sites that contributed normals, the tumour–normal difference is unchanged (β = −4.04, AUC 0.999).
- **B — C-index comparators.** With age/sex/stage/site (+ MSI), adding GLP2R raises C-index 0.7625 → 0.7670 (LRT P = 0.18).
- **C — Hartung–Knapp.** The random-effects estimate is HR 0.856 (0.705–1.039); the Hartung–Knapp CI widens to 0.559–1.310 (P = 0.256).
- **D — MSI-stratified heterogeneity.** Immune correlations are examined separately within MSI strata.
- **E — leave-GLP2R-out eigengene.** Removing GLP2R leaves the turquoise-module eigengene essentially unchanged (partial ρ 0.252 → 0.251); GLP2R's kME rank is 1516/1747, i.e. it is not a hub.
- **F — MLH1 promoter-restricted control.** MLH1 methylation–expression coupling is preserved in the promoter window (ρ = −0.43), validating the negative GLP2R result.

---

## What is in this repository

```
GLP2R/
├── README.md                  this file
├── RUNNING.md                 step-by-step instructions to re-run the scripts
├── UPLOAD_STEPS.md            upload instructions (Chinese) for GitHub + Zenodo
├── CITATION.cff               citation metadata (CFF 1.2.0)
├── LICENSE                    MIT licence (code)
├── LICENSE-DATA.md            CC BY 4.0 and third-party data terms (results/tables)
├── .gitignore                 excludes data cache, R library and intermediate objects
├── scripts/                   43 files: portable R analysis scripts (no PowerShell helper)
├── results/                   47 files: result ledger, data manifest, session info, 43 result CSVs + 4 pre-submission outputs
└── tables/                    35 files: 34 table CSVs + the combined Excel workbook
```

### `scripts/` (43 files)

Portable copies of the analysis scripts. Two kinds of change were made relative to the original analysis workspace:

1. **Configurable paths.** Hard-coded `W:/Re_build`, `W:/Supplementary` and `V:/` drive paths are replaced by environment variables with the original values as defaults: `GLP2R_RB` (default `W:/Re_build`), `GLP2R_SCRIPTS` (default `<GLP2R_RB>/05_Analysis/scripts`), `GLP2R_OUT` (default `<GLP2R_RB>`), `GLP2R_INPUT` (default `W:/Supplementary`), `GLP2R_ORIG` (default `V:/`).
2. **R ≥ 4.5 compatibility.** Every script that prepends the project-local library now branches on `getRversion()`: under R < 4.5 it uses `_work/Rlib` (compiled for R 4.4.1); under R ≥ 4.5 it uses the user library `R_LIBS_USER` instead. This avoids the `rlang.dll`/`SETLENGTH` load failure that occurs when an R 4.4.1-compiled library is loaded by a newer R.

| Stage | Files | Purpose |
|---|---|---|
| Environment | `S00_install_r46.R`, `00_setup.R`, `00b_install_extra.R`, `00c_install_annotation.R`, `00d_install_dose.R`, `00e_check_pkgs.R` | Install CRAN/Bioconductor packages (`S00` is the minimal set; `00_setup`+`00b`/`00c`/`00d` is the full set) and check them. |
| Shared code | `00_common.R` | Paths, seed (20240617), plotting theme and helper functions; sourced by every analysis script. |
| Inventory / exploration | `01a_rdata_inventory.R`, `01b_explore_combined.R`, `01b2_explore_fields.R` | Inventory and exploration of the author-supplied `V:/` R objects (not needed to reproduce the reported models). |
| External data extraction | `02_extract_external.R`, `02b_extract_cbio_thorsson.R`, `03_geo_prepare.R` | Build the Xena/cBioPortal/Thorsson and GEO caches. |
| Cohort build | `10_cohort.R` | Harmonise clinical variables and build `cohort.rds`. |
| Differential expression | `20_deseq.R`, `22a_gdc_paired_download.R`, `22b_paired_deseq.R`, `21_expr_tests.R` | Tumour-vs-normal and paired DESeq2, VST matrices, GLP2R expression statistics. |
| Survival | `30_survival.R`, `32_geo_glp2r.R` | TCGA Cox models, time-dependent AUC, subgroups, GEO Cox models and meta-analysis. |
| Immune | `40_immune.R` | ssGSEA/GSVA (28 Charoentong signatures), MCP-counter, ESTIMATE, checkpoints, Thorsson subtypes. |
| Pathways / modules | `50_pathway.R`, `80_wgcna.R` | fgsea against MSigDB collections, co-expression correlations, WGCNA modules. |
| Epigenetics / genomics | `60_epi_cnv.R`, `62_epi_extra.R`, `63_cimp.R` | HM450 methylation, GISTIC copy number, DNMT/TET, TMB, mutations and CIMP. |
| Signature model | `70_lasso.R` | LASSO-Cox signature with bootstrap stability, permutation null and external validation. |
| Figures | `90_fig_helpers.R`, `91_fig1.R` … `97_fig7.R`, `98_figS.R` | Main Figures 1–7 and Supplementary Figures S1–S6 (PNG + PDF, 300 dpi). |
| Result ledger and tables | `99a_results_master.R`, `99b_tables.R`, `99d_crossref_refs.R` | Aggregate the result ledger, build the table CSVs/Excel workbook and check reference DOIs. |
| Manuscript / revision build | `100a_refdb.R`, `100_build_manuscript.R`, `110_build_response_docs.R` | Build the manuscript docx/PDF and response documents from text sources in the full analysis workspace. |
| Pre-submission sensitivity | `S01_presubmission_sensitivity.R` | The six sensitivity blocks A–F described above (outputs `results/S01_presubmission_sensitivity.txt` and the `A_*`/`D_*`/`E_*` CSVs). |

The PowerShell helper `100b_word_finish.ps1` from the original workspace is **not** included in this release (the portable package ships R scripts only; see `RUNNING.md`).

### `results/` (47 files)

- `results_master.csv` — **the single source of every number reported in the manuscript's main analysis** (360 key–value rows in 13 modules; columns `module,key,value,description`). Note: this ledger reflects the 2026-10-11 production run; the pre-submission sensitivity analyses are reported in the separate files below and are not folded into this ledger.
- `data_manifest.csv` — 65 entries describing every input in the data cache: file name, location, source, public URL, size, modification date and MD5 checksum.
- `sessionInfo.txt` — the R session and package versions used for the reported main results.
- 40 result CSVs produced by the analysis scripts (expression statistics, Cox models, immune/pathway/methylation/mutation tables, LASSO and WGCNA outputs, and the `rm_*.csv` fragments from which `results_master.csv` is assembled).
- 4 pre-submission outputs: `A_tumour_normal_by_TSS.csv`, `D_MSI_strata_difference.csv`, `E_TableS8h_turquoise_top100_hubs.csv`, `S01_presubmission_sensitivity.txt` (written by `S01_presubmission_sensitivity.R`).

### `tables/` (35 files)

The manuscript's tables: `Table1_cohort_characteristics.csv`, `Table2_cox_models.csv` (+ long format), the supplementary tables `TableS1…TableS8` as individual CSVs (34 CSVs in total), and `Tables_Table1-2_and_S1-S8_presubmission.xlsx`, the combined workbook. The LASSO tables (`TableS3b`, `TableS3c`, `TableS3d`) encode the eight-gene signature and `TableS8h_WGCNA_turquoise_top100_hubs.csv` the top-100 hub definition used in the manuscript.

Figures are **not** included in this repository (they are vector/raster artefacts of the analysis); the code that produces them is in `scripts/91_fig1.R`–`98_figS.R`.

---

## Data used in the study

All data are public except one author-assembled object (the TCGA count/clinical matrix). The URL recorded for each source and the expected cache file name are listed in `results/data_manifest.csv`. The scripts expect these files to be placed in `05_Analysis/data_cache/` under the names used below; **the scripts do not download most of them** (see the reproduction caveat).

| Data | Accession / dataset | Source (as recorded in `data_manifest.csv`) | Cache file(s) expected by the scripts |
|---|---|---|---|
| TCGA-COADREAD RNA-seq counts + clinical (primary input) | TCGA-COADREAD via GDC (STAR-Counts, GENCODE v36) | Author-assembled `TCGA_CRC_combined.rds` (GDC-derived; not a public download; size and MD5 in `data_manifest.csv`) | `V:/TCGA_CRC_combined.rds` |
| TCGA paired normal/tumour counts | GDC API (`api.gdc.cancer.gov/files`, `STAR - Counts`, open) | https://api.gdc.cancer.gov | `gdc_paired/*.tsv`, `gdc_paired_manifest.csv` (downloaded by `22a_gdc_paired_download.R`) |
| TCGA-COADREAD clinical, GISTIC2-thresholded copy number, HM450 methylation, 450k probe map | UCSC Xena legacy TCGA hub (`TCGA.COAD.sampleMap`, `TCGA.READ.sampleMap`) | https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/ | `COAD_clinicalMatrix.tsv`, `READ_clinicalMatrix.tsv`, `COAD_gistic_thresholded.tsv.gz`, `READ_gistic_thresholded.tsv.gz`, `COAD_HM450.tsv.gz`, `READ_HM450.tsv.gz`, `illuminaMethyl450_probeMap.tsv` |
| PanCancer Atlas clinical and MC3 mutation calls | cBioPortal study `coadread_tcga_pan_can_atlas_2018` | https://www.cbioportal.org/api | `cbio_pancan_sample.json`, `cbio_pancan_patient.json`, `cbio_mutations_selected.json` |
| TCGA-CDR survival endpoints and immune subtypes | PanCanAtlas immune landscape (Thorsson et al. 2018; TCGA-CDR, Liu et al. 2018) | https://gdc.cancer.gov/about-data/publications/panimmune | `Thorsson_2018_immune_landscape.xlsx`, `Thorsson_Subtype_Immune_Model_Based.txt.gz` |
| ABSOLUTE purity/ploidy master calls | GDC PanCanAtlas supplement (Taylor et al. 2018) | https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5 | `ABSOLUTE_mastercalls.txt` |
| GSE39582 expression + phenotype | GEO GSE39582 (Affymetrix HG-U133 Plus 2.0) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz | `GSE39582_series_matrix.txt.gz` |
| GSE17536 expression + phenotype | GEO GSE17536 (Affymetrix HG-U133 Plus 2.0) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17536/matrix/GSE17536_series_matrix.txt.gz | `GSE17536_series_matrix.txt.gz` |
| Affymetrix GPL570 platform annotation | GEO platform GPL570 | NCBI GEO (https://ftp.ncbi.nlm.nih.gov/geo) | `GPL570_table.txt` (needs `ID` and `Gene Symbol` columns) |
| 28 immune cell-type signatures | Charoentong et al. 2017, Table S6 gene sets | https://github.com/msfuji/immunoduct | `Charoentong2017.gmt` |
| MCP-counter marker signatures | Becht et al. 2016 | https://github.com/ebecht/MCPcounter | `MCPcounter_genes.txt`, `MCPcounter_probesets.txt` |
| MSigDB collections (Hallmark, KEGG, GO-BP, C8) | MSigDB via R package `msigdbr` | R package data (no cache file) | — |

Notes:

- The cBioPortal clinical JSON files, Xena downloads, GEO series matrices, the GPL570 table, the Thorsson workbook, the ABSOLUTE calls, the Charoentong gene sets, the MCP-counter files and the 32 GDC paired count files were placed in the data cache **before** the analysis scripts were run; the archive does not contain the download commands for the non-GDC files. `22a_gdc_paired_download.R` (GDC paired counts), `62_epi_extra.R` (cBioPortal mutation calls) and `99d_crossref_refs.R` (Crossref) are the only scripts in the archive that fetch data from the internet.
- `TCGA_CRC_combined.rds` is the primary input to the cohort build (`10_cohort.R`) and is listed in the data manifest with its MD5 (`431ae1e443d3438c3f6e398a70e0fb58`). The code that assembled this object from GDC files is **not** part of the archive.
- The TCGA/GEO/cBioPortal/Xena data remain subject to their original terms of use (see `LICENSE-DATA.md`).

**Re-creating the cache.** Only three archived scripts fetch data from the internet (`22a_gdc_paired_download.R` for the GDC paired counts, `62_epi_extra.R` for cBioPortal mutation calls, `99d_crossref_refs.R` for reference metadata); everything else must be assembled with the providers' tools. For TCGA-COADREAD, export the `STAR - Counts` gene-expression quantification files (GENCODE v36, open access) for project TCGA-COADREAD through the GDC Data Portal/API (`https://api.gdc.cancer.gov/`), then build the gene × sample count matrix (the assembly step itself is not archived); the exact GDC API query used for the 32 matched normal/tumour pairs is shown verbatim in `scripts/22a_gdc_paired_download.R`. The Xena files are the `TCGA.COAD.sampleMap` and `TCGA.READ.sampleMap` datasets (clinical matrix, GISTIC2 thresholded copy number, HM450 methylation) available from https://xenabrowser.net; the cBioPortal files are the clinical attributes and MC3 mutation calls of study `coadread_tcga_pan_can_atlas_2018`; the GEO files are the two series matrices linked above plus the GPL570 platform table; the Thorsson and ABSOLUTE files come from the PanCanAtlas pages/UUID listed above. `data_manifest.csv` records the expected size and MD5 of every cache file, which is the practical checklist for this step.

---

## Software environment

- **Main analysis (results, tables, figures 1–4, 6):** R 4.4.1 (2024-06-14 ucrt) on Windows 11 x64 (build 26200), locale `Chinese (Simplified)_China.utf8`, time zone `Asia/Shanghai`, random seed **20240617**. Full session information is in `results/sessionInfo.txt`.
- **Pre-submission sensitivity + regenerated figures 5, 7, S1–S6:** R 4.6.1.

| Package | Version | Package | Version |
|---|---|---|---|
| DESeq2 | 1.46.0 | survival | 3.6.4 |
| glmnet | 4.1-10 | GSVA | 2.0.7 |
| fgsea | 1.32.4 | WGCNA | 1.74 |
| timeROC | 0.4.1 | survminer | 0.5.2 |
| mice | 3.19.0 | maxstat | 0.7-26 |
| msigdbr | 26.1.0 | tidyestimate | 1.1.1 |
| ggplot2 | 4.0.3 | patchwork | 1.3.2 |
| officer | 0.7.3 | flextable | 0.9.11 |
| GEOquery | 2.74.0 | limma | 3.62.2 |
| clusterProfiler | 4.14.6 | | |

Also used (per `sessionInfo.txt`): `Matrix` 1.7-0, `fastcluster` 1.3.0, `dynamicTreeCut` 1.63-1, `ggpubr` 0.6.3, `readxl`, `DOSE` 4.0.1 and `org.Hs.eg.db`. The scripts install these from CRAN/Bioconductor (`00_setup.R`, `00b_install_extra.R`, `00c_install_annotation.R`, `00d_install_dose.R`, or minimally `S00_install_r46.R`).

External executables required by individual scripts:

- **`gzip`** on `PATH` — `02_extract_external.R` streams the GISTIC files through `gzip -dc` (a pure-R `gzfile` fallback exists).
- **`curl` / `curl.exe`** — GDC queries and downloads (`22a_gdc_paired_download.R`), cBioPortal mutation calls (`62_epi_extra.R`), Crossref lookups (`99d_crossref_refs.R`).

**Verification performed for this release:** all 43 R scripts were parsed with `Rscript` (R 4.6.1) and all 43 parse without error. This is a syntax check only — the analyses themselves were not re-executed for this archive. The production run that generated the included main results used R 4.4.1 as recorded in `sessionInfo.txt`, so exact numeric output on a different R/package version is not guaranteed.

---

## How to re-run the analysis

Full commands, directory layout, required inputs and troubleshooting are in **`RUNNING.md`**. In short:

1. Install R 4.4.1 (for the main analysis) or R ≥ 4.5 (which the portable scripts support via the user library). Install the packages: `S00_install_r46.R` for the minimal set, or `00_setup.R` + `00b`/`00c`/`00d` for the full set.
2. Set the environment variables (`GLP2R_RB`, `GLP2R_SCRIPTS`, `GLP2R_OUT`, `GLP2R_INPUT`, `GLP2R_ORIG`) or map the `W:`/`V:` drives to match the defaults (see `RUNNING.md`).
3. Place the public downloads under `05_Analysis/data_cache/` using the names in the table above.
4. Run the scripts in the order given in `RUNNING.md` (numeric order, with three documented dependency deviations).
5. Compare the outputs with the included `results/` CSVs and `tables/` files; `results_master.csv` is the numeric reference.

**Expected runtime.** The archived production logs show the whole main analysis on 2026-10-11 running from 02:46 to 04:55 (about 2 h 10 min). A fresh run additionally needs package installation and the manual download of the multi-hundred-MB cache, so allow several hours on a clean machine. (These times are indicative, not a benchmark.)

---

## Script execution order

The analysis workspace README says the scripts are "全部分析与构建脚本（按编号顺序运行）" — all analysis and build scripts, **to be run in numeric order**. Three deviations are required by the code itself:

1. **Install before analysis.** `00_setup.R` (and `00b`/`00c`/`00d`) install the package library. On a fresh machine, install first (or run `S00_install_r46.R` for the minimal set).
2. **`21_expr_tests.R` must run after `22a`/`22b`.** `21_expr_tests.R` reads `results/paired_glp2r_vst.csv`, which is written by `22b_paired_deseq.R`; `22b` in turn needs the GDC files downloaded by `22a`. The correct sub-order for steps 20–21 is: **20 → 22a → 22b → 21**.
3. **`99a`/`99b` need one extra pass.** `99b_tables.R` writes `rm_table1.csv` and `rm_table_pfi.csv` and reads `results_master.csv` at start; `99a_results_master.R` aggregates every `rm_*.csv` into `results_master.csv`. Run **99a → 99b → 99a again** so the final ledger also contains the `table1`/`table_pfi` modules. `91`–`98` must run before the first `99a`, because `98_figS.R` writes the `rm_figS.csv` fragment.

The pre-submission sensitivity script `S01_presubmission_sensitivity.R` runs after the main pipeline (it reads the main result objects) and writes `results/S01_presubmission_sensitivity.txt` plus the `A_*`/`D_*`/`E_*` CSVs.

All other scripts are independent given their inputs; `RUNNING.md` gives the recommended full sequence.

---

## What can and cannot be reproduced from this archive

This section is deliberately explicit. The package was assembled from the author's analysis workspace; the archive contains code and reported results, not the full runnable data environment.

### What can be reproduced from this archive

1. **Every reported number can be checked against the included tables.** `results/results_master.csv` is the ledger from which the manuscript's main-analysis numbers were taken, and `tables/` contains the final table files (including the combined Excel workbook). Each key–value row is labelled (`module`, `key`, `description`).
2. **Every number can be traced to the code.** The `rm_*.csv` fragments and the scripts that write them are included, and the script names, function calls and seed (20240617) are unchanged (modulo the portability edits described above, which touch only path resolution and `.libPaths()`).
3. **The analysis code can be audited and re-executed step by step** if the inputs are supplied (see "Cannot" below). The 43 R scripts are plain text.
4. **Individual downstream steps can be re-run** once the missing inputs are supplied. `RUNNING.md` lists the prerequisite of each stage.

### What cannot be reproduced from this archive

1. **A fresh end-to-end run from this repository alone.** The archive excludes:
   - the analysis data cache (`05_Analysis/data_cache/`, **156 files, ≈1.1 GB**), which holds the Xena/GEO/cBioPortal/Thorsson/GDC downloads *and* the derived intermediates (`cohort.rds`, `vst_all_624.rds`, `vst_tumor_591.rds`, `de_tumor_normal.rds`, `geo_cohorts.rds`, `hm450_subset.rds`, …) that later scripts read;
   - the R package library (`_work/Rlib`, ≈1.16 GB), which can be rebuilt from `sessionInfo.txt` but not without time and network access;
   - the intermediate result `.rds` objects and the figure files.
2. **The primary TCGA matrix cannot be regenerated from this archive.** `TCGA_CRC_combined.rds` (MD5 `431ae1e443d3438c3f6e398a70e0fb58`) is an author-assembled object read at `10_cohort.R`. It is not distributed here, and the code that assembled it from GDC STAR-Counts is not among the scripts (only the 32-pair downloader `22a_gdc_paired_download.R` is).
3. **One author-supplied supplementary CSV is missing**, so some steps cannot run as archived: `W:/Supplementary/Table_S1_GLP_genes.csv` (read by `21_expr_tests.R` and `70_lasso.R`). Only its derivatives are included.
4. **Most external inputs are not downloaded by the scripts.** Only `22a_gdc_paired_download.R`, `62_epi_extra.R` and `99d_crossref_refs.R` fetch data. The Xena, GEO, GPL570, Thorsson, ABSOLUTE, cBioPortal clinical, Charoentong and MCP-counter files must be obtained manually.
5. **The manuscript and response documents cannot be rebuilt from this repository** (`100_build_manuscript.R`, `110_build_response_docs.R`): they read text sources that belong to the full analysis workspace.
6. **Bit-identical numeric output is not guaranteed** on a different R or package version. The main results were produced under R 4.4.1; the pre-submission sensitivity outputs under R 4.6.1.

**Bottom line:** this archive supports *verification and audit* of the reported results and *re-execution of the code* where inputs are available, but it is not a self-contained reproduction of the entire study.

---

## Licence

- **Code** (`scripts/`): MIT — see `LICENSE`.
- **Derived results and tables** (`results/`, `tables/`, documentation): CC BY 4.0 — see `LICENSE-DATA.md`.
- **Third-party data** are not redistributed here and remain under their original terms (TCGA/GDC data use terms, GEO, cBioPortal, UCSC Xena, MSigDB, Charoentong 2017 and MCP-counter licences). See `LICENSE-DATA.md`.

## Citation

If you use this code or these results, please cite the manuscript and this archive:

> Lu X, Sun X, Li J, Yang X, Tao T, Cao H. *Reduced GLP2R expression in colorectal cancer: a confounder-adjusted multi-cohort analysis of its prognostic, immune and epigenetic correlates*. Code and results archive: https://github.com/LXC27/GLP2R (DOI: https://doi.org/ZENODO_DOI_PLACEHOLDER).

The machine-readable citation is in `CITATION.cff`.

The DOI above identifies this archive (code plus derived results and tables); it does not replace the citations of the underlying data. Cite TCGA through the GDC/PanCanAtlas publications, GSE39582 and GSE17536 through their original GEO series publications, and cBioPortal, Thorsson et al., MSigDB, Charoentong et al. and MCP-counter through their respective references. The Zenodo concept DOI always resolves to the latest version; cite the version DOI when the exact release matters.
