# Re-running the analysis

This document explains how to reconstruct the analysis workspace that the archived scripts expect, how to supply the input data, and in what order to run the scripts.

The scripts in `scripts/` are **portable**: the hard-coded `W:/Re_build`, `W:/Supplementary` and `V:/` drive paths of the original workspace have been replaced by environment variables with the original values as defaults, and the package library is selected by R version. You can either set the environment variables, or recreate the `W:`/`V:` drive mapping (which the defaults still assume). **Do not edit the scripts if you want to match the archived files.**

The archive does not contain the data cache, the R library or the author-supplied `V:` objects. See the section "What can and cannot be reproduced from this archive" in `README.md` for the full list.

---

## 0. Path configuration (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `GLP2R_RB` | `W:/Re_build` | Root of the analysis workspace (the folder that contains `_work/`, `05_Analysis/`, `03_Tables/`, …). |
| `GLP2R_SCRIPTS` | `<GLP2R_RB>/05_Analysis/scripts` | Where the `.R` scripts live. |
| `GLP2R_OUT` | `<GLP2R_RB>` | Output root (manuscript, response, claims outputs). |
| `GLP2R_INPUT` | `W:/Supplementary` | Folder holding the author-supplied gene-list CSV (`Table_S1_GLP_genes.csv`). |
| `GLP2R_ORIG` | `V:/` | Folder holding the author-supplied `.rds` objects (e.g. `TCGA_CRC_combined.rds`). |

On Windows PowerShell, set them for the current session with:

```powershell
$env:GLP2R_RB     = "D:\GLP2R_rebuild\Re_build"
$env:GLP2R_INPUT  = "D:\GLP2R_rebuild\Supplementary"
$env:GLP2R_ORIG   = "D:\GLP2R_inputs\Rdata"
$env:GLP2R_SCRIPTS = "$env:GLP2R_RB\05_Analysis\scripts"
$env:GLP2R_OUT    = $env:GLP2R_RB
```

If you prefer not to set variables, the defaults expect the original layout, which you can create with `subst` (see §1). Use **ASCII-only** paths — R and `data.table` handle non-ASCII inconsistently on Windows.

---

## 1. Recreate the workspace layout

The scripts assume a workspace with this layout (whether under `W:/Re_build` or a custom `GLP2R_RB`):

```
<GLP2R_RB>/
├── _work/Rlib/                  <- project-local R library (R < 4.5; created by 00_setup.R)
├── 01_Manuscript/src/           <- manuscript text sources (optional; only for step 100)
├── 02_Figures/
│   ├── Main/                    <- Figure1 … Figure7 (created by the scripts)
│   └── Supplementary/           <- FigureS1 … FigureS6 (created by 98_figS.R)
├── 03_Tables/                   <- table CSVs + xlsx (created by 99b_tables.R)
├── 04_Response/src/             <- revision text sources (optional; only for 110/120)
├── 05_Analysis/
│   ├── data_cache/              <- input data + derived objects (you must supply)
│   ├── logs/                    <- console logs (created by the scripts)
│   ├── results/                 <- result CSVs and .rds objects (created by the scripts)
│   └── scripts/                 <- this repository's scripts/
├── 06_References/               <- reference CSVs (optional; only for 99d/100)
└── README.md
<GLP2R_INPUT>/                   <- author-supplied gene-list CSVs (not in this archive)
```

`00_common.R` creates `results/`, `logs/`, `02_Figures/Main`, `02_Figures/Supplementary` and `03_Tables` on first run, but **not** `data_cache/` — create it yourself.

To use the default `W:`/`V:` mapping (the original layout), from a non-elevated Command Prompt:

```bat
subst W: D:\GLP2R_rebuild          :: W:\Re_build must exist
subst V: D:\GLP2R_inputs\Rdata     :: V:\TCGA_CRC_combined.rds must exist
```

`subst` is per-user and disappears on reboot; re-run after a restart.

---

## 2. Install R and the packages

Two supported environments:

- **R 4.4.1** — reproduces the main analysis bit-for-bit (closest). The setup scripts install binary packages into `<GLP2R_RB>/_work/Rlib`.
- **R ≥ 4.5** — the portable scripts load packages from the user library (`R_LIBS_USER`) instead of the R 4.4.1-compiled `_work/Rlib` (which would fail to load under a newer R). Install the packages into the user library first.

**Minimal install (figures + sensitivity only), any R ≥ 4.5:**

```powershell
Rscript.exe <GLP2R_SCRIPTS>\S00_install_r46.R
```

`S00_install_r46.R` installs the CRAN packages needed by the figure and sensitivity scripts (ggplot2, dplyr, tidyr, tibble, patchwork, ggrepel, cowplot, ragg, scales, glmnet, survival) from the TUNA (Tsinghua) and cloud.r-project.org mirrors.

**Full install (the whole pipeline), matching the reported run:**

```bat
set RSCRIPT="C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
%RSCRIPT% W:/Re_build/05_Analysis/scripts/00_setup.R
%RSCRIPT% W:/Re_build/05_Analysis/scripts/00b_install_extra.R
%RSCRIPT% W:/Re_build/05_Analysis/scripts/00c_install_annotation.R
%RSCRIPT% W:/Re_build/05_Analysis/scripts/00d_install_dose.R
%RSCRIPT% W:/Re_build/05_Analysis/scripts/00e_check_pkgs.R
```

`00_setup.R` installs the main CRAN and Bioconductor sets (including BiocManager); `00b`–`00d` install the additional packages used by individual steps. Under R ≥ 4.5, the install scripts automatically target `R_LIBS_USER` instead of `_work/Rlib`, so the full set is installed into your user library. `00e_check_pkgs.R` prints what is present and what is missing.

Additional requirements outside the package library:

- **`gzip`** on `PATH` — `02_extract_external.R` streams the GISTIC files through `gzip -dc` (a pure-R fallback exists).
- **`curl` / `curl.exe`** on `PATH` — `22a_gdc_paired_download.R` (GDC), `62_epi_extra.R` (cBioPortal), `99d_crossref_refs.R` (Crossref). Windows 10/11 ships `curl.exe`.
- The manuscript-finishing PowerShell helper `100b_word_finish.ps1` is **not** in this release (it requires Microsoft Word COM automation and is not needed for any analysis step).

---

## 3. Supply the input data

All inputs go into `<GLP2R_RB>/05_Analysis/data_cache/`. `results/data_manifest.csv` lists every expected file with its size, modification date and MD5 checksum — use it as the checklist. The full source table is in `README.md` → "Data used in the study"; in summary:

**Public downloads (manual):** Xena `COAD_*/READ_*` clinical/GISTIC/HM450 matrices + `illuminaMethyl450_probeMap.tsv`; `GSE39582_series_matrix.txt.gz`; `GSE17536_series_matrix.txt.gz`; `GPL570_table.txt`; `Thorsson_2018_immune_landscape.xlsx` + `Thorsson_Subtype_Immune_Model_Based.txt.gz`; `ABSOLUTE_mastercalls.txt`; `Charoentong2017.gmt`; `MCPcounter_genes.txt` + `MCPcounter_probesets.txt`; cBioPortal `cbio_pancan_sample.json` / `cbio_pancan_patient.json`.

**Auto-downloaded:** `cbio_mutations_selected.json` (by `62_epi_extra.R`), the 32 GDC paired count files (by `22a_gdc_paired_download.R`).

**Author-supplied (not in the archive):** `V:/TCGA_CRC_combined.rds` (primary input), `W:/Supplementary/Table_S1_GLP_genes.csv` and `06_References/reference_verification.csv` (read by `99d_crossref_refs.R`).

If you only want to audit the reported numbers, you do not need these files: `results/results_master.csv` and `tables/` contain the reported outputs.

---

## 4. Run the scripts

Replace `%RSCRIPT%` with the full path to `Rscript.exe`. The recommended full sequence reflects the actual code dependencies (see `README.md` → "Script execution order"):

```bat
:: ---- stage 0: environment (once) -------------------------------------
%RSCRIPT% 00_setup.R
%RSCRIPT% 00b_install_extra.R
%RSCRIPT% 00c_install_annotation.R
%RSCRIPT% 00d_install_dose.R
%RSCRIPT% 00e_check_pkgs.R

:: ---- stage 1: external data caches ------------------------------------
%RSCRIPT% 02_extract_external.R
%RSCRIPT% 02b_extract_cbio_thorsson.R
%RSCRIPT% 03_geo_prepare.R

:: ---- stage 2: cohort, expression, paired analysis ----------------------
%RSCRIPT% 10_cohort.R
%RSCRIPT% 20_deseq.R
%RSCRIPT% 22a_gdc_paired_download.R
%RSCRIPT% 22b_paired_deseq.R
%RSCRIPT% 21_expr_tests.R

:: ---- stage 3: downstream analyses (order within the stage is flexible) -
%RSCRIPT% 30_survival.R
%RSCRIPT% 32_geo_glp2r.R
%RSCRIPT% 40_immune.R
%RSCRIPT% 50_pathway.R
%RSCRIPT% 60_epi_cnv.R
%RSCRIPT% 62_epi_extra.R
%RSCRIPT% 63_cimp.R
%RSCRIPT% 70_lasso.R
%RSCRIPT% 80_wgcna.R

:: ---- stage 4: figures (must precede the first 99a, for rm_figS.csv) ---
%RSCRIPT% 91_fig1.R
%RSCRIPT% 92_fig2.R
%RSCRIPT% 93_fig3.R
%RSCRIPT% 94_fig4.R
%RSCRIPT% 95_fig5.R
%RSCRIPT% 96_fig6.R
%RSCRIPT% 97_fig7.R
%RSCRIPT% 98_figS.R

:: ---- stage 5: result ledger and tables (99a twice) --------------------
%RSCRIPT% 99a_results_master.R
%RSCRIPT% 99b_tables.R
%RSCRIPT% 99a_results_master.R

:: ---- stage 6: audit outputs (optional) ---------------------------------
%RSCRIPT% 99d_crossref_refs.R

:: ---- stage 7: pre-submission sensitivity (after the main pipeline) -----
%RSCRIPT% S01_presubmission_sensitivity.R
```

Why the order deviates from plain numeric order:

1. **Install first.** `00_setup.R` creates `_work/Rlib`; running an analysis script before the library exists fails at `library()`.
2. **20 → 22a → 22b → 21.** `21_expr_tests.R` reads `results/paired_glp2r_vst.csv`, which `22b_paired_deseq.R` writes after reading the GDC download manifest produced by `22a`.
3. **99a → 99b → 99a.** `99b_tables.R` reads `results_master.csv` but also writes `rm_table1.csv`/`rm_table_pfi.csv`, which `99a_results_master.R` aggregates. Run 99a again after 99b.
4. `98_figS.R` writes `rm_figS.csv`, so the figure stage must precede the first `99a`.
5. `S01_presubmission_sensitivity.R` reads the main result objects (`cox_*.rds`, immune/wgcna/epi objects) and must run after the main pipeline.

Optional manuscript/revision steps need text sources not in this repository — do not attempt unless you have the full analysis workspace.

---

## 5. Check the results

After a successful run, compare against this archive:

```powershell
Compare-Object (Import-Csv W:\Re_build\05_Analysis\results_master.csv).key `
               (Import-Csv results\results_master.csv).key |
  Format-Table -AutoSize
```

Exact equality is expected only with the same R and package versions and the same input files (MD5s in `data_manifest.csv`). Small numeric differences are acceptable evidence of a version effect; structural differences (missing `key`s, different sample sizes) indicate an input problem.

Approximate reference values to spot-check: `expr` module → `glp2r_log2fc` ≈ −4.073, `n_pairs` = 32; `survival` → `m2_glp2r_hr` ≈ "0.86 (0.69-1.07)"; `cimp` → `tcga_cimp_n` = 207. The pre-submission sensitivity outputs are in `results/S01_presubmission_sensitivity.txt` (β ≈ −4.04 in block A, HK CI 0.559–1.310 in block C, MLH1 promoter ρ ≈ −0.43 in block F).

---

## 6. Runtime and resources

- The archived production run (logs from 2026-10-11) spans 02:46–04:55 for the main pipeline (stages 1–5), i.e. **about 2 h 10 min**. A fresh run additionally needs package installation (tens of minutes to hours) and manual downloads of the multi-hundred-MB cache (the Xena HM450 files alone are ≈ 400 MB; GEO matrices ≈ 220 MB).
- RAM: 16 GB is sufficient. Disk: ≈ 3–4 GB for cache + R library + outputs.

---

## 7. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `cannot open file 'W:/Re_build/...` or `cannot open the connection` | The path variable or drive mapping is wrong. Set `GLP2R_RB` (and `GLP2R_INPUT`/`GLP2R_ORIG`) to the correct folders, or `subst W:` so that `W:\Re_build` exists. |
| `there is no package called 'X'` | Package library incomplete. Re-run `S00_install_r46.R` (minimal) or `00_setup.R`/`00b`–`00d` (full), then `00e_check_pkgs.R`. |
| `rlang.dll` / `SETLENGTH` load error under R ≥ 4.5 | A script still points at the R 4.4.1 `_work/Rlib`. The portable scripts branch on `getRversion()`; make sure you are running the scripts from this release, and that `R_LIBS_USER` exists and contains the installed packages. |
| `cannot open file '.../data_cache/xxx'` | The input is missing from the cache. Look up its source and MD5 in `results/data_manifest.csv` and place it under the exact name. |
| `gzip` errors from `fread` | `gzip` is not on `PATH`; the fallback in `02_extract_external.R` is used only for GISTIC, so install `gzip` to be safe. |
| `curl` failures / TLS errors | Network/proxy issue; retry. `22a` and `62` skip existing downloads, so re-running resumes. |
| `00_setup.R` install failures (SSL or 404) | Transient CRAN/Bioc mirror issue; re-run the script. |
| `Error: cannot allocate vector of size ...` | Insufficient RAM or non-64-bit R; use 64-bit R with ≥16 GB. |
| `10_cohort.R` fails on `V:/TCGA_CRC_combined.rds` | The author-assembled primary matrix is not in the archive; set `GLP2R_ORIG` to the folder holding it, or see README "What can and cannot be reproduced". |
| Numbers differ slightly from the released ledger | Different R/package versions. Match `sessionInfo.txt` (R 4.4.1) for exact reproduction. |
