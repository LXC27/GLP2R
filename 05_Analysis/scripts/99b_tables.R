## 99b_tables.R - rebuild Table 1, Table 2 and Supplementary Tables S1-S8 from the result ledgers
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages(library(writexl))
set.seed(SEED); options(scipen = 0); dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
M <- read.csv(file.path(RB, "05_Analysis", "results_master.csv"), stringsAsFactors = FALSE); R <- setNames(M$value, M$key)
rd <- function(f) read.csv(file.path(RES, f), stringsAsFactors = FALSE)
pp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
r3 <- function(x) ifelse(is.na(x), NA, signif(x, 4))
SHEETS <- list()   # sheet name -> data.frame
add_sheet <- function(name, df, csv = TRUE) { df[] <- lapply(df, function(z) if (is.numeric(z) && !is.integer(z)) signif(z, 4) else z); SHEETS[[name]] <<- df; if (csv) write_csv_utf8(df, file.path(TABD, paste0(name, ".csv"))); invisible(df) }

## ---------------------------------------------------------------- Table 1
ch <- rd("cohort_clinical_harmonised.csv"); vt <- readRDS(file.path(DC, "vst_tumor_591.rds")); so <- readRDS(file.path(RES, "survival_objects.rds"))
ch <- ch[ch$sample_type == "Primary Tumor", ]; stopifnot(nrow(ch) == 591, all(ch$sample %in% colnames(vt)))
ch$glp2r <- vt["GLP2R", ch$sample]; med <- median(ch$glp2r); ch$grp <- factor(ifelse(ch$glp2r >= med, "High", "Low"), c("High", "Low"))
for (v in c("sex", "stage_grp", "msi")) ch[[v]][ch[[v]] == ""] <- NA
ch$in_os <- ch$sample %in% so$dd$barcode; ch$os_ev <- so$dd$event[match(ch$sample, so$dd$barcode)]
ch$subtype <- c(C1 = "C1 Wound healing", C2 = "C2 IFN-gamma dominant", C3 = "C3 Inflammatory", C4 = "C4 Lymphocyte depleted", C6 = "C6 TGF-beta dominant")[sub(".*(C[0-9]).*", "\\1", ifelse(ch$immune_subtype == "", NA, ch$immune_subtype))]
nH <- sum(ch$grp == "High"); nL <- sum(ch$grp == "Low"); N <- nrow(ch)
cont_row <- function(lab, x) {
  ok <- !is.na(x); g <- ch$grp
  f <- function(z) sprintf("%.1f (%.1f-%.1f)", median(z), quantile(z, .25), quantile(z, .75))
  data.frame(Variable = lab, Level = "median (IQR)", High = f(x[ok & g == "High"]), Low = f(x[ok & g == "Low"]), All = f(x[ok]),
             P = pp(wilcox.test(x[ok] ~ g[ok], exact = FALSE)$p.value), Test = "Wilcoxon rank-sum",
             Missing = sprintf("%d / %d", sum(!ok & g == "High"), sum(!ok & g == "Low")), stringsAsFactors = FALSE)
}
cat_rows <- function(lab, x, levs = NULL, sim = FALSE) {
  x <- if (is.null(levs)) factor(x) else factor(x, levs); ok <- !is.na(x); g <- ch$grp
  tb <- table(x[ok], g[ok]); pv <- if (sim || any(chisq.test(tb)$expected < 5)) fisher.test(tb, simulate.p.value = nrow(tb) > 2 || ncol(tb) > 2, B = 1e5)$p.value else chisq.test(tb)$p.value
  if (nrow(tb) == 2 && ncol(tb) == 2) pv <- fisher.test(tb)$p.value
  out <- data.frame(Variable = c(lab, rep("", nrow(tb) - 1)), Level = rownames(tb),
                    High = sprintf("%d (%.1f%%)", tb[, "High"], 100 * tb[, "High"] / sum(tb[, "High"])),
                    Low = sprintf("%d (%.1f%%)", tb[, "Low"], 100 * tb[, "Low"] / sum(tb[, "Low"])),
                    All = sprintf("%d (%.1f%%)", rowSums(tb), 100 * rowSums(tb) / sum(tb)),
                    P = c(pp(pv), rep("", nrow(tb) - 1)), Test = c(ifelse(nrow(tb) == 2 && ncol(tb) == 2, "Fisher exact", ifelse(nrow(tb) > 2 && (sim || any(chisq.test(tb)$expected < 5)), "Fisher (simulated)", "Chi-square")), rep("", nrow(tb) - 1)),
                    Missing = c(sprintf("%d / %d", sum(!ok & g == "High"), sum(!ok & g == "Low")), rep("", nrow(tb) - 1)), stringsAsFactors = FALSE)
  out
}
t1 <- rbind(
  data.frame(Variable = "Patients (primary tumours)", Level = "n", High = as.character(nH), Low = as.character(nL), All = as.character(N), P = "", Test = "", Missing = ""),
  cont_row("Age at diagnosis, years", ch$age_y),
  cat_rows("Sex", ch$sex),
  cat_rows("AJCC stage", ch$stage_grp, c("I", "II", "III", "IV")),
  cat_rows("MSI status", ch$msi, c("MSS/MSI-L", "MSI-H")),
  cat_rows("Primary site", ch$site, c("COAD", "READ")),
  cont_row("Tumour mutation burden, mutations/Mb", ch$tmb),
  cont_row("ABSOLUTE tumour purity", ch$purity_absolute),
  cat_rows("Thorsson immune subtype", ch$subtype, sim = TRUE),
  {
    d <- ch[ch$in_os, ]; tb <- table(factor(d$os_ev, c(1, 0)), d$grp)
    data.frame(Variable = "Deaths in OS analysis set", Level = "n/N (%)", High = sprintf("%d/%d (%.1f%%)", tb[1, "High"], sum(tb[, "High"]), 100 * tb[1, "High"] / sum(tb[, "High"])),
               Low = sprintf("%d/%d (%.1f%%)", tb[1, "Low"], sum(tb[, "Low"]), 100 * tb[1, "Low"] / sum(tb[, "Low"])),
               All = sprintf("%d/%d (%.1f%%)", sum(tb[1, ]), sum(tb), 100 * sum(tb[1, ]) / sum(tb)), P = pp(fisher.test(tb)$p.value),
               Test = "Fisher exact (crude; not time-adjusted)", Missing = sprintf("%d / %d", sum(!ch$in_os & ch$grp == "High"), sum(!ch$in_os & ch$grp == "Low")))
  })
names(t1)[3:4] <- c(sprintf("GLP2R-high (n=%d)", nH), sprintf("GLP2R-low (n=%d)", nL)); names(t1)[5] <- sprintf("All (n=%d)", N); names(t1)[8] <- "Missing (high / low)"
add_sheet("Table1_cohort_characteristics", t1)
cat("Table 1\n"); print(t1, right = FALSE)
rm_write(data.frame(key = c("t1_median_vst_cut", "t1_n_high", "t1_n_low", "t1_age_p", "t1_sex_p", "t1_stage_p", "t1_msi_p", "t1_site_p", "t1_tmb_p", "t1_purity_p"),
                    value = c(round(med, 3), nH, nL, t1$P[t1$Variable == "Age at diagnosis, years"], t1$P[t1$Variable == "Sex"], t1$P[t1$Variable == "AJCC stage"], t1$P[t1$Variable == "MSI status"], t1$P[t1$Variable == "Primary site"],
                              t1$P[t1$Variable == "Tumour mutation burden, mutations/Mb"], t1$P[t1$Variable == "ABSOLUTE tumour purity"]), description = ""), "table1")

## ---------------------------------------------------------------- Table 2
mm <- so$mods; lab <- c(glp2r_z = "GLP2R expression (per 1 SD)", age10 = "Age (per 10 years)", sexmale = "Sex (male vs female)", stage_grpII = "Stage II vs I", stage_grpIII = "Stage III vs I",
                         stage_grpIV = "Stage IV vs I", siteREAD = "Primary site (rectum vs colon)", `msiMSI-H` = "MSI-H vs MSS/MSI-L")
tt <- lapply(c(M0 = 1, M1 = 3, M2 = 4), function(i) { d <- cox_tab(mm[[i]]); d$term <- lab[d$term]; d })
t2 <- data.frame(Variable = unname(lab), stringsAsFactors = FALSE)
for (k in names(tt)) { d <- tt[[k]]; i <- match(t2$Variable, d$term); t2[[paste0(k, "_HR (95% CI)")]] <- ifelse(is.na(i), "", fmt_hr(d$HR[i], d$lo[i], d$hi[i])); t2[[paste0(k, "_P")]] <- ifelse(is.na(i), "", pp(d$p[i])) }
ns <- sapply(c(1, 3, 4), function(i) c(mm[[i]]$n, mm[[i]]$nevent))
t2 <- rbind(t2, data.frame(Variable = c("Patients", "Deaths"), `M0_HR (95% CI)` = as.character(ns[1, c(1)]), M0_P = "", `M1_HR (95% CI)` = "", M1_P = "", `M2_HR (95% CI)` = "", M2_P = "", check.names = FALSE))
t2[nrow(t2) - 1, c(2, 4, 6)] <- as.character(ns[1, ]); t2[nrow(t2), c(2, 4, 6)] <- as.character(ns[2, ])
names(t2) <- c("Variable", "M0 HR (95% CI)", "M0 P", "M1 HR (95% CI)", "M1 P", "M2 HR (95% CI)", "M2 P")
add_sheet("Table2_cox_models", t2); cat("Table 2\n"); print(t2, right = FALSE)
t2long <- do.call(rbind, lapply(names(tt), function(k) cbind(model = k, n = mm[[c(M0 = 1, M1 = 3, M2 = 4)[[k]]]]$n, tt[[k]])))
add_sheet("Table2_cox_models_long", t2long, csv = TRUE)

## ---------------------------------------------------------------- S1 candidate genes
s1 <- rd("TableS1_candidate_genes_DE.csv")
names(s1) <- c("Gene", "baseMean", "log2FC_tumour_vs_normal", "lfcSE", "P_unpaired", "FDR_unpaired", "log2FC_paired", "FDR_paired", "Detection_status", "Mean_raw_count")
add_sheet("TableS1_candidate_gene_DE", s1)

## ---------------------------------------------------------------- S2 immune
s2a <- rd("immune_celltypes_charoentong.csv"); s2a <- s2a[order(s2a$FDR_partial_p), ]
names(s2a) <- c("Cell type", "rho (unadjusted)", "P (unadjusted)", "partial rho (purity, MSI, site)", "P (partial)", "n (partial)", "partial rho (RNA-based purity)", "P (RNA-based purity)", "Median ssGSEA, GLP2R-high", "Median ssGSEA, GLP2R-low",
                "AUC high vs low", "P (Wilcoxon, high vs low)", "FDR (unadjusted)", "FDR (partial)", "FDR (RNA-based purity)", "FDR (Wilcoxon)")
add_sheet("TableS2a_immune_cell_types", s2a)
s2b <- rd("immune_msi_stratified.csv"); add_sheet("TableS2b_immune_MSI_stratified", s2b)
s2c <- rd("immune_scores_mcp_estimate.csv"); add_sheet("TableS2c_MCPcounter_ESTIMATE", s2c)

## ---------------------------------------------------------------- S3 LASSO
s3a <- rd("lasso_feature_filter.csv"); names(s3a) <- c("Gene", "Fraction of tumours with >=10 reads", "In TCGA VST matrix", "Probe present on both arrays", "Retained for LASSO"); add_sheet("TableS3a_LASSO_feature_filter", s3a)
s3b <- rd("lasso_bootstrap_selection.csv"); names(s3b) <- c("Gene", "Bootstrap selection frequency (B=200)", "Coefficient at lambda.min (training set)", "Selected at lambda.min"); s3b[[3]] <- signif(s3b[[3]], 4); add_sheet("TableS3b_LASSO_selection", s3b)
sf <- rd("lasso_performance.csv")
hrtxt <- c(TCGA_train = R[["lasso_km_TCGA_train"]], TCGA_test = R[["lasso_km_TCGA_test"]], GSE39582 = R[["lasso_km_GSE39582"]], GSE17536 = R[["lasso_km_GSE17536"]])
adj <- c(TCGA_train = "not applicable (training set)", TCGA_test = R[["lasso_adjHR_TCGA_test"]], GSE39582 = R[["lasso_adjHR_GSE39582"]], GSE17536 = R[["lasso_adjHR_GSE17536"]])
cclin <- c(TCGA_train = "", TCGA_test = R[["lasso_Cindex_clin_TCGA_test"]], GSE39582 = R[["lasso_Cindex_clin_GSE39582"]], GSE17536 = R[["lasso_Cindex_clin_GSE17536"]])
s3c <- data.frame(Cohort = c("TCGA training", "TCGA internal test", "GSE39582", "GSE17536"), n = sf$n, Events = sf$events, `C-index` = sprintf("%.3f (%.3f-%.3f)", sf$C, sf$C_lo, sf$C_hi),
                  `AUC 12 / 36 / 60 months` = sprintf("%.3f / %.3f / %.3f", sf$AUC12, sf$AUC36, sf$AUC60), `High vs low risk (cohort-median split)` = unname(hrtxt),
                  `Adjusted HR per SD of risk score (age, sex, stage)` = unname(adj), `C-index of clinical model -> + risk score` = unname(cclin), check.names = FALSE)
add_sheet("TableS3c_LASSO_performance", s3c)
lo <- readRDS(file.path(RES, "lasso_objects.rds")); smn <- lo$sel_min[lo$sel_min != 0]; sm <- smn[abs(smn) > 1e-6]; sz <- smn[abs(smn) <= 1e-6]
gtxt <- paste0(paste(sprintf("%s(%.3f)", names(sm), sm), collapse = ";"), if (length(sz)) sprintf(" [%s coefficient %.1e, numerically zero, not counted]", names(sz), sz) else "")
s3d <- data.frame(Item = c("Split", "Steps computed before the split", "Steps using training data only", "Lambda.min", sprintf("Genes selected at lambda.min (|coefficient| > 1e-6; n = %d)", length(sm)),
                           "Genes selected at lambda.1se", "C-index of GLP2R alone, test set", "Bootstrap 95% CI of test-set C-index",
                           "Permutation reference (100 label-shuffled refits, each tuned by one 10-fold CV; descriptive): median C / one-sided P", "Share of permuted fits with an empty model"),
                  Value = c("566 patients with OS, event-stratified, approximately 70/30: training n=396 (80 events), test n=170 (34 events); seed 20240617",
                            "Detection filter (566 patients with OS) and VST (591 tumours)", "Gene standardisation (means, SDs), lambda, coefficients, risk-score median threshold",
                            R[["lasso_lambda_min"]], gtxt, R[["lasso_nsel_1se"]], R[["glp2r_alone_test_C"]], R[["lasso_test_C_boot_ci"]],
                            paste(R[["lasso_perm_null_median"]], "/", R[["lasso_perm_p"]]), R[["lasso_perm_prop_empty"]]))
add_sheet("TableS3d_LASSO_model_details", s3d)

## ---------------------------------------------------------------- S4 GEO
go <- readRDS(file.path(RES, "geo_objects.rds")); t39 <- go$t39; t17 <- go$t17
charac <- function(d, sexv, stagev, extra) {
  r <- list(c("Patients with overall survival", nrow(d)), c("Deaths", sum(d$event)), c("Follow-up, months, median (IQR; observed)", sprintf("%.1f (%.1f-%.1f)", median(d$time), quantile(d$time, .25), quantile(d$time, .75))),
            c("Age, years, mean (SD)", sprintf("%.1f (%.1f)", mean(d$age10 * 10, na.rm = TRUE), sd(d$age10 * 10, na.rm = TRUE))), c("Female, n (%)", sprintf("%d (%.1f%%)", sum(tolower(sexv) == "female", na.rm = TRUE), 100 * mean(tolower(sexv) == "female", na.rm = TRUE))))
  for (s in c("I", "II", "III", "IV")) r[[length(r) + 1]] <- c(paste0("Stage ", s, ", n (% of known)"), sprintf("%d (%.1f%%)", sum(stagev == s, na.rm = TRUE), 100 * sum(stagev == s, na.rm = TRUE) / sum(!is.na(stagev))))
  r[[length(r) + 1]] <- c("Stage missing, n", sum(is.na(stagev))); do.call(rbind, c(r, extra))
}
x39 <- charac(t39, t39$sex2, as.character(t39$stage), list(c("dMMR, n (% of known)", sprintf("%d (%.1f%%)", sum(t39$mmr == "dMMR", na.rm = TRUE), 100 * mean(t39$mmr == "dMMR", na.rm = TRUE))), c("dMMR unknown, n", sum(is.na(t39$mmr))),
                                                       c("Proximal location, n (%)", sprintf("%d (%.1f%%)", sum(t39$loc == "proximal"), 100 * mean(t39$loc == "proximal")))))
x17 <- charac(t17, t17$sex2, as.character(t17$stage), list(c("dMMR, n (% of known)", "not available"), c("dMMR unknown, n", "not available"), c("Proximal location, n (%)", "not available")))
s4a <- data.frame(Characteristic = x39[, 1], GSE39582 = x39[, 2], GSE17536 = x17[, 2], check.names = FALSE)
s4a <- rbind(s4a, data.frame(Characteristic = c("Platform", "GLP2R probeset", "GLP2R probeset mean log2 expression (percentile among all probesets)", "Non-tumoral mucosa samples (not in survival analyses)"),
                             GSE39582 = c("Affymetrix HG-U133 Plus 2.0", "221312_at", sprintf("%s (%s)", R[["GSE39582_glp2r_probe_mean_log2"]], R[["GSE39582_glp2r_probe_percentile"]]), "19"),
                             GSE17536 = c("Affymetrix HG-U133 Plus 2.0", "221312_at", sprintf("%s (%s)", R[["GSE17536_glp2r_probe_mean_log2"]], R[["GSE17536_glp2r_probe_percentile"]]), "0")))
add_sheet("TableS4a_GEO_cohort_characteristics", s4a)
gl <- function(nm, x, cohort, model, n, e) data.frame(Cohort = cohort, Model = model, n = n, Events = e, HR = x$HR, lo = x$lo, hi = x$hi, P = x$p)
s4b <- rbind(gl("a", go$r39$unadj, "GSE39582", "Unadjusted, per SD", go$r39$n0, go$r39$e0), gl("a", go$r39$adj, "GSE39582", "Adjusted for age, sex, stage, per SD", go$r39$n1, go$r39$e1), gl("a", go$r39$med, "GSE39582", "High vs low (median split)", go$r39$n0, go$r39$e0),
             gl("a", go$r17$unadj, "GSE17536", "Unadjusted, per SD", go$r17$n0, go$r17$e0), gl("a", go$r17$adj, "GSE17536", "Adjusted for age, sex, stage, per SD", go$r17$n1, go$r17$e1), gl("a", go$r17$med, "GSE17536", "High vs low (median split)", go$r17$n0, go$r17$e0))
s4b[, 5:7] <- lapply(s4b[, 5:7], function(z) round(z, 3)); s4b$P <- signif(s4b$P, 3)
me <- rd("glp2r_meta_cohorts.csv"); me <- me[, c("cohort", "n", "events", "HR", "lo", "hi")]; me[, 4:6] <- lapply(me[, 4:6], round, 3)
names(me) <- c("Cohort", "n", "Events", "Adjusted HR per SD", "Lower 95% CI", "Upper 95% CI")
me$Note <- c("", "", "", sprintf("Random-effects (DerSimonian-Laird); I2 = %s%%; P = %s; heterogeneity P(Q) = %s", R[["meta_I2_pct"]], round(go$p_re, 3), signif(go$p_Q, 2)))
add_sheet("TableS4b_GEO_cox", s4b); add_sheet("TableS4c_cross_cohort_synthesis", me)

## ---------------------------------------------------------------- S5 checkpoints
s5 <- rd("TableS5_checkpoints_recomputed.csv"); add_sheet("TableS5_checkpoints", s5)

## ---------------------------------------------------------------- S6 survival sensitivity
pf <- function(s) { g <- function(r) { z <- regmatches(s, regexec(r, s))[[1]]; if (length(z) > 1) z[2] else NA }
  data.frame(HR = as.numeric(g("^([0-9.]+) \\(")), lo = as.numeric(g("\\(([0-9.]+)-")), hi = as.numeric(g("-([0-9.]+)\\)")), p = as.numeric(g("p=([0-9.eE+-]+)")), n = as.numeric(g("n=([0-9]+)")), events = as.numeric(g("events=([0-9]+)"))) }
dd <- so$dd; dp <- dd[!is.na(dd$cdr_PFI_time) & dd$cdr_PFI_time > 0 & !is.na(dd$cdr_PFI), ]
dp$hi <- as.integer(dp$glp2r > median(dp$glp2r)); cp <- cox_tab(coxph(Surv(cdr_PFI_time / 30.4375, cdr_PFI) ~ hi, dp))
rm_write(data.frame(key = c("pfi_median_split_hr", "pfi_median_split_n"), value = c(sprintf("%.2f (%.2f-%.2f); p=%.3f", cp$HR, cp$lo, cp$hi, cp$p), nrow(dp)), description = "PFI, GLP2R high vs low (cohort median)"), "table_pfi")
sens <- rbind(
  cbind(Analysis = "M0: unadjusted, all patients with OS", pf(R[["cox_glp2r_hr_M0_unadjusted_all_OS_"]])),
  cbind(Analysis = "M0: unadjusted, restricted to the M2 sample", pf(R[["cox_glp2r_hr_M0_unadjusted_same_n_as_M2_"]])),
  cbind(Analysis = "M1: age, sex, stage, site", pf(R[["cox_glp2r_hr_M1_age_sex_stage_site"]])),
  cbind(Analysis = "M1 refitted on the M2 sample (no MSI)", pf(R[["cox_glp2r_hr_M1_on_M2_sample_no_MSI_"]])),
  cbind(Analysis = "M2: + MSI (primary model)", pf(R[["cox_glp2r_hr_M2_MSI_primary_model_"]])),
  cbind(Analysis = "M2, follow-up truncated at 60 months", pf(R[["cox_m2_trunc60_glp2r"]])),
  cbind(Analysis = "M2, stratified by stage", pf(R[["cox_m2_stage_strata_glp2r"]])),
  cbind(Analysis = "M2, log2CPM instead of VST", pf(R[["cox_m2_log2cpm_glp2r"]])),
  cbind(Analysis = "M2, MSS/MSI-L tumours only", pf(R[["cox_glp2r_in_MSS"]])),
  cbind(Analysis = "Multiple imputation of stage and MSI (m=30), all patients with OS", pf(R[["mi_glp2r_hr"]])),
  cbind(Analysis = "M2, GLP2R high vs low (cohort median)", pf(R[["cox_m2_median_binary_glp2r"]])),
  cbind(Analysis = "Progression-free interval, unadjusted", pf(R[["pfi_unadj_glp2r"]])),
  cbind(Analysis = "Progression-free interval, adjusted (M2 covariates)", pf(R[["pfi_adj_glp2r"]])))
cc5 <- so$cc5; ev60 <- sum(cc5$event == 1 & cc5$time <= 60)
setn <- function(a, n, e) { i <- sens$Analysis == a; sens$n[i] <<- n; sens$events[i] <<- e }
setn("M0: unadjusted, all patients with OS", 566, 114); setn("M0: unadjusted, restricted to the M2 sample", 451, 77); setn("M1: age, sex, stage, site", 459, 77); setn("M1 refitted on the M2 sample (no MSI)", 451, 77)
setn("M2: + MSI (primary model)", 451, 77); setn("M2, follow-up truncated at 60 months", 451, ev60); setn("M2, stratified by stage", 451, 77); setn("M2, log2CPM instead of VST", 451, 77)
setn("Multiple imputation of stage and MSI (m=30), all patients with OS", 566, 114); setn("M2, GLP2R high vs low (cohort median)", 451, 77)
setn("Progression-free interval, unadjusted", 564, 142); setn("Progression-free interval, adjusted (M2 covariates)", 451, NA)
names(sens) <- c("Analysis", "HR", "Lower 95% CI", "Upper 95% CI", "P", "n", "Events")
add_sheet("TableS6a_survival_sensitivity", sens)
s6b <- rd("cox_subgroups.csv"); s6b[, c("HR", "lo", "hi")] <- lapply(s6b[, c("HR", "lo", "hi")], round, 3); add_sheet("TableS6b_survival_subgroups", s6b)
s6c <- rd("glp2r_timeROC.csv"); names(s6c) <- c("Time (months)", "AUC", "SE", "Lower 95% CI", "Upper 95% CI"); add_sheet("TableS6c_GLP2R_timeROC", s6c)
s6d <- data.frame(Item = c("Median follow-up (reverse Kaplan-Meier), months", "Patients at risk at 12 / 36 / 60 / 96 months", "Median split, high vs low GLP2R", "Median split log-rank P", "Maxstat-corrected P (cut at the 52.1th percentile; exploratory)", "Tertile log-rank P",
                          "Proportional-hazards test, global P / GLP2R P", "Non-linearity (spline vs linear) P", "Interaction GLP2R x MSI P", "Harrell C, clinical model -> + GLP2R", "Delta C (bootstrap 95% CI)", "Likelihood-ratio test for adding GLP2R, P"),
                  Value = c(R[["followup_median_months_revKM"]], R[["n_at_risk_12_36_60_96"]], R[["km_median_hr_high_vs_low"]], R[["km_median_logrank_p"]], R[["maxstat_corrected_p"]], R[["km_tertile_logrank_p"]], paste(R[["zph_global_p"]], "/", R[["zph_glp2r_p"]]), R[["glp2r_nonlinearity_p"]],
                            R[["interaction_glp2r_msi_p"]], paste(R[["cindex_base"]], "->", R[["cindex_with_glp2r"]]), R[["cindex_delta_ci"]], R[["lrt_glp2r_added_p"]]))
add_sheet("TableS6d_survival_diagnostics", s6d)

## ---------------------------------------------------------------- S7 epigenetic / genomic
m1 <- rd("glp2r_methylation_tumour_vs_normal.csv"); m2 <- rd("glp2r_meth_vs_expr_probes.csv")
s7a <- merge(m1, m2, by = "probe", suffixes = c("_tn", "_expr")); s7a <- s7a[order(s7a$pos), ]; add_sheet("TableS7a_GLP2R_promoter_methylation", s7a)
s7b <- rd("dnmt_tet_regulators.csv"); add_sheet("TableS7b_DNMT_TET", s7b)
s7c <- rd("mutation_frequency.csv"); add_sheet("TableS7c_mutation_frequency", s7c)
s7d <- rd("glp2r_mutations.csv"); add_sheet("TableS7d_GLP2R_mutations", s7d)
s7e <- rd("mutation_vs_glp2r.csv"); add_sheet("TableS7e_mutation_vs_GLP2R", s7e)
cn <- strsplit(R[["cnv_table"]], "; ")[[1]]; cnt <- data.frame(`GISTIC state` = sub("=.*", "", cn), n = as.integer(sub(".*=", "", cn)), check.names = FALSE)
s7f <- rbind(cnt, data.frame(`GISTIC state` = c("Loss (shallow or deep) %", "Gain %", "Deep deletion %", "Spearman rho, copy-number state vs expression", "Kruskal-Wallis P", "Loss vs non-loss P (Wilcoxon)", "Adjusted rho (purity, MSI, site)"),
                             n = NA, check.names = FALSE))
s7f$Value <- c(rep("", nrow(cnt)), R[["cnv_loss_pct"]], R[["cnv_gain_pct"]], R[["cnv_deep_del_pct"]], R[["cnv_spearman_rho"]], R[["cnv_kw_p"]], R[["cnv_loss_vs_nonloss_p"]], R[["cnv_rho_adj"]])
add_sheet("TableS7f_GLP2R_copy_number", s7f)
ci <- rd("cimp_associations.csv")
s7g <- data.frame(Cohort = ci$cohort, Variable = ci$variable, Contrast = ci$contrast, n_case = ci$n_case, n_reference = ci$n_ref,
                  Median_case = ci$median_case, Median_reference = ci$median_ref, AUC = ci$auc_case_gt_ref, P_Wilcoxon = ci$p_wilcoxon, FDR = ci$p_fdr,
                  Stratum = ci$stratum, n_stratum = ci$n_stratum, P_stratum = ci$p_wilcoxon_stratum, Beta_adjusted = ci$beta_adjusted, P_adjusted = ci$p_adjusted)
add_sheet("TableS7g_CIMP_BRAF_MMR_associations", s7g)

## ---------------------------------------------------------------- S8 pathway, co-expression, compartments, WGCNA
fg <- rd("fgsea_all.csv"); fg <- fg[fg$ranking == "adjusted", ]; fg$padj <- fg$FDR_collection
s8a <- fg[fg$collection == "Hallmark", c("pathway", "size", "ES", "NES", "pval", "padj")]; s8a <- s8a[order(s8a$padj), ]; names(s8a)[6] <- "FDR"; add_sheet("TableS8a_GSEA_Hallmark", s8a)
fk <- fg[fg$collection == "KEGG" & fg$padj < 0.05, ]; fk <- fk[order(-abs(fk$NES)), ]
s8b <- rbind(head(fk[fk$NES > 0, ], 25), head(fk[fk$NES < 0, ], 25))[, c("pathway", "size", "ES", "NES", "pval", "padj")]; names(s8b)[6] <- "FDR"; add_sheet("TableS8b_GSEA_KEGG_top", s8b)
cg <- rd("glp2r_correlated_genes.csv"); cg <- cg[order(-cg$rho_partial), ]
s8c <- rbind(cbind(Direction = "Positive", head(cg, 50)), cbind(Direction = "Negative", head(cg[order(cg$rho_partial), ], 50))); names(s8c) <- c("Direction", "Gene", "partial rho", "raw rho", "P (partial)", "FDR (partial)"); add_sheet("TableS8c_top_correlated_genes", s8c)
s8d <- rd("candidate_genes_correlation.csv"); add_sheet("TableS8d_candidate_gene_correlation", s8d)
s8e <- rd("glp2r_compartment_markers.csv"); add_sheet("TableS8e_compartment_markers", s8e)
s8f <- rd("wgcna_module_table.csv"); add_sheet("TableS8f_WGCNA_modules", s8f)
en <- rd("wgcna_module_enrichment.csv"); en <- en[en$module == "turquoise", ]; en <- en[order(en$collection, en$FDR), ]
s8g <- do.call(rbind, lapply(split(en, en$collection), head, 10)); rownames(s8g) <- NULL; add_sheet("TableS8g_WGCNA_turquoise_enrichment", s8g)
wo <- readRDS(file.path(RES, "wgcna_objects.rds")); s8h <- data.frame(Gene = names(wo$hub), kME_turquoise = round(as.numeric(wo$hub), 3)); add_sheet("TableS8h_WGCNA_turquoise_hubs", s8h)

## ---------------------------------------------------------------- workbook
nm <- names(SHEETS); nm2 <- substr(sub("^Table", "", nm), 1, 31)
write_xlsx(setNames(SHEETS, nm2), file.path(TABD, "Tables_Table1-2_and_S1-S8.xlsx"))
cat("\nSheets written:", length(SHEETS), "\n"); print(names(SHEETS))
