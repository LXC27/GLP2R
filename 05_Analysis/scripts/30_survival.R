## 30_survival.R : GLP2R and overall survival in TCGA-COAD/READ (pre-specified, covariate-adjusted, same-n comparisons)
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressPackageStartupMessages({library(maxstat); library(timeROC); library(mice)})
sink(file.path(LOGD, "30_survival.txt"), split = TRUE)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d
co <- readRDS(file.path(DC, "cohort.rds")); cl <- co$clin
RM <- list(); add <- function(k, v, d = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = d)
d$time <- d$os_time_m; d$event <- d$os_event
dd <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
dd$glp2r_z <- as.numeric(scale(dd$glp2r)); dd$age10 <- dd$age_y / 10
dd$grp <- factor(ifelse(dd$glp2r > median(dd$glp2r), "High", "Low"), levels = c("Low", "High"))
cat("OS set n =", nrow(dd), " events =", sum(dd$event), " median OS-follow-up months (reverse KM):\n")
rk <- survfit(Surv(time, 1 - event) ~ 1, dd); print(summary(rk)$table[c("records", "median", "0.95LCL", "0.95UCL")])
add("os_n", nrow(dd)); add("os_events", sum(dd$event)); add("followup_median_months_revKM", round(summary(rk)$table["median"], 1)); add("followup_max_months", round(max(dd$time), 1))
nar <- summary(survfit(Surv(time, event) ~ 1, dd), times = c(12, 36, 60, 96))$n.risk; cat("number at risk 12/36/60/96 months:", nar, "\n"); add("n_at_risk_12_36_60_96", paste(nar, collapse = "/"))
## ---- 1. KM / log-rank (median split, pre-specified) + maxstat-adjusted p (exploratory)
lr <- survdiff(Surv(time, event) ~ grp, dd); p_lr <- 1 - pchisq(lr$chisq, 1)
cx0 <- coxph(Surv(time, event) ~ grp, dd); t0 <- cox_tab(cx0)[1, ]; cat(sprintf("Median split: log-rank p=%s; HR(high vs low)=%s\n", fmt_p(p_lr), fmt_hr(t0$HR, t0$lo, t0$hi)))
add("km_median_logrank_p", signif(p_lr, 4)); add("km_median_hr_high_vs_low", fmt_hr(t0$HR, t0$lo, t0$hi))
ms <- maxstat.test(Surv(time, event) ~ glp2r, data = dd, smethod = "LogRank", pmethod = "Lau92", minprop = 0.1, maxprop = 0.9)
cat(sprintf("maxstat: cutpoint=%.3f (percentile %.1f%%) corrected p=%s\n", ms$estimate, 100 * mean(dd$glp2r <= ms$estimate), fmt_p(ms$p.value)))
add("maxstat_cut_percentile", round(100 * mean(dd$glp2r <= ms$estimate), 1)); add("maxstat_corrected_p", signif(ms$p.value, 4), "Lausen-Schumacher adjusted; exploratory")
terc <- cut(dd$glp2r, quantile(dd$glp2r, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("T1", "T2", "T3")); lt <- survdiff(Surv(time, event) ~ terc, dd); add("km_tertile_logrank_p", signif(1 - pchisq(lt$chisq, 2), 4))
## ---- 2. Cox models (continuous GLP2R per 1 SD)
f_base <- c("age10", "sex", "stage_grp", "site"); fm <- function(v) as.formula(paste("Surv(time, event) ~", paste(v, collapse = " + ")))
cc4 <- dd[complete.cases(dd[, c(f_base)]), ]; cc5 <- dd[complete.cases(dd[, c(f_base, "msi")]), ]
cat("\ncomplete-case n: 4-var (age,sex,stage,site) =", nrow(cc4), "; +MSI =", nrow(cc5), "\n"); add("cox_n_cc_no_msi", nrow(cc4)); add("cox_n_cc_msi", nrow(cc5)); add("cox_events_cc_msi", sum(cc5$event))
mods <- list(
  `M0 unadjusted (all OS)`      = coxph(fm("glp2r_z"), dd),
  `M0 unadjusted (same n as M2)`= coxph(fm("glp2r_z"), cc5),
  `M1 + age, sex, stage, site`  = coxph(fm(c("glp2r_z", f_base)), cc4),
  `M2 + MSI (primary model)`    = coxph(fm(c("glp2r_z", f_base, "msi")), cc5),
  `M2 without GLP2R (same n)`   = coxph(fm(c(f_base, "msi")), cc5),
  `M1 on M2 sample (no MSI)`    = coxph(fm(c("glp2r_z", f_base)), cc5))
tabm <- do.call(rbind, lapply(names(mods), function(n) { t <- cox_tab(mods[[n]]); t <- t[t$term == "glp2r_z", ]; if (nrow(t)) data.frame(model = n, n = mods[[n]]$n, events = mods[[n]]$nevent, t[, c("HR", "lo", "hi", "p")]) }))
print(tabm); write_csv_utf8(tabm, file.path(RES, "cox_glp2r_models.csv"))
for (i in seq_len(nrow(tabm))) add(paste0("cox_glp2r_hr_", gsub("[^A-Za-z0-9]+", "_", tabm$model[i])), sprintf("%s; p=%s; n=%d", fmt_hr(tabm$HR[i], tabm$lo[i], tabm$hi[i]), fmt_p(tabm$p[i]), tabm$n[i]))
m2 <- mods[["M2 + MSI (primary model)"]]; full_tab <- cox_tab(m2); cat("\nPrimary model M2 coefficient table:\n"); print(full_tab, digits = 3)
saveRDS(list(full_tab = full_tab, tabm = tabm), file.path(RES, "cox_primary.rds"))
## assumption checks
zp <- cox.zph(m2); cat("\ncox.zph (M2):\n"); print(zp$table); add("zph_global_p", signif(zp$table["GLOBAL", "p"], 4)); add("zph_glp2r_p", signif(zp$table["glp2r_z", "p"], 4))
pm <- coxph(Surv(time, event) ~ splines::ns(glp2r_z, df = 3) + age10 + sex + stage_grp + site + msi, cc5); lrt_nl <- anova(m2, pm); pnl <- lrt_nl[2, "Pr(>|Chi|)"]
cat("natural-spline (3 df) non-linearity LRT p =", fmt_p(pnl), "\n"); add("glp2r_nonlinearity_p", signif(pnl, 4), "LRT: ns(df=3) vs linear")
## ---- 3. Stage-stratified sensitivity, 5-year truncation, log2CPM measure
cc5t <- cc5; cc5t$event[cc5t$time > 60] <- 0; cc5t$time <- pmin(cc5t$time, 60); m_t <- coxph(fm(c("glp2r_z", f_base, "msi")), cc5t); t <- cox_tab(m_t)[1, ]; add("cox_m2_trunc60_glp2r", sprintf("%s; p=%s", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p)))
m_s <- coxph(Surv(time, event) ~ glp2r_z + age10 + sex + strata(stage_grp) + site + msi, cc5); t <- cox_tab(m_s)[1, ]; add("cox_m2_stage_strata_glp2r", sprintf("%s; p=%s", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p)))
cc5$glp2r_lcpm_z <- as.numeric(scale(cc5$glp2r_lcpm)); m_l <- coxph(fm(c("glp2r_lcpm_z", f_base, "msi")), cc5); t <- cox_tab(m_l)[1, ]; add("cox_m2_log2cpm_glp2r", sprintf("%s; p=%s", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p)))
cc5$grp2 <- factor(ifelse(cc5$glp2r > median(cc5$glp2r), "High", "Low"), levels = c("Low", "High")); m_g <- coxph(fm(c("grp2", f_base, "msi")), cc5); t <- cox_tab(m_g)[1, ]; add("cox_m2_median_binary_glp2r", sprintf("%s; p=%s", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p)))
## ---- 4. interaction with MSI and stage
mi <- coxph(fm(c("glp2r_z * msi", "age10", "sex", "stage_grp", "site")), cc5); lrt <- anova(m2, mi); add("interaction_glp2r_msi_p", signif(lrt[2, "Pr(>|Chi|)"], 4))
cat("GLP2R x MSI interaction LRT p =", fmt_p(lrt[2, "Pr(>|Chi|)"]), "\n")
mss <- cc5[cc5$msi == "MSS/MSI-L", ]; mh <- cc5[cc5$msi == "MSI-H", ]
t <- cox_tab(coxph(fm(c("glp2r_z", "age10", "sex", "stage_grp", "site")), mss))[1, ]; add("cox_glp2r_in_MSS", sprintf("%s; p=%s; n=%d events=%d", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p), nrow(mss), sum(mss$event)))
cat("MSS/MSI-L only:", sprintf("%s p=%s n=%d ev=%d\n", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p), nrow(mss), sum(mss$event))); add("n_msih_events", sum(mh$event)); add("n_msih_in_cox", nrow(mh))
## ---- 5. Subgroup forest (adjusted), P-interaction
sub <- list(Age = list(var = "agegrp", lv = c("<65 y", ">=65 y")), Sex = list(var = "sex", lv = c("female", "male")), Stage = list(var = "stagegrp", lv = c("I-II", "III-IV")), Site = list(var = "site", lv = c("COAD", "READ")), MSI = list(var = "msi", lv = c("MSS/MSI-L", "MSI-H")))
cc5$agegrp <- factor(ifelse(cc5$age_y >= 65, ">=65 y", "<65 y"), levels = c("<65 y", ">=65 y")); cc5$stagegrp <- factor(ifelse(cc5$stage_grp %in% c("I", "II"), "I-II", "III-IV"), levels = c("I-II", "III-IV"))
adjv <- function(excl) setdiff(c("age10", "sex", "stage_grp", "site", "msi"), excl)
excl_map <- list(Age = "age10", Sex = "sex", Stage = "stage_grp", Site = "site", MSI = "msi")
fr <- list()
for (nm in names(sub)) { v <- sub[[nm]]$var; lv <- sub[[nm]]$lv; adj <- adjv(excl_map[[nm]])
  f_int <- fm(c(paste0("glp2r_z*", v), adj)); f_noint <- fm(c("glp2r_z", v, adj)); pint <- anova(coxph(f_noint, cc5), coxph(f_int, cc5))[2, "Pr(>|Chi|)"]
  for (l in lv) { z <- cc5[cc5[[v]] == l, ]; if (sum(z$event) >= 5) { t <- cox_tab(coxph(fm(c("glp2r_z", adj)), z))[1, ]; fr[[length(fr) + 1]] <- data.frame(Subgroup = nm, Level = l, n = nrow(z), events = sum(z$event), HR = t$HR, lo = t$lo, hi = t$hi, p = t$p, p_interaction = pint) } } }
fr <- do.call(rbind, fr); pu <- unique(fr[, c("Subgroup", "p_interaction")]); pu$fdr <- p.adjust(pu$p_interaction, "BH"); fr$p_int_FDR <- pu$fdr[match(fr$Subgroup, pu$Subgroup)]
print(fr, digits = 3); write_csv_utf8(fr, file.path(RES, "cox_subgroups.csv")); saveRDS(fr, file.path(RES, "cox_subgroups.rds"))
## ---- 6. Multiple imputation sensitivity on all OS patients
imp_d <- dd[, c("time", "event", "glp2r_z", "age10", "sex", "stage_grp", "site", "msi")]; imp_d$cumhaz <- nelsonaalen(imp_d, time, event); imp_d$stage_grp <- factor(imp_d$stage_grp, ordered = FALSE)
cat("\nmissingness in OS set:\n"); print(colSums(is.na(imp_d)))
M <- 30; imp <- mice(imp_d, m = M, maxit = 10, seed = SEED, printFlag = FALSE)
fits <- lapply(seq_len(M), function(i) coxph(fm(c("glp2r_z", "age10", "sex", "stage_grp", "site", "msi")), complete(imp, i)))
Q <- sapply(fits, function(f) coef(f)["glp2r_z"]); U <- sapply(fits, function(f) vcov(f)["glp2r_z", "glp2r_z"]); Qb <- mean(Q); Ub <- mean(U); B <- var(Q); Tt <- Ub + (1 + 1 / M) * B
nu <- (M - 1) * (1 + Ub / ((1 + 1 / M) * B))^2; ci <- Qb + c(-1, 1) * qt(.975, nu) * sqrt(Tt); pv <- 2 * pt(-abs(Qb / sqrt(Tt)), nu)
cat(sprintf("MI (m=%d, n=%d): GLP2R HR per SD = %s, p=%s\n", M, nrow(imp_d), fmt_hr(exp(Qb), exp(ci[1]), exp(ci[2])), fmt_p(pv)))
add("mi_glp2r_hr", sprintf("%s; p=%s; n=%d", fmt_hr(exp(Qb), exp(ci[1]), exp(ci[2])), fmt_p(pv), nrow(imp_d))); add("mi_m", M)
add("n_missing_stage_os", sum(is.na(dd$stage_grp))); add("n_missing_msi_os", sum(is.na(dd$msi)))
## ---- 7. Discrimination: time-dependent AUC of GLP2R, and incremental C-index
set.seed(SEED); tr <- timeROC(T = dd$time, delta = dd$event, marker = -dd$glp2r, cause = 1, times = c(12, 36, 60), iid = TRUE)
ta <- data.frame(time_m = c(12, 36, 60), AUC = tr$AUC, se = tr$inference$vect_sd_1, lo = tr$AUC - 1.96 * tr$inference$vect_sd_1, hi = tr$AUC + 1.96 * tr$inference$vect_sd_1)
print(ta); write_csv_utf8(ta, file.path(RES, "glp2r_timeROC.csv")); saveRDS(tr, file.path(RES, "glp2r_timeROC.rds"))
for (i in 1:3) add(paste0("glp2r_auc_", c(12, 36, 60)[i], "m"), sprintf("%.3f (%.3f-%.3f)", ta$AUC[i], ta$lo[i], ta$hi[i]), "GLP2R alone (low = high risk), IPCW")
cb <- concordance(coxph(fm(c(f_base, "msi")), cc5))$concordance; cf <- concordance(m2)$concordance
bs <- replicate(500, { s <- cc5[sample.int(nrow(cc5), replace = TRUE), ]; tryCatch(concordance(coxph(fm(c("glp2r_z", f_base, "msi")), s))$concordance - concordance(coxph(fm(c(f_base, "msi")), s))$concordance, error = function(e) NA) })
cat(sprintf("C-index base %.3f; +GLP2R %.3f; delta %.4f (bootstrap 95%% CI %.4f to %.4f)\n", cb, cf, cf - cb, quantile(bs, .025, na.rm = TRUE), quantile(bs, .975, na.rm = TRUE)))
add("cindex_base", round(cb, 3)); add("cindex_with_glp2r", round(cf, 3)); add("cindex_delta_ci", sprintf("%.4f (%.4f to %.4f)", cf - cb, quantile(bs, .025, na.rm = TRUE), quantile(bs, .975, na.rm = TRUE)))
lr_add <- anova(mods[["M2 without GLP2R (same n)"]], m2); add("lrt_glp2r_added_p", signif(lr_add[2, "Pr(>|Chi|)"], 4))
## ---- 8. other endpoint: progression-free interval (CDR)
dp <- d[!is.na(d$cdr_PFI_time) & d$cdr_PFI_time > 0 & !is.na(d$cdr_PFI), ]; dp$glp2r_z <- as.numeric(scale(dp$glp2r)); dp$age10 <- dp$age_y / 10; dp$time <- dp$cdr_PFI_time / 30.4375; dp$event <- dp$cdr_PFI
cp <- dp[complete.cases(dp[, c(f_base, "msi")]), ]; tp <- cox_tab(coxph(fm(c("glp2r_z", f_base, "msi")), cp))[1, ]; tu <- cox_tab(coxph(fm("glp2r_z"), dp))[1, ]
add("pfi_n", nrow(dp)); add("pfi_events", sum(dp$event)); add("pfi_unadj_glp2r", sprintf("%s; p=%s", fmt_hr(tu$HR, tu$lo, tu$hi), fmt_p(tu$p))); add("pfi_adj_glp2r", sprintf("%s; p=%s; n=%d", fmt_hr(tp$HR, tp$lo, tp$hi), fmt_p(tp$p), nrow(cp)))
cat("PFI unadj:", fmt_hr(tu$HR, tu$lo, tu$hi), fmt_p(tu$p), " adj:", fmt_hr(tp$HR, tp$lo, tp$hi), fmt_p(tp$p), "\n")
## ---- 9. GLP1R KM (negative control) and reconciliation with the original 199-patient Cox
dd$g1 <- factor(ifelse(dd$glp1r > median(dd$glp1r), "High", "Low"), levels = c("Low", "High")); l1 <- survdiff(Surv(time, event) ~ g1, dd); add("glp1r_km_logrank_p", signif(1 - pchisq(l1$chisq, 1), 4))
f_cox <- file.path(Sys.getenv("GLP2R_ORIG", "V:"), "cox_full.rds")   # optional: object from an earlier analysis, not distributed
if (file.exists(f_cox)) {
orig <- readRDS(f_cox)$data; cat("\noriginal cox_full$data columns:", paste(colnames(orig), collapse = ", "), " n =", nrow(orig), "\n")
add("orig_cox_n", nrow(orig))
}
saveRDS(list(dd = dd, cc5 = cc5, m2 = m2, mods = mods), file.path(RES, "survival_objects.rds"))
rm_write(do.call(rbind, RM), "survival"); log_msg("done"); sink()
