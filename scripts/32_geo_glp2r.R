## 32_geo_glp2r.R : external (microarray) assessment of GLP2R in GSE39582 and GSE17536 + random-effects synthesis with TCGA
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
sink(file.path(LOGD, "32_geo_glp2r.txt"), split = TRUE)
g <- readRDS(file.path(DC, "geo_cohorts.rds")); so <- readRDS(file.path(RES, "survival_objects.rds")); dd <- so$dd
RM <- list(); add <- function(k, v, d = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = d)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
prep39 <- function() { p <- g$GSE39582$pheno; x <- as.numeric(g$GSE39582$glp2r_probe_expr[1, ]); names(x) <- colnames(g$GSE39582$glp2r_probe_expr)
  p$glp2r <- x[rownames(p)]; if (all(is.na(p$glp2r))) p$glp2r <- x[p$geo_accession]
  p$tumour <- p$dataset != "Non Tumoral"; p }
p39 <- prep39(); cat("GSE39582: n =", nrow(p39), " tumour =", sum(p39$tumour), " non-tumoral mucosa =", sum(!p39$tumour), "\n")
## ---- probe signal relative to the whole array (is GLP2R reliably measured on HG-U133 Plus 2?)
for (nm in c("GSE39582", "GSE17536")) { e <- g[[nm]]$expr; pm <- rowMeans(e); gp <- g[[nm]]$glp2r_probe_expr; pr <- mean(pm <= mean(gp)); cat(sprintf("%s: GLP2R probe 221312_at mean log2 = %.2f; percentile among %d probesets = %.1f%%; array 10th/25th/50th percentile of probe means = %.2f/%.2f/%.2f\n", nm, mean(gp), nrow(e), 100 * pr, quantile(pm, .1), quantile(pm, .25), quantile(pm, .5)))
  add(paste0(nm, "_glp2r_probe_percentile"), round(100 * pr, 1), "percentile of GLP2R probe mean among all probesets"); add(paste0(nm, "_glp2r_probe_mean_log2"), round(mean(gp), 2)) }
nt <- p39[!p39$tumour, "glp2r"]; tt <- p39[p39$tumour, "glp2r"]; wt <- wilcox.test(tt, nt); cat(sprintf("GSE39582 GLP2R probe: tumour median %.3f (n=%d) vs non-tumoral %.3f (n=%d); Wilcoxon p=%s; AUC(non-tumour > tumour)=%.3f\n", median(tt), length(tt), median(nt), length(nt), fmt_p(wt$p.value), auc_mw(nt, tt)))
add("gse39582_tumour_vs_normal_p", signif(wt$p.value, 4)); add("gse39582_n_tumour", length(tt)); add("gse39582_n_normal", length(nt)); add("gse39582_glp2r_auc_normal_gt_tumour", round(auc_mw(nt, tt), 3))
## ---- survival data
t39 <- p39[p39$tumour, ]; t39$time <- num(t39$`os.delay (months)`); t39$event <- num(t39$`os.event`); t39$age10 <- num(t39$`age.at.diagnosis (year)`) / 10; t39$sex2 <- factor(t39$Sex)
t39$stage <- factor(num(t39$`tnm.stage`), levels = 1:4, labels = c("I", "II", "III", "IV")); t39$mmr <- factor(t39$`mmr.status`, levels = c("pMMR", "dMMR")); t39$loc <- factor(ifelse(t39$tumor.location %in% c("distal", "proximal"), t39$tumor.location, NA))
t39 <- t39[!is.na(t39$time) & t39$time > 0 & !is.na(t39$event), ]; cat("GSE39582 tumours with OS:", nrow(t39), " events:", sum(t39$event), "\n"); print(table(t39$stage, useNA = "ifany"))
p17 <- g$GSE17536$pheno; x17 <- as.numeric(g$GSE17536$glp2r_probe_expr[1, ]); names(x17) <- colnames(g$GSE17536$glp2r_probe_expr); p17$glp2r <- x17[rownames(p17)]; if (all(is.na(p17$glp2r))) p17$glp2r <- x17[p17$geo_accession]
t17 <- p17; t17$time <- num(t17$`overall survival follow-up time`); t17$event <- as.numeric(t17$`overall_event (death from any cause)` == "death"); t17$age10 <- num(t17$age) / 10; t17$sex2 <- factor(t17$gender)
t17$stage <- factor(num(t17$ajcc_stage), levels = 1:4, labels = c("I", "II", "III", "IV")); t17 <- t17[!is.na(t17$time) & t17$time > 0, ]; cat("GSE17536 with OS:", nrow(t17), " events:", sum(t17$event), "\n"); print(table(t17$stage, useNA = "ifany"))
run <- function(df, lab, extra = NULL) { df$z <- as.numeric(scale(df$glp2r)); df$grp <- factor(ifelse(df$glp2r > median(df$glp2r), "High", "Low"), levels = c("Low", "High"))
  f0 <- coxph(Surv(time, event) ~ z, df); vars <- c("z", "age10", "sex2", "stage", extra); cc <- df[complete.cases(df[, setdiff(vars, "z")]), ]; f1 <- coxph(as.formula(paste("Surv(time, event) ~", paste(vars, collapse = "+"))), cc)
  fg <- coxph(Surv(time, event) ~ grp, df); lr <- survdiff(Surv(time, event) ~ grp, df); t0 <- cox_tab(f0)[1, ]; t1 <- cox_tab(f1)[1, ]; tg <- cox_tab(fg)[1, ]
  cat(sprintf("\n%s: unadj per-SD HR %s p=%s (n=%d,ev=%d) | adj %s p=%s (n=%d,ev=%d) | median-split HR(high vs low) %s, log-rank p=%s\n", lab, fmt_hr(t0$HR, t0$lo, t0$hi), fmt_p(t0$p), f0$n, f0$nevent, fmt_hr(t1$HR, t1$lo, t1$hi), fmt_p(t1$p), f1$n, f1$nevent, fmt_hr(tg$HR, tg$lo, tg$hi), fmt_p(1 - pchisq(lr$chisq, 1))))
  list(unadj = t0, adj = t1, med = tg, p_lr = 1 - pchisq(lr$chisq, 1), n0 = f0$n, e0 = f0$nevent, n1 = f1$n, e1 = f1$nevent, df = df, fit_adj = f1) }
r39 <- run(t39, "GSE39582 (tumours)"); r39m <- run(t39, "GSE39582 (+ MMR)", extra = "mmr"); r17 <- run(t17, "GSE17536")
## same-covariate (age, sex, stage) TCGA estimate for the synthesis
dd$sex2 <- dd$sex; ct <- dd[complete.cases(dd[, c("age10", "sex", "stage_grp")]), ]; ft <- coxph(Surv(time, event) ~ glp2r_z + age10 + sex + stage_grp, ct); tt3 <- cox_tab(ft)[1, ]; cat(sprintf("\nTCGA (age, sex, stage; n=%d): %s p=%s\n", ft$n, fmt_hr(tt3$HR, tt3$lo, tt3$hi), fmt_p(tt3$p)))
est <- data.frame(cohort = c("TCGA-COADREAD (RNA-seq)", "GSE39582 (HG-U133 Plus 2)", "GSE17536 (HG-U133 Plus 2)"), n = c(ft$n, r39$n1, r17$n1), events = c(ft$nevent, r39$e1, r17$e1),
  logHR = c(log(tt3$HR), log(r39$adj$HR), log(r17$adj$HR)), se = c((log(tt3$hi) - log(tt3$lo)) / 3.92, (log(r39$adj$hi) - log(r39$adj$lo)) / 3.92, (log(r17$adj$hi) - log(r17$adj$lo)) / 3.92))
w <- 1 / est$se^2; fe <- sum(w * est$logHR) / sum(w); Q <- sum(w * (est$logHR - fe)^2); df_ <- nrow(est) - 1; C <- sum(w) - sum(w^2) / sum(w); tau2 <- max(0, (Q - df_) / C); I2 <- max(0, (Q - df_) / Q) * 100
wr <- 1 / (est$se^2 + tau2); re <- sum(wr * est$logHR) / sum(wr); se_re <- sqrt(1 / sum(wr)); prc <- 2 * pnorm(-abs(re / se_re))
cat(sprintf("Random-effects (DL) per-SD HR = %s, p=%s; fixed-effect HR=%.3f; I2=%.0f%%; tau2=%.4f; Q p=%s\n", fmt_hr(exp(re), exp(re - 1.96 * se_re), exp(re + 1.96 * se_re)), fmt_p(prc), exp(fe), I2, tau2, fmt_p(1 - pchisq(Q, df_))))
est$HR <- exp(est$logHR); est$lo <- exp(est$logHR - 1.96 * est$se); est$hi <- exp(est$logHR + 1.96 * est$se); est$w_re <- round(100 * wr / sum(wr), 1)
est <- rbind(est, data.frame(cohort = "Random-effects summary", n = sum(est$n), events = sum(est$events), logHR = re, se = se_re, HR = exp(re), lo = exp(re - 1.96 * se_re), hi = exp(re + 1.96 * se_re), w_re = 100))
write_csv_utf8(est, file.path(RES, "glp2r_meta_cohorts.csv")); saveRDS(list(est = est, I2 = I2, tau2 = tau2, p_re = prc, p_Q = 1 - pchisq(Q, df_), r39 = r39[c("unadj", "adj", "med", "p_lr", "n0", "e0", "n1", "e1")], r39m = r39m[c("adj", "n1", "e1")], r17 = r17[c("unadj", "adj", "med", "p_lr", "n0", "e0", "n1", "e1")], t39 = t39, t17 = t17), file.path(RES, "geo_objects.rds"))
for (lab in c("r39", "r39m", "r17")) { r <- get(lab); nm <- c(r39 = "gse39582", r39m = "gse39582_mmr", r17 = "gse17536")[lab]
  add(paste0(nm, "_adj_hr"), sprintf("%s; p=%s; n=%d; events=%d", fmt_hr(r$adj$HR, r$adj$lo, r$adj$hi), fmt_p(r$adj$p), r$n1, r$e1)); if (lab != "r39m") { add(paste0(nm, "_unadj_hr"), sprintf("%s; p=%s; n=%d; events=%d", fmt_hr(r$unadj$HR, r$unadj$lo, r$unadj$hi), fmt_p(r$unadj$p), r$n0, r$e0)); add(paste0(nm, "_median_logrank_p"), signif(r$p_lr, 4)) } }
add("meta_re_hr", sprintf("%s; p=%s", fmt_hr(exp(re), exp(re - 1.96 * se_re), exp(re + 1.96 * se_re)), fmt_p(prc))); add("meta_I2_pct", round(I2)); add("tcga_age_sex_stage_hr", sprintf("%s; p=%s; n=%d", fmt_hr(tt3$HR, tt3$lo, tt3$hi), fmt_p(tt3$p), ft$n))
rm_write(do.call(rbind, RM), "geo"); log_msg("done"); sink()
