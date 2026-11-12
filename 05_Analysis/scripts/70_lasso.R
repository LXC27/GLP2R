## 70_lasso.R : re-built GLP-pathway LASSO-Cox signature with pre-specified filtering, repeated CV, bootstrap stability, permutation null and external validation
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressMessages({ library(glmnet); library(timeROC) })
sink(file.path(LOGD, "70_lasso.txt"), split = TRUE)
RM <- list(); add <- function(k, v, d = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = d)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
so <- readRDS(file.path(RES, "survival_objects.rds")); dd <- so$dd; coh <- readRDS(file.path(DC, "cohort.rds")); vt <- readRDS(file.path(DC, "vst_tumor_591.rds")); g <- readRDS(file.path(DC, "geo_cohorts.rds"))
t1 <- read.csv(file.path(Sys.getenv("GLP2R_INPUT", "W:/Supplementary"), "Table_S1_GLP_genes.csv"), stringsAsFactors = FALSE); cand <- unique(t1$Gene); cat("candidate genes:", length(cand), "\n")
bar <- if (all(rownames(dd) %in% colnames(vt))) rownames(dd) else dd$barcode; stopifnot(all(bar %in% colnames(vt))); dd$bar <- bar
## ---- pre-specified feature filter: detected (>=10 raw counts) in >=20% of tumours, and measurable on both microarrays
cnt <- coh$counts; gm <- coh$gmap; cn_sym <- gm$symbol[match(rownames(cnt), gm$ensembl_full)]
if (all(is.na(cn_sym))) cn_sym <- gm$symbol[match(sub("\\..*", "", rownames(cnt)), gm$ensembl)]
stopifnot(all(bar %in% colnames(cnt)))
det <- sapply(cand, function(s) { r <- which(cn_sym == s); if (!length(r)) return(NA_real_); mean(cnt[r[which.max(rowMeans(cnt[r, bar, drop = FALSE]))], bar] >= 10) })
in_vt <- cand %in% rownames(vt); in_geo <- cand %in% rownames(g$GSE39582$expr) & cand %in% rownames(g$GSE17536$expr)
feat <- data.frame(gene = cand, frac_detected = round(det, 3), in_tcga_vst = in_vt, in_both_arrays = in_geo); feat$keep <- feat$in_tcga_vst & feat$in_both_arrays & !is.na(feat$frac_detected) & feat$frac_detected >= 0.20
write_csv_utf8(feat, file.path(RES, "lasso_feature_filter.csv")); print(feat[!feat$keep, ])
genes <- feat$gene[feat$keep]; cat("genes retained:", length(genes), "of", length(cand), "\n"); add("lasso_n_candidates", length(cand)); add("lasso_n_retained", length(genes)); add("lasso_excluded", paste(feat$gene[!feat$keep], collapse = ","))
## ---- stratified 70/30 split (seed fixed)
set.seed(SEED); idx_e <- which(dd$event == 1); idx_c <- which(dd$event == 0)
tr <- c(sample(idx_e, round(0.7 * length(idx_e))), sample(idx_c, round(0.7 * length(idx_c)))); trn <- dd[tr, ]; tst <- dd[-tr, ]
cat(sprintf("train n=%d (events %d); test n=%d (events %d)\n", nrow(trn), sum(trn$event), nrow(tst), sum(tst$event))); add("lasso_train_n", nrow(trn)); add("lasso_train_events", sum(trn$event)); add("lasso_test_n", nrow(tst)); add("lasso_test_events", sum(tst$event))
Xr <- t(vt[genes, dd$bar]); rownames(Xr) <- dd$bar; mu <- colMeans(Xr[trn$bar, ]); sdv <- apply(Xr[trn$bar, ], 2, sd)
Xtr <- scale(Xr[trn$bar, ], mu, sdv); Xte <- scale(Xr[tst$bar, ], mu, sdv); ytr <- Surv(trn$time, trn$event); yte <- Surv(tst$time, tst$event)
## ---- repeated 10-fold CV (20 fold assignments) for lambda
lam_min <- numeric(20); lam_1se <- numeric(20); cvs <- list()
for (r in 1:20) { set.seed(SEED + r); cv <- cv.glmnet(Xtr, ytr, family = "cox", alpha = 1, nfolds = 10); lam_min[r] <- cv$lambda.min; lam_1se[r] <- cv$lambda.1se; if (r == 1) cvs <- cv }
lmin <- exp(median(log(lam_min))); l1se <- exp(median(log(lam_1se))); cat(sprintf("median lambda.min = %.4f; median lambda.1se = %.4f\n", lmin, l1se))
fit <- glmnet(Xtr, ytr, family = "cox", alpha = 1); co_min <- as.matrix(coef(fit, s = lmin))[, 1]; co_1se <- as.matrix(coef(fit, s = l1se))[, 1]
sel_min <- co_min[co_min != 0]; sel_1se <- co_1se[co_1se != 0]; cat("lambda.min genes (", length(sel_min), "):\n"); print(round(sort(sel_min), 3)); cat("lambda.1se genes (", length(sel_1se), "):\n"); print(round(sort(sel_1se), 3))
add("lasso_lambda_min", signif(lmin, 4)); add("lasso_nsel_min", length(sel_min)); add("lasso_sel_min", paste(sprintf("%s(%.3f)", names(sel_min), sel_min), collapse = ";")); add("lasso_nsel_1se", length(sel_1se)); add("lasso_sel_1se", paste(sprintf("%s(%.3f)", names(sel_1se), sel_1se), collapse = ";"))
## ---- bootstrap selection frequency at the (fixed) lambda.min
B <- 200; selmat <- matrix(0, B, length(genes), dimnames = list(NULL, genes)); set.seed(SEED)
for (b in 1:B) { i <- sample(nrow(Xtr), replace = TRUE); cvb <- tryCatch(cv.glmnet(Xtr[i, ], Surv(trn$time[i], trn$event[i]), family = "cox", nfolds = 5), error = function(e) NULL); if (is.null(cvb)) next
  cb <- as.matrix(coef(cvb, s = "lambda.min"))[, 1]; selmat[b, ] <- as.numeric(cb != 0) }
sf <- colMeans(selmat); sf_tab <- data.frame(gene = genes, selection_freq = sf, coef_min = co_min[genes], selected_min = co_min[genes] != 0)[order(-sf), ]; print(head(sf_tab, 15)); write_csv_utf8(sf_tab, file.path(RES, "lasso_bootstrap_selection.csv"))
add("lasso_boot_n_genes_ge_80pct", sum(sf >= 0.8)); add("lasso_boot_n_genes_ge_50pct", sum(sf >= 0.5)); add("lasso_boot_top", paste(sprintf("%s(%.2f)", sf_tab$gene[1:5], sf_tab$selection_freq[1:5]), collapse = ";")); add("lasso_boot_mean_size", round(mean(rowSums(selmat)), 1))
## ---- scores
risk <- function(X, co) as.numeric(X[, names(co), drop = FALSE] %*% co)
cidx <- function(time, event, s) { cc <- concordance(Surv(time, event) ~ s, reverse = TRUE); c(C = unname(cc$concordance), se = sqrt(cc$var)) }
ci_str <- function(x) sprintf("%.3f (%.3f-%.3f)", x["C"], x["C"] - 1.96 * x["se"], x["C"] + 1.96 * x["se"])
auc_t <- function(time, event, s, times = c(12, 36, 60)) { r <- tryCatch(timeROC(T = time, delta = event, marker = s, cause = 1, weighting = "marginal", times = times, iid = FALSE), error = function(e) NULL); if (is.null(r)) rep(NA_real_, length(times)) else as.numeric(r$AUC) }
cohorts <- list()
cohorts$TCGA_train <- data.frame(time = trn$time, event = trn$event, s = risk(Xtr, sel_min))
cohorts$TCGA_test <- data.frame(time = tst$time, event = tst$event, s = risk(Xte, sel_min), age10 = tst$age10, sex = tst$sex, stage = tst$stage_grp)
## GEO cohorts: genes z-scored within cohort
p9 <- g$GSE39582$pheno; p9$tumour <- p9$dataset != "Non Tumoral"; p9$time <- num(p9$`os.delay (months)`); p9$event <- num(p9$`os.event`); p9$age10 <- num(p9$`age.at.diagnosis (year)`) / 10; p9$sex <- factor(p9$Sex); p9$stage <- factor(num(p9$`tnm.stage`), levels = 1:4, labels = c("I", "II", "III", "IV"))
col9 <- if (all(rownames(p9) %in% colnames(g$GSE39582$expr))) rownames(p9) else p9$geo_accession; p9 <- p9[p9$tumour & !is.na(p9$time) & p9$time > 0 & !is.na(p9$event), ]; col9 <- if (all(rownames(p9) %in% colnames(g$GSE39582$expr))) rownames(p9) else p9$geo_accession
p7 <- g$GSE17536$pheno; p7$time <- num(p7$`overall survival follow-up time`); p7$event <- as.numeric(p7$`overall_event (death from any cause)` == "death"); p7$age10 <- num(p7$age) / 10; p7$sex <- factor(p7$gender); p7$stage <- factor(num(p7$ajcc_stage), levels = 1:4, labels = c("I", "II", "III", "IV")); p7 <- p7[!is.na(p7$time) & p7$time > 0, ]
col7 <- if (all(rownames(p7) %in% colnames(g$GSE17536$expr))) rownames(p7) else p7$geo_accession
zg <- function(e, cols) { m <- t(e[genes, cols, drop = FALSE]); scale(m) }
cohorts$GSE39582 <- data.frame(time = p9$time, event = p9$event, s = risk(zg(g$GSE39582$expr, col9), sel_min), age10 = p9$age10, sex = p9$sex, stage = p9$stage)
cohorts$GSE17536 <- data.frame(time = p7$time, event = p7$event, s = risk(zg(g$GSE17536$expr, col7), sel_min), age10 = p7$age10, sex = p7$sex, stage = p7$stage)
## ---- performance table
perf <- do.call(rbind, lapply(names(cohorts), function(n) { x <- cohorts[[n]]; ci <- cidx(x$time, x$event, x$s); a <- auc_t(x$time, x$event, x$s); data.frame(cohort = n, n = nrow(x), events = sum(x$event), C = ci["C"], C_lo = ci["C"] - 1.96 * ci["se"], C_hi = ci["C"] + 1.96 * ci["se"], AUC12 = a[1], AUC36 = a[2], AUC60 = a[3], row.names = NULL) }))
cat("\nPerformance of the lambda.min signature (risk score; higher = higher hazard):\n"); print(perf, digits = 3); write_csv_utf8(perf, file.path(RES, "lasso_performance.csv"))
for (i in seq_len(nrow(perf))) add(paste0("lasso_C_", perf$cohort[i]), sprintf("%.3f (%.3f-%.3f); n=%d events=%d; AUC12/36/60=%s", perf$C[i], perf$C_lo[i], perf$C_hi[i], perf$n[i], perf$events[i], paste(round(c(perf$AUC12[i], perf$AUC36[i], perf$AUC60[i]), 3), collapse = "/")))
## ---- binary stratification (train median; GEO: within-cohort median) + adjusted HR per SD
cut_tr <- median(cohorts$TCGA_train$s); km <- list(); hr_adj <- list()
for (n in names(cohorts)) { x <- cohorts[[n]]; x$grp <- factor(ifelse(x$s > (if (n %in% c("TCGA_train", "TCGA_test")) cut_tr else median(x$s)), "High risk", "Low risk"), levels = c("Low risk", "High risk")); cohorts[[n]] <- x
  lr <- survdiff(Surv(time, event) ~ grp, x); p <- 1 - pchisq(lr$chisq, 1); h <- cox_tab(coxph(Surv(time, event) ~ grp, x))[1, ]; cat(sprintf("%s: median split HR(high vs low risk) %s; log-rank p=%s\n", n, fmt_hr(h$HR, h$lo, h$hi), fmt_p(p))); add(paste0("lasso_km_", n), sprintf("HR %s; logrank p=%s", fmt_hr(h$HR, h$lo, h$hi), fmt_p(p)))
  if ("age10" %in% names(x)) { x$sz <- as.numeric(scale(x$s)); m <- coxph(Surv(time, event) ~ sz + age10 + sex + stage, x); t <- cox_tab(m); t <- t[t$term == "sz", ]; hr_adj[[n]] <- t; cat(sprintf("   adjusted (age, sex, stage; n=%d): signature per SD HR %s, p=%s\n", m$n, fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p))); add(paste0("lasso_adjHR_", n), sprintf("%s; p=%s; n=%d", fmt_hr(t$HR, t$lo, t$hi), fmt_p(t$p), m$n))
    m0 <- coxph(Surv(time, event) ~ age10 + sex + stage, x[complete.cases(x[, c("age10", "sex", "stage")]), ]); m1 <- coxph(Surv(time, event) ~ sz + age10 + sex + stage, x[complete.cases(x[, c("age10", "sex", "stage")]), ]); add(paste0("lasso_Cindex_clin_", n), sprintf("%.3f -> %.3f", concordance(m0)$concordance, concordance(m1)$concordance)); cat(sprintf("   C-index clinical %.3f -> +signature %.3f\n", concordance(m0)$concordance, concordance(m1)$concordance)) } }
## ---- bootstrap CI of test C-index and comparison with GLP2R alone
set.seed(SEED); xt <- cohorts$TCGA_test; bc <- replicate(1000, { i <- sample(nrow(xt), replace = TRUE); tryCatch(cidx(xt$time[i], xt$event[i], xt$s[i])["C"], error = function(e) NA) }); bci <- quantile(bc, c(.025, .975), na.rm = TRUE); cat(sprintf("test C-index bootstrap 95%% CI %.3f-%.3f\n", bci[1], bci[2])); add("lasso_test_C_boot_ci", sprintf("%.3f-%.3f", bci[1], bci[2]))
g2_te <- cidx(tst$time, tst$event, -tst$glp2r); cat(sprintf("GLP2R alone (low expression = high risk) in the same test set: C=%s\n", ci_str(g2_te))); add("glp2r_alone_test_C", ci_str(g2_te))
## ---- permutation null for the whole pipeline (outcome permuted in training set; test C evaluated)
nperm <- 100; pc <- numeric(nperm); set.seed(SEED)
for (k in 1:nperm) { j <- sample(nrow(Xtr)); cvp <- tryCatch(cv.glmnet(Xtr, Surv(trn$time[j], trn$event[j]), family = "cox", nfolds = 10), error = function(e) NULL); if (is.null(cvp)) { pc[k] <- NA; next }
  cp <- as.matrix(coef(cvp, s = "lambda.min"))[, 1]; cp <- cp[cp != 0]; pc[k] <- if (length(cp) == 0) 0.5 else cidx(tst$time, tst$event, risk(Xte, cp))["C"] }
obs <- perf$C[perf$cohort == "TCGA_test"]; pperm <- (1 + sum(pc >= obs, na.rm = TRUE)) / (1 + sum(!is.na(pc))); cat(sprintf("permutation null test C-index: median %.3f (95%% range %.3f-%.3f); observed %.3f; one-sided permutation p=%.3f\n", median(pc, na.rm = TRUE), quantile(pc, .025, na.rm = TRUE), quantile(pc, .975, na.rm = TRUE), obs, pperm))
add("lasso_perm_null_median", round(median(pc, na.rm = TRUE), 3)); add("lasso_perm_p", round(pperm, 3)); add("lasso_perm_prop_empty", mean(pc == 0.5, na.rm = TRUE))
## ---- re-evaluation of the ORIGINAL 16-gene signature (as supplied) to document what was reported
## optional: needs lasso_combined.rds from an earlier analysis (not distributed); skipped when absent
f_lc <- file.path(Sys.getenv("GLP2R_ORIG", "V:"), "lasso_combined.rds")
if (file.exists(f_lc)) {
lc <- readRDS(f_lc); st <- lc$surv_test; cat("\noriginal surv_test columns:", paste(names(st), collapse = ","), "\n"); print(summary(st[[intersect(c("time", "OS.time", "os_time"), names(st))[1]]]))
tcol <- intersect(c("time", "OS.time", "os_time"), names(st))[1]; ecol <- intersect(c("status", "event", "OS", "os_event"), names(st))[1]
if (!is.na(tcol) && !is.na(ecol)) { ro <- as.numeric(lc$X_test[, lc$sel_genes] %*% lc$sel_coefs); tm <- st[[tcol]]; ev <- st[[ecol]]; mscale <- if (max(tm, na.rm = TRUE) > 400) 30.44 else 1
  co <- cidx(tm / mscale, ev, ro); au <- auc_t(tm / mscale, ev, ro); rtr <- as.numeric(lc$X_train[, lc$sel_genes] %*% lc$sel_coefs); trs <- lc$surv_train; ctr <- cidx(trs[[tcol]] / mscale, trs[[ecol]], rtr); atr <- auc_t(trs[[tcol]] / mscale, trs[[ecol]], rtr)
  cat(sprintf("ORIGINAL signature (supplied): train C=%s AUC12/36/60=%s; test (n=%d) C=%s AUC12/36/60=%s\n", ci_str(ctr), paste(round(atr, 3), collapse = "/"), length(tm), ci_str(co), paste(round(au, 3), collapse = "/")))
  add("orig_signature_train", sprintf("C=%s; AUC12/36/60=%s", ci_str(ctr), paste(round(atr, 3), collapse = "/"))); add("orig_signature_test", sprintf("C=%s; AUC12/36/60=%s; n=%d", ci_str(co), paste(round(au, 3), collapse = "/"), length(tm))); add("orig_signature_time_scale_divisor", mscale) }
}
saveRDS(list(genes = genes, feat = feat, trn = trn, tst = tst, sel_min = sel_min, sel_1se = sel_1se, lmin = lmin, l1se = l1se, cvs = cvs, fit = fit, sf_tab = sf_tab, perf = perf, cohorts = cohorts, perm = pc, mu = mu, sdv = sdv, hr_adj = hr_adj, test_boot = bc), file.path(RES, "lasso_objects.rds")); rm_write(do.call(rbind, RM), "lasso"); log_msg("done"); sink()
