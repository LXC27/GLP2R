## S01_presubmission_sensitivity.R : additional sensitivity analyses
## A tissue-source-site restricted tumour/normal comparison; B C-index comparators (M1, M2); C Hartung-Knapp interval; D MSI-stratum heterogeneity; E leave-GLP2R-out eigengene and top-100 hub table; F promoter-restricted MLH1 control.
## Run on 2026-11-11 with R 4.6.1; output: results_presubmission/S01_presubmission_sensitivity.txt.
## Reads the existing cache (GLP2R_RB) and writes only to <GLP2R_OUT>/results_presubmission.
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
OUTD <- file.path(OUT, "results_presubmission"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUTD, "S01_presubmission_sensitivity.txt"), split = TRUE)
blk <- function(name, expr) { cat("\n=====", name, "=====\n"); tryCatch(expr, error = function(e) cat("FAILED:", conditionMessage(e), "\n")) }
pcor_s <- function(x, y, Z) { ok <- complete.cases(x, y, Z); rx <- resid(lm(rank(x[ok]) ~ ., data = Z[ok, , drop = FALSE])); ry <- resid(lm(rank(y[ok]) ~ ., data = Z[ok, , drop = FALSE])); ct <- cor.test(rx, ry); c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok)) }

## ---- A (M2). Tumour vs normal within tissue source sites that contributed normal samples; TSS as covariate
blk("A. tissue-source-site sensitivity", {
  co <- readRDS(file.path(DC, "cohort.rds")); cl <- co$clin; vsa <- readRDS(file.path(DC, "vst_all_624.rds"))
  bc <- colnames(vsa); type <- substr(bc, 14, 15); keep <- type %in% c("01", "11") & bc %in% rownames(cl)[cl$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")]
  df <- data.frame(bc = bc[keep], glp2r = as.numeric(vsa["GLP2R", keep]), cond = ifelse(type[keep] == "11", "Normal", "Tumour"), tss = substr(bc[keep], 6, 7), plate = sapply(strsplit(bc[keep], "-"), `[`, 6))
  df$site <- cl$project_id[match(df$bc, rownames(cl))]; if (all(is.na(df$site))) df$site <- cl$project_id[match(df$bc, cl$barcode)]
  cat("normal samples per TSS:\n"); print(table(df$tss[df$cond == "Normal"]))
  sub <- df[df$tss %in% unique(df$tss[df$cond == "Normal"]), ]
  cat("restricted to TSS with normals: tumours", sum(sub$cond == "Tumour"), " normals", sum(sub$cond == "Normal"), "\n")
  w <- wilcox.test(glp2r ~ cond, sub); cat("Wilcoxon P =", signif(w$p.value, 3), "; AUC(normal > tumour) =", round(auc_mw(sub$glp2r[sub$cond == "Normal"], sub$glp2r[sub$cond == "Tumour"]), 3), "\n")
  f <- lm(glp2r ~ cond + site + tss, sub); print(round(summary(f)$coefficients["condTumour", ], 4))
  per <- do.call(rbind, lapply(split(sub, sub$tss), function(z) if (all(table(factor(z$cond, c("Normal", "Tumour"))) >= 2)) data.frame(tss = z$tss[1], n_normal = sum(z$cond == "Normal"), n_tumour = sum(z$cond == "Tumour"), median_diff = median(z$glp2r[z$cond == "Tumour"]) - median(z$glp2r[z$cond == "Normal"]))))
  print(per); write_csv_utf8(per, file.path(OUTD, "A_tumour_normal_by_TSS.csv"))
  cat("plates of normal samples:\n"); print(table(df$plate[df$cond == "Normal"]))
})

## ---- B (M3). Apparent C-index, base vs +GLP2R, for M1 (n=459) and M2 (n=451)
blk("B. C-index comparators", {
  so <- readRDS(file.path(RES, "survival_objects.rds")); dd <- so$dd
  for (cov in list(M1 = c("age10", "sex", "stage_grp", "site"), M2 = c("age10", "sex", "stage_grp", "site", "msi"))) {
    z <- dd[complete.cases(dd[, cov]), ]; f0 <- coxph(as.formula(paste("Surv(time, event) ~", paste(cov, collapse = "+"))), z); f1 <- update(f0, . ~ . + glp2r_z)
    cat(sprintf("%s: n=%d events=%d; C base %.4f; +GLP2R %.4f; LRT P=%.4f\n", paste(cov, collapse = "+"), nrow(z), sum(z$event), concordance(f0)$concordance, concordance(f1)$concordance, anova(f0, f1)[2, "Pr(>|Chi|)"])) }
})

## ---- C (M3). Hartung-Knapp interval for the three-cohort summary (DL tau2, variance factor truncated at 1)
blk("C. Hartung-Knapp", {
  m <- read.csv(file.path(RES, "glp2r_meta_cohorts.csv")); m <- m[1:3, ]; y <- m$logHR; v <- m$se^2; w <- 1 / v; k <- 3
  Q <- sum(w * (y - sum(w * y) / sum(w))^2); tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (v + tau2); mu <- sum(wr * y) / sum(wr); se <- sqrt(1 / sum(wr))
  qf <- sum(wr * (y - mu)^2) / (k - 1); se_hk <- sqrt(max(1, qf)) * se; ci <- mu + c(-1, 1) * qt(.975, k - 1) * se_hk
  cat(sprintf("DL HR %.3f (%.3f-%.3f); tau2 %.4f; HK %.3f-%.3f; HK P=%.3f\n", exp(mu), exp(mu - 1.96 * se), exp(mu + 1.96 * se), tau2, exp(ci[1]), exp(ci[2]), 2 * pt(-abs(mu / se_hk), k - 1)))
})

## ---- D (M5). MSI-stratum comparison of immune correlations (Fisher z, Fieller SE)
blk("D. MSI-stratified heterogeneity", {
  s <- read.csv(file.path(RES, "immune_msi_stratified.csv")); print(names(s))
  r1 <- s[[grep("MSS.*rho", names(s))[1]]]; r2 <- s[[grep("MSIH.*rho|MSI.H.*rho", names(s))[1]]]; n1 <- s[[grep("MSS.*n$", names(s))[1]]]; n2 <- s[[grep("MSIH.*n$|MSI.H.*n$", names(s))[1]]]
  z <- (atanh(r1) - atanh(r2)) / sqrt(1.06 / (n1 - 3) + 1.06 / (n2 - 3)); p <- 2 * pnorm(-abs(z)); out <- data.frame(cell_type = s[[1]], z = z, p = p, fdr = p.adjust(p, "BH")); print(out); write_csv_utf8(out, file.path(OUTD, "D_MSI_strata_difference.csv"))
})

## ---- E (M5). Turquoise eigengene without GLP2R; full hub table (top 100 by kME)
blk("E. leave-GLP2R-out eigengene and hub table", {
  wo <- readRDS(file.path(RES, "wgcna_objects.rds")); vt <- readRDS(file.path(DC, "vst_tumor_591.rds")); ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d
  mc <- wo$mc; gm <- wo$gl_mod; genes <- intersect(names(mc)[mc == gm], rownames(vt)); cat("module", gm, "genes:", length(genes), "\n")
  X <- t(vt[genes, ]); pc1 <- function(M) { M <- scale(M); u <- svd(M, nu = 1, nv = 0)$u[, 1]; if (cor(u, rowMeans(M)) < 0) u <- -u; u }
  me_all <- pc1(X); me_loo <- pc1(X[, colnames(X) != "GLP2R"]); gx <- d[match(rownames(X), rownames(d)), ]; if (all(is.na(gx$glp2r))) gx <- d[match(rownames(X), d$barcode), ]
  for (nm in c("with GLP2R", "without GLP2R")) { me <- if (nm == "with GLP2R") me_all else me_loo; s <- cor.test(me, gx$glp2r, method = "spearman", exact = FALSE); a <- pcor_s(me, gx$glp2r, gx[, c("purity_absolute", "msi", "site")])
    cat(sprintf("eigengene %s: rho=%.3f (P=%.2g); partial rho=%.3f (P=%.2g, n=%d)\n", nm, s$estimate, s$p.value, a["rho"], a["p"], a["n"])) }
  kme <- cor(X, me_all); hub <- data.frame(gene = rownames(kme), kME = kme[, 1]); hub <- hub[order(-hub$kME), ]; hub$rank <- seq_len(nrow(hub))
  hub$rho_GLP2R_tumour <- sapply(hub$gene, function(g) cor(X[, g], X[, "GLP2R"], method = "spearman"))
  write_csv_utf8(hub[1:100, ], file.path(OUTD, "E_TableS8h_turquoise_top100_hubs.csv")); cat("GLP2R kME rank:", hub$rank[hub$gene == "GLP2R"], "of", nrow(hub), "\n")
})

## ---- F (M6). MLH1 control restricted to promoter-proximal probes (hg19 chr3:37,033,841-37,035,341; verify the TSS coordinate before use)
blk("F. MLH1 promoter-restricted control", {
  ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; vsa <- readRDS(file.path(DC, "vst_all_624.rds")); h <- readRDS(file.path(DC, "hm450_subset.rds")); B <- h$beta; pr <- h$probes
  d$pid <- substr(rownames(d), 1, 12); typ <- substr(colnames(B), 14, 15); bt <- B[, typ == "01", drop = FALSE]; colnames(bt) <- substr(colnames(bt), 1, 12); bt <- bt[, !duplicated(colnames(bt)), drop = FALSE]
  ml <- pr[grepl("(^|,)MLH1(,|$)", pr$gene), ]; cat("MLH1-annotated probes:", nrow(ml), "\n"); print(ml[order(ml$chromStart), c("id", "chrom", "chromStart")])
  prom <- ml[ml$chrom == "chr3" & ml$chromStart >= 37033841 & ml$chromStart <= 37035341, ]; cat("promoter-window probes:", nrow(prom), "\n")
  common <- intersect(colnames(bt), d$pid); dt <- d[match(common, d$pid), ]; dt$mlh1 <- as.numeric(vsa["MLH1", rownames(dt)])
  for (nm in c("all annotated", "promoter window")) { ids <- intersect(if (nm == "all annotated") ml$id else prom$id, rownames(bt)); mm <- colMeans(bt[ids, common, drop = FALSE], na.rm = TRUE); a <- cor.test(mm, dt$mlh1, method = "spearman", exact = FALSE)
    cat(sprintf("MLH1 %s (%d probes): rho=%.3f P=%.2g\n", nm, length(ids), a$estimate, a$p.value)) }
})
sink()
