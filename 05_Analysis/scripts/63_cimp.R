## 63_cimp.R -- GLP2R vs CIMP / BRAF / MMR / CIT subtype
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
log_msg("63_cimp start")
set.seed(SEED)

res <- list()
add <- function(cohort, variable, contrast, n1, n0, med1, med0, auc, p, p_strat = NA_real_, stratum = NA_character_, n_strat = NA_integer_, beta_adj = NA_real_, p_adj = NA_real_) {
  res[[length(res) + 1]] <<- data.frame(cohort = cohort, variable = variable, contrast = contrast, n_case = n1, n_ref = n0,
    median_case = med1, median_ref = med0, auc_case_gt_ref = auc, p_wilcoxon = p, stratum = stratum, n_stratum = n_strat,
    p_wilcoxon_stratum = p_strat, beta_adjusted = beta_adj, p_adjusted = p_adj, stringsAsFactors = FALSE)
}
two <- function(y, g, case, ref) {
  ok <- !is.na(y) & g %in% c(case, ref)
  y <- y[ok]; g <- g[ok]
  a <- y[g == case]; b <- y[g == ref]
  list(n1 = length(a), n0 = length(b), m1 = median(a), m0 = median(b), auc = auc_mw(a, b),
       p = suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value))
}

## ---- TCGA (HM450-based CIMP clusters from the TCGA 2012 CRC paper, as distributed by UCSC Xena legacy)
co <- readRDS(file.path(DC, "cohort.rds"))
cl <- co$clin
vst <- readRDS(file.path(DC, "vst_tumor_591.rds"))
hc <- read.csv(file.path(RES, "cohort_clinical_harmonised.csv"), stringsAsFactors = FALSE)
hc <- hc[hc$sample_type == "Primary Tumor", ]
hc$glp2r <- vst["GLP2R", hc$sample]
hc$cimp <- cl[hc$sample, "xena_CIMP"]
hc$cimp[hc$cimp == ""] <- NA
hc$cimp_grp <- ifelse(is.na(hc$cimp), NA, ifelse(hc$cimp == "CIMP.H", "CIMP-high", "other"))
hc$msi[hc$msi == ""] <- NA
n_cimp <- sum(!is.na(hc$cimp))
kw_cimp <- kruskal.test(glp2r ~ factor(cimp), data = hc[!is.na(hc$cimp), ])$p.value
tc <- two(hc$glp2r, hc$cimp_grp, "CIMP-high", "other")
mss <- !is.na(hc$msi) & hc$msi == "MSS/MSI-L"
tcm <- two(hc$glp2r[mss], hc$cimp_grp[mss], "CIMP-high", "other")
d <- hc[!is.na(hc$cimp_grp) & !is.na(hc$msi), ]
fit <- lm(glp2r ~ I(cimp_grp == "CIMP-high") + msi + site, data = d)
sm <- summary(fit)$coefficients
add("TCGA", "CIMP cluster (TCGA 2012)", "CIMP-high vs other clusters", tc$n1, tc$n0, tc$m1, tc$m0, tc$auc, tc$p,
    p_strat = tcm$p, stratum = "MSS/MSI-L only", n_strat = tcm$n1 + tcm$n0, beta_adj = sm[2, 1], p_adj = sm[2, 4])
cat("TCGA CIMP n=", n_cimp, " KW p=", kw_cimp, "\n"); print(tc); print(tcm); print(sm)
d2 <- d[!is.na(d$purity_absolute), ]
fit2 <- lm(glp2r ~ I(cimp_grp == "CIMP-high") + msi + site + purity_absolute, data = d2)
sm2 <- summary(fit2)$coefficients
cimp_pur <- sprintf("beta=%.3f; P=%s; n=%d", sm2[2, 1], fmt_p(sm2[2, 4]), nrow(d2))
cat("purity-adjusted:", cimp_pur, "\n")
cim_tab <- table(hc$cimp)

## ---- GSE39582 (probe 221312_at; near background, interpret with caution)
g <- readRDS(file.path(DC, "geo_cohorts.rds"))$GSE39582
ph <- g$pheno
ph$glp2r <- as.numeric(g$glp2r_probe_expr["221312_at", ph$geo_accession])
ph <- ph[ph$dataset != "Non Tumoral", ]
clean <- function(x) { x <- as.character(x); x[x %in% c("N/A", "", "NA")] <- NA; x }
ph$cimp <- clean(ph$cimp.status); ph$braf <- clean(ph$braf.mutation); ph$mmr <- clean(ph$mmr.status)
ph$kras <- clean(ph$kras.mutation); ph$tp53 <- clean(ph$tp53.mutation); ph$cin <- clean(ph$cin.status); ph$cit <- clean(ph$cit.molecularsubtype)
ph$loc <- clean(ph$tumor.location); ph$sex <- clean(ph$Sex)
ph$age <- suppressWarnings(as.numeric(ph[["age.at.diagnosis (year)"]]))
cat("GSE39582 tumours:", nrow(ph), "\n")
specs <- list(list("CIMP status", "cimp", "+", "-"), list("BRAF mutation", "braf", "M", "WT"), list("MMR status", "mmr", "dMMR", "pMMR"),
              list("KRAS mutation", "kras", "M", "WT"), list("TP53 mutation", "tp53", "M", "WT"), list("CIN status", "cin", "+", "-"))
for (s in specs) {
  v <- s[[2]]; tt <- two(ph$glp2r, ph[[v]], s[[3]], s[[4]])
  ps <- NA; ns <- NA; stratum <- NA; b <- NA; pa <- NA
  if (v != "mmr") {
    k <- !is.na(ph$mmr) & ph$mmr == "pMMR"; t2 <- two(ph$glp2r[k], ph[[v]][k], s[[3]], s[[4]])
    ps <- t2$p; ns <- t2$n1 + t2$n0; stratum <- "pMMR only"
    dd <- ph[!is.na(ph[[v]]) & ph[[v]] %in% c(s[[3]], s[[4]]) & !is.na(ph$mmr) & !is.na(ph$loc) & !is.na(ph$sex) & !is.na(ph$age), ]
    ff <- lm(glp2r ~ I(dd[[v]] == s[[3]]) + mmr + loc + sex + age, data = dd)
    cf <- summary(ff)$coefficients; b <- cf[2, 1]; pa <- cf[2, 4]
  }
  add("GSE39582", s[[1]], paste(s[[3]], "vs", s[[4]]), tt$n1, tt$n0, tt$m1, tt$m0, tt$auc, tt$p, ps, stratum, ns, b, pa)
}
kw_cit <- kruskal.test(glp2r ~ factor(cit), data = ph[!is.na(ph$cit), ])$p.value
kw_cit_pmmr <- kruskal.test(glp2r ~ factor(cit), data = ph[!is.na(ph$cit) & !is.na(ph$mmr) & ph$mmr == "pMMR", ])$p.value
cit_med <- tapply(ph$glp2r[!is.na(ph$cit)], ph$cit[!is.na(ph$cit)], median)
cat("CIT KW p=", kw_cit, " pMMR p=", kw_cit_pmmr, "\n"); print(round(cit_med, 3)); print(table(ph$cit))

out <- do.call(rbind, res)
out$p_fdr <- p.adjust(out$p_wilcoxon, "BH")
num <- sapply(out, is.numeric); out[num] <- lapply(out[num], function(z) signif(z, 4))
write_csv_utf8(out, file.path(RES, "cimp_associations.csv"))
print(out)

fm <- function(r) sprintf("median %.2f vs %.2f; AUC %.2f; P=%s", r$median_case, r$median_ref, r$auc_case_gt_ref, fmt_p(r$p_wilcoxon))
led <- data.frame(key = c("tcga_cimp_n", "tcga_cimp_kw_p", "tcga_cimp_high_vs_other", "tcga_cimp_high_mss_p", "tcga_cimp_adj", "tcga_cimp_adj_purity",
                          "g39582_n_tumours", "g39582_cit_kw_p", "g39582_cit_kw_p_pmmr"),
                  value = c(n_cimp, signif(kw_cimp, 4), fm(out[1, ]), fmt_p(out$p_wilcoxon_stratum[1]),
                            sprintf("beta=%.3f; P=%s", out$beta_adjusted[1], fmt_p(out$p_adjusted[1])), cimp_pur,
                            nrow(ph), signif(kw_cit, 4), signif(kw_cit_pmmr, 4)),
                  description = c("TCGA primary tumours with a CIMP cluster label (UCSC Xena legacy)", "Kruskal-Wallis across the four TCGA methylation clusters",
                                  "GLP2R VST in CIMP-high vs other clusters", "MSS/MSI-L only", "adjusted for MSI and site (linear model on VST)", "adjusted for MSI, site and ABSOLUTE purity",
                                  "GSE39582 tumours (non-tumoral excluded)", "GLP2R probe by CIT subtype C1-C6", "same, pMMR only"))
for (i in 2:nrow(out)) {
  r <- out[i, ]; nm <- tolower(gsub("[^A-Za-z0-9]+", "_", r$variable))
  led <- rbind(led, data.frame(key = paste0("g39582_", nm), value = fm(r), description = paste0("GSE39582 ", r$variable, ": ", r$contrast, " (n=", r$n_case, " vs ", r$n_ref, ")")))
  if (!is.na(r$p_wilcoxon_stratum)) led <- rbind(led, data.frame(key = paste0("g39582_", nm, "_strat_p"), value = fmt_p(r$p_wilcoxon_stratum), description = paste0(r$stratum, "; n=", r$n_stratum)))
  if (!is.na(r$p_adjusted)) led <- rbind(led, data.frame(key = paste0("g39582_", nm, "_adj"), value = sprintf("beta=%.3f; P=%s", r$beta_adjusted, fmt_p(r$p_adjusted)), description = "adjusted for MMR, tumour location, sex, age (linear model)"))
}
rm_write(led, "cimp")
log_msg("63_cimp done")
