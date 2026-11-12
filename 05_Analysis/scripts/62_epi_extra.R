## 62_epi_extra.R : DNMT/TET regulators, TMB, somatic-mutation context (cBioPortal TCGA PanCancer Atlas MC3 calls) and the data frames needed for Figure 6
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
sink(file.path(LOGD, "62_epi_extra.txt"), split = TRUE)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; vsa <- readRDS(file.path(DC, "vst_all_624.rds"))
de <- readRDS(file.path(DC, "de_tumor_normal.rds")); dp <- readRDS(file.path(DC, "de_paired.rds"))
h <- readRDS(file.path(DC, "hm450_subset.rds")); B <- h$beta; pr <- h$probes; gi <- readRDS(file.path(DC, "gistic_thresholded_coadread.rds"))
RM <- list(); add <- function(k, v, desc = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = desc)
d$pid <- substr(rownames(d), 1, 12)
pc <- function(x, y, Z) { ok <- complete.cases(x, y, Z); Zm <- model.matrix(~ ., as.data.frame(Z[ok, , drop = FALSE]))[, -1, drop = FALSE]; rx <- resid(lm(rank(x[ok]) ~ Zm)); ry <- resid(lm(rank(y[ok]) ~ Zm)); r <- cor(rx, ry); n <- sum(ok); tt <- r * sqrt((n - 2 - ncol(Zm)) / (1 - r^2)); c(rho = r, p = 2 * pt(-abs(tt), n - 2 - ncol(Zm)), n = n) }
adj <- c("purity_absolute", "msi", "site")

## ---- 1. DNMT/TET regulators: tumour vs normal and correlation with GLP2R in tumours
reg <- c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3")
rr <- do.call(rbind, lapply(reg, function(g) {
  i <- which(de$symbol == g)[1]; j <- which(dp$symbol == g)[1]
  x <- vsa[g, rownames(d)]; a <- cor.test(d$glp2r, x, method = "spearman", exact = FALSE); q <- pc(d$glp2r, x, d[, adj])
  data.frame(gene = g, log2FC_unpaired = de$log2FoldChange[i], padj_unpaired = de$padj[i], log2FC_paired = dp$log2FoldChange[j], padj_paired = dp$padj[j], rho = unname(a$estimate), p = a$p.value, rho_adj = unname(q["rho"]), p_adj = unname(q["p"]), n_adj = unname(q["n"]))
}))
rr$FDR_rho <- p.adjust(rr$p, "BH"); rr$FDR_rho_adj <- p.adjust(rr$p_adj, "BH"); print(rr, digits = 3); write_csv_utf8(rr, file.path(RES, "dnmt_tet_regulators.csv"))
for (k in seq_len(nrow(rr))) add(paste0("reg_", rr$gene[k]), sprintf("log2FC=%.2f (padj=%s); paired log2FC=%.2f (padj=%s); rho=%.3f (p=%s); adj rho=%.3f (p=%s)", rr$log2FC_unpaired[k], fmt_p(rr$padj_unpaired[k]), rr$log2FC_paired[k], fmt_p(rr$padj_paired[k]), rr$rho[k], fmt_p(rr$p[k]), rr$rho_adj[k], fmt_p(rr$p_adj[k])))

## ---- 2. TMB
tm <- d[!is.na(d$tmb), c("glp2r", "tmb", "msi", "purity_absolute", "site")]; tm$ltmb <- log10(tm$tmb + 0.1)
a <- cor.test(tm$glp2r, tm$ltmb, method = "spearman", exact = FALSE); qa <- pc(tm$glp2r, tm$ltmb, tm[, c("msi", "site")])
ms <- tm[which(tm$msi == "MSS/MSI-L"), ]; b <- cor.test(ms$glp2r, ms$ltmb, method = "spearman", exact = FALSE)
cat(sprintf("TMB all: n=%d rho=%.3f p=%s; MSI/site-adjusted rho=%.3f p=%s; MSS-only n=%d rho=%.3f p=%s\n", nrow(tm), a$estimate, fmt_p(a$p.value), qa["rho"], fmt_p(qa["p"]), nrow(ms), b$estimate, fmt_p(b$p.value)))
add("tmb_adj_msi_site", sprintf("rho=%.3f p=%s n=%d", qa["rho"], fmt_p(qa["p"]), qa["n"]), "TMB vs GLP2R, adjusted for MSI and site")

## ---- 3. somatic mutations (cBioPortal API, TCGA COADREAD PanCancer Atlas, 'sequenced' sample list)
genes <- c(GLP2R = 9248, APC = 324, TP53 = 7157, KRAS = 3845, PIK3CA = 5290, TTN = 7273, SYNE1 = 23345, FBXW7 = 55294, BRAF = 673, SMAD4 = 4089, ARID1A = 8289, TCF7L2 = 6934)
body <- file.path(tempdir(), "mut_body.json"); outf <- file.path(DC, "cbio_mutations_selected.json")
writeLines(sprintf('{"entrezGeneIds":[%s],"sampleListId":"coadread_tcga_pan_can_atlas_2018_sequenced"}', paste(genes, collapse = ",")), body)
if (!file.exists(outf) || file.size(outf) < 5000) system2("curl.exe", c("-s", "-X", "POST", "-H", "Content-Type:application/json", "-d", paste0("@", body), "--max-time", "300", "-o", outf, "https://www.cbioportal.org/api/molecular-profiles/coadread_tcga_pan_can_atlas_2018_mutations/mutations/fetch?projection=SUMMARY"))
if (file.exists(outf) && file.size(outf) < 5000) file.remove(outf); mu <- jsonlite::fromJSON(outf); cat("mutation records:", nrow(mu), "\n"); mu$symbol <- names(genes)[match(mu$entrezGeneId, genes)]; mu$pid <- substr(mu$sampleId, 1, 12)
sl <- system2("curl.exe", c("-s", "https://www.cbioportal.org/api/sample-lists/coadread_tcga_pan_can_atlas_2018_sequenced"), stdout = TRUE); sl <- jsonlite::fromJSON(paste(sl, collapse = ""))
seq_pid <- unique(substr(sl$sampleIds, 1, 12)); cat("sequenced patients:", length(seq_pid), "\n")
tum <- d[d$pid %in% seq_pid, ]; cat("tumours in our cohort with sequencing:", nrow(tum), "\n")
nonsil <- !mu$mutationType %in% c("Silent", "Intron", "RNA", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "Splice_Region")
mu <- mu[nonsil, ]
freq <- do.call(rbind, lapply(names(genes), function(g) { m <- unique(mu$pid[mu$symbol == g]); m <- m[m %in% tum$pid]; data.frame(gene = g, n_mut = length(m), n = nrow(tum), pct = 100 * length(m) / nrow(tum)) }))
print(freq); write_csv_utf8(freq, file.path(RES, "mutation_frequency.csv"))
gm <- mu[mu$symbol == "GLP2R" & mu$pid %in% tum$pid, c("sampleId", "proteinChange", "mutationType")]; print(gm); write_csv_utf8(gm, file.path(RES, "glp2r_mutations.csv"))
add("glp2r_mut_n", sprintf("%d/%d (%.2f%%)", freq$n_mut[freq$gene == "GLP2R"], nrow(tum), freq$pct[freq$gene == "GLP2R"]), "non-silent GLP2R mutations, cBioPortal PanCancer Atlas sequenced tumours also in the RNA-seq set")
add("glp2r_mut_types", paste(names(table(gm$mutationType)), table(gm$mutationType), sep = "=", collapse = "; "))
## GLP2R expression by driver-gene mutation status (descriptive; all genes and the MSS-only subset)
dr <- c("APC", "TP53", "KRAS", "PIK3CA", "BRAF", "SMAD4", "FBXW7", "ARID1A", "TCF7L2")
mo <- do.call(rbind, lapply(dr, function(g) { m <- unique(mu$pid[mu$symbol == g]); tum$mut <- tum$pid %in% m
  w <- wilcox.test(glp2r ~ mut, tum); s <- tum[which(tum$msi == "MSS/MSI-L"), ]; w2 <- wilcox.test(glp2r ~ mut, s)
  q <- pc(tum$glp2r, as.numeric(tum$mut), tum[, c("purity_absolute", "msi", "site")])
  data.frame(gene = g, n_mut = sum(tum$mut), n_wt = sum(!tum$mut), med_mut = median(tum$glp2r[tum$mut]), med_wt = median(tum$glp2r[!tum$mut]), p_all = w$p.value, p_mss = w2$p.value, n_mut_mss = sum(s$pid %in% m), rho_adj = unname(q["rho"]), p_adj = unname(q["p"])) }))
mo$FDR_all <- p.adjust(mo$p_all, "BH"); mo$FDR_mss <- p.adjust(mo$p_mss, "BH"); mo$FDR_adj <- p.adjust(mo$p_adj, "BH"); print(mo, digits = 3); write_csv_utf8(mo, file.path(RES, "mutation_vs_glp2r.csv"))
for (k in seq_len(nrow(mo))) add(paste0("mut_", mo$gene[k]), sprintf("mut n=%d wt n=%d; median VST %.2f vs %.2f; p=%s (FDR %s); MSS-only p=%s (FDR %s); adj rho=%.3f p=%s", mo$n_mut[k], mo$n_wt[k], mo$med_mut[k], mo$med_wt[k], fmt_p(mo$p_all[k]), fmt_p(mo$FDR_all[k]), fmt_p(mo$p_mss[k]), fmt_p(mo$FDR_mss[k]), mo$rho_adj[k], fmt_p(mo$p_adj[k])))
add("mut_cohort", sprintf("%d tumours with MC3 mutation calls (PanCancer Atlas sequenced list)", nrow(tum)))

## ---- 4. data frames for Figure 6 (MLH1 positive control; CNV)
ml <- pr[grepl("(^|,)MLH1(,|$)", pr$gene), ]; typ <- substr(colnames(B), 14, 15)
bt <- B[, typ == "01", drop = FALSE]; colnames(bt) <- substr(colnames(bt), 1, 12); bt <- bt[, !duplicated(colnames(bt)), drop = FALSE]
common <- intersect(colnames(bt), d$pid); dtum <- d[match(common, d$pid), ]; dtum$mlh1 <- vsa["MLH1", rownames(dtum)]; dtum$mlh1_meth <- colMeans(bt[ml$id, common, drop = FALSE], na.rm = TRUE)
colnames(gi) <- substr(colnames(gi), 1, 12); cn_ok <- intersect(colnames(gi), d$pid); dc <- d[match(cn_ok, d$pid), ]; dc$cnv <- as.numeric(gi["GLP2R", cn_ok]); dc$cn <- factor(dc$cnv, levels = -2:2, labels = c("Deep loss", "Shallow loss", "Diploid", "Gain", "Amp"))
saveRDS(list(reg = rr, tm = tm, freq = freq, mo = mo, gm = gm, mlh1 = dtum[, c("pid", "mlh1", "mlh1_meth", "msi", "glp2r")], cnv = dc[, c("pid", "cnv", "cn", "glp2r", "purity_absolute", "msi", "site")], n_seq = nrow(tum)), file.path(RES, "epi_extra_objects.rds"))
rm_write(do.call(rbind, RM), "mut"); log_msg("done"); sink()
