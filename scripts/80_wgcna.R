## 80_wgcna.R : signed co-expression network in 591 primary tumours (sex-chromosome genes removed) and GLP2R module annotation
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
suppressMessages({ library(WGCNA); library(msigdbr) })
sink(file.path(LOGD, "80_wgcna.txt"), split = TRUE)
options(stringsAsFactors = FALSE); disableWGCNAThreads()
RM <- list(); add <- function(k, v, d = "") RM[[length(RM) + 1]] <<- data.frame(key = k, value = as.character(v), description = d)
ex <- readRDS(file.path(RES, "expr_tests.rds")); d <- ex$d; vt <- readRDS(file.path(DC, "vst_tumor_591.rds")); coh <- readRDS(file.path(DC, "cohort.rds"))
cnt <- coh$counts; gm <- coh$gmap; tum <- colnames(vt); cat("tumours:", length(tum), "\n")
sym <- gm$symbol[match(rownames(cnt), gm$ensembl_full)]; keepg <- rowMeans(cnt[, tum] >= 10) >= 0.20 & !is.na(sym); detg <- unique(sym[keepg]); cat("genes detected (>=10 counts in >=20% tumours):", length(detg), "\n")
c1 <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C1")); ygenes <- unique(c1$gene_symbol[grepl("^chrY", c1$gs_name)]); cat("chrY genes from MSigDB C1 positional sets:", length(ygenes), "\n")
excl <- unique(c(ygenes, "XIST", "TSIX", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y")); use <- setdiff(intersect(detg, rownames(vt)), excl); cat("genes after removal of sex-chromosome genes:", length(use), " (excluded Y/XIST genes present among detected:", length(intersect(detg, excl)), ")\n")
add("wgcna_n_detected", length(detg)); add("wgcna_n_sexchr_removed", length(intersect(detg, excl)))
mad_v <- apply(vt[use, ], 1, mad); top <- names(sort(mad_v, decreasing = TRUE))[1:8000]; if (!"GLP2R" %in% top) top <- c(top, "GLP2R"); datExpr <- t(vt[top, ]); cat("WGCNA input:", dim(datExpr), " GLP2R included:", "GLP2R" %in% colnames(datExpr), "\n")
gsg <- goodSamplesGenes(datExpr, verbose = 0); cat("goodSamplesGenes OK:", gsg$allOK, "\n"); if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
sdt <- hclust(dist(datExpr), "average"); ht <- cutreeStatic(sdt, cutHeight = max(sdt$height) * 0.98, minSize = 10)
cat("sample dendrogram: outlier clusters (kept; no samples removed):", table(ht), "\n")
pw <- 1:30; sft <- pickSoftThreshold(datExpr, powerVector = pw, networkType = "signed", verbose = 0); fi <- sft$fitIndices; print(fi[, c("Power", "SFT.R.sq", "slope", "mean.k.")], digits = 3)
rsq <- -sign(fi$slope) * fi$SFT.R.sq; pick <- fi$Power[which(rsq >= 0.85)[1]]; if (is.na(pick)) { pick <- 12; cat("scale-free R2>=0.85 not reached; WGCNA default for signed networks (n>=40) = 12\n") }; cat("soft-threshold power:", pick, "\n"); add("wgcna_power", pick); add("wgcna_scale_free_R2_at_power", round(rsq[fi$Power == pick], 3)); add("wgcna_input_genes", ncol(datExpr)); add("wgcna_input_samples", nrow(datExpr))
net <- blockwiseModules(datExpr, power = pick, networkType = "signed", TOMType = "signed", minModuleSize = 30, mergeCutHeight = 0.25, maxBlockSize = ncol(datExpr) + 1, numericLabels = FALSE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 1, randomSeed = SEED)
mc <- net$colors; cat("modules:\n"); print(sort(table(mc), decreasing = TRUE)); nmod <- length(setdiff(unique(mc), "grey")); add("wgcna_n_modules", nmod); add("wgcna_n_grey", sum(mc == "grey")); add("wgcna_module_sizes", paste(names(sort(table(mc), decreasing = TRUE)), sort(table(mc), decreasing = TRUE), sep = "=", collapse = "; "))
ME <- net$MEs; colnames(ME) <- sub("^ME", "", colnames(ME)); ME <- ME[, colnames(ME) != "grey", drop = FALSE]
## ---- GLP2R module and module-GLP2R / trait associations
gl_mod <- mc["GLP2R"]; cat("GLP2R module:", gl_mod, "(size", sum(mc == gl_mod), ")\n"); add("wgcna_glp2r_module", gl_mod); add("wgcna_glp2r_module_size", sum(mc == gl_mod))
gx <- d[match(rownames(datExpr), rownames(d)), ]; if (all(is.na(gx$glp2r))) gx <- d[match(rownames(datExpr), d$barcode), ]
if (!all(rownames(datExpr) %in% rownames(d))) stop("sample mapping")
pc <- function(x, y, Z) { ok <- complete.cases(x, y, Z); Zm <- model.matrix(~ ., as.data.frame(Z[ok, , drop = FALSE]))[, -1, drop = FALSE]; rx <- resid(lm(rank(x[ok]) ~ Zm)); ry <- resid(lm(rank(y[ok]) ~ Zm)); r <- cor(rx, ry); n <- sum(ok); tt <- r * sqrt((n - 2 - ncol(Zm)) / (1 - r^2)); c(rho = r, p = 2 * pt(-abs(tt), n - 2 - ncol(Zm)), n = n) }
modtab <- do.call(rbind, lapply(colnames(ME), function(m) { s <- cor.test(ME[, m], gx$glp2r, method = "spearman", exact = FALSE); a <- pc(ME[, m], gx$glp2r, gx[, c("purity_absolute", "msi", "site")])
  data.frame(module = m, size = sum(mc == m), rho_GLP2R = unname(s$estimate), p = s$p.value, rho_partial = unname(a["rho"]), p_partial = unname(a["p"]), rho_MSI_H = unname(cor(ME[, m], as.numeric(gx$msi == "MSI-H"), use = "complete.obs", method = "spearman")), rho_purity = unname(cor(ME[, m], gx$purity_absolute, use = "complete.obs", method = "spearman")), stringsAsFactors = FALSE) }))
modtab$FDR <- p.adjust(modtab$p, "BH"); modtab <- modtab[order(-abs(modtab$rho_GLP2R)), ]; print(modtab, digits = 3); write_csv_utf8(modtab, file.path(RES, "wgcna_module_table.csv"))
cat(sprintf("GLP2R (gene) is assigned to module '%s'; module eigengene vs GLP2R rho=%.3f\n", gl_mod, modtab$rho_GLP2R[modtab$module == gl_mod])); add("wgcna_glp2r_module_ME_rho", round(modtab$rho_GLP2R[modtab$module == gl_mod], 3)); add("wgcna_top_modules_by_rho", paste(sprintf("%s(%.2f)", head(modtab$module, 5), head(modtab$rho_GLP2R, 5)), collapse = "; "))
## ---- hub genes (kME) in the GLP2R module and the strongest-associated module
kme <- cor(datExpr, ME[, colnames(ME)], use = "pairwise.complete.obs"); topmod <- modtab$module[1]
hub <- function(m, n = 15) { gs <- names(mc)[mc == m]; k <- kme[gs, m]; head(sort(k, decreasing = TRUE), n) }
hub_g <- hub(gl_mod); cat("hub genes of GLP2R module:", paste(names(hub_g), collapse = ", "), "\n"); add("wgcna_glp2r_module_hubs", paste(names(hub_g), collapse = ", ")); add("wgcna_glp2r_kME", round(kme["GLP2R", gl_mod], 3)); add("wgcna_glp2r_kME_rank", sum(kme[names(mc)[mc == gl_mod], gl_mod] > kme["GLP2R", gl_mod]) + 1)
## ---- module annotation: hypergeometric enrichment against Hallmark, GO BP and C8 cell-type sets (background = WGCNA input genes)
bg <- colnames(datExpr); H <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H")); GO <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")); C8 <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C8"))
gl <- function(df) lapply(split(df$gene_symbol, df$gs_name), function(x) intersect(unique(x), bg))
sets <- list(Hallmark = gl(H), GO_BP = gl(GO), C8 = gl(C8)); sets <- lapply(sets, function(s) s[lengths(s) >= 10 & lengths(s) <= 500])
enrich <- function(genes, S, collname) { N <- length(bg); k <- length(genes); res <- do.call(rbind, lapply(names(S), function(n) { s <- S[[n]]; x <- length(intersect(genes, s)); data.frame(collection = collname, set = n, overlap = x, set_size = length(s), p = phyper(x - 1, length(s), N - length(s), k, lower.tail = FALSE)) })); res$FDR <- p.adjust(res$p, "BH"); res[order(res$p), ] }
mods_ann <- unique(c(gl_mod, modtab$module[1:3])); en_all <- list()
for (m in mods_ann) { gs <- names(mc)[mc == m]; for (cn in names(sets)) { e <- enrich(gs, sets[[cn]], cn); e$module <- m; en_all[[paste(m, cn)]] <- e; cat(sprintf("\n[%s module, %d genes] top %s sets:\n", m, length(gs), cn)); print(head(e[, c("set", "overlap", "set_size", "p", "FDR")], 5), digits = 3) } }
en <- do.call(rbind, en_all); write_csv_utf8(en, file.path(RES, "wgcna_module_enrichment.csv"))
for (m in mods_ann) { e <- en[en$module == m & en$collection == "Hallmark", ]; add(paste0("wgcna_", m, "_hallmark_top"), paste(head(e$set[e$FDR < 0.05], 4), collapse = "; ")); e <- en[en$module == m & en$collection == "GO_BP", ]; add(paste0("wgcna_", m, "_GO_top"), paste(head(e$set[e$FDR < 0.05], 4), collapse = "; ")); e <- en[en$module == m & en$collection == "C8", ]; add(paste0("wgcna_", m, "_C8_top"), paste(head(e$set[e$FDR < 0.05], 4), collapse = "; ")) }
## ---- module preservation-like check: do the module's top hub genes keep GLP2R correlation in the normal samples? (descriptive)
vsa <- readRDS(file.path(DC, "vst_all_624.rds")); nrm <- rownames(d)[d$condition == "Normal"]; if (length(nrm) < 10) nrm <- colnames(vsa)[substr(colnames(vsa), 14, 15) == "11"]
hg <- intersect(names(hub_g), rownames(vsa)); rn <- sapply(hg, function(g) cor(vsa[g, nrm], vsa["GLP2R", nrm], method = "spearman")); cat("hub-gene Spearman with GLP2R in the normal samples (n =", length(nrm), "): median", round(median(rn), 2), "\n"); add("wgcna_hub_rho_in_normals_median", round(median(rn), 2))
saveRDS(list(mc = mc, ME = ME, modtab = modtab, kme_glp2r = kme["GLP2R", ], hub = hub_g, power = pick, fi = fi, sdt = sdt, en = en, topmod = topmod, gl_mod = gl_mod, sft_rsq = rsq), file.path(RES, "wgcna_objects.rds")); rm_write(do.call(rbind, RM), "wgcna"); log_msg("done"); sink()
