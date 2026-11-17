## 00_common.R : shared paths, packages, theme and helpers (sourced by every analysis script)
options(width = 200, warn = 1, stringsAsFactors = FALSE, scipen = 999)
Sys.setenv(LANGUAGE = "en")
RB   <- Sys.getenv("GLP2R_RB", "W:/Re_build")    # root holding 05_Analysis/data_cache and results (read)
OUT  <- Sys.getenv("GLP2R_OUT", RB)              # root for figures and tables (write)
LIB  <- file.path(RB, "_work", "Rlib")
## the bundled library was built for R 4.4; under newer R use the user library instead (see S00_install_r46.R)
if (getRversion() < "4.5") { .libPaths(c(LIB, .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
DC   <- file.path(RB, "05_Analysis", "data_cache")
RES  <- file.path(RB, "05_Analysis", "results")
LOGD <- file.path(RB, "05_Analysis", "logs")
FIGM <- file.path(OUT, "02_Figures", "Main")
FIGS <- file.path(OUT, "02_Figures", "Supplementary")
TABD <- file.path(OUT, "03_Tables")
RMD  <- RES
if (!identical(OUT, RB)) { LOGD <- file.path(OUT, "logs"); RMD <- file.path(OUT, "results_presubmission") }   # keep RB read-only when OUT differs
for (d in c(LOGD, RMD, FIGM, FIGS, TABD)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
SEED <- 20240617L

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(survival)
})

## ---- publication theme / palettes (colour-blind safe; min text 6 pt)
col_tn   <- c(Normal = "#4C78A8", Tumor = "#D95F02")
col_grp  <- c(Low = "#3B6FB6", High = "#C8472B")
col_msi  <- c("MSS/MSI-L" = "#7F7F7F", "MSI-H" = "#B2182B")
col_site <- c(COAD = "#1B9E77", READ = "#7570B3")
theme_pub <- function(base = 7) {
  theme_classic(base_size = base) %+replace%
    theme(text = element_text(colour = "black"),
          axis.text = element_text(colour = "black", size = base - 0.5),
          axis.title = element_text(size = base),
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_text(size = base + 1, face = "bold", hjust = 0, margin = margin(b = 2)),
          legend.title = element_text(size = base - 0.5), legend.text = element_text(size = base - 1),
          legend.key.size = unit(3, "mm"), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = base, face = "bold"),
          plot.tag = element_text(size = base + 3, face = "bold"),
          plot.margin = margin(3, 4, 3, 3))
}
theme_set(theme_pub())

## save a figure as 300-dpi PNG + vector PDF; widths in mm (single 85, double 180)
save_fig <- function(p, name, w = 180, h = 120, dir = FIGM) {
  f_png <- file.path(dir, paste0(name, ".png")); f_pdf <- file.path(dir, paste0(name, ".pdf"))
  ragg::agg_png(f_png, width = w, height = h, units = "mm", res = 300, background = "white")
  print(p); invisible(dev.off())
  grDevices::cairo_pdf(f_pdf, width = w / 25.4, height = h / 25.4)
  print(p); invisible(dev.off())
  invisible(c(f_png, f_pdf))
}

## ---- formatting helpers
fmt_p <- function(p, digits = 2) {
  vapply(p, function(z) if (is.na(z)) "NA" else if (z < 1e-300) "<1e-300" else if (z < 0.001) formatC(z, format = "e", digits = digits) else formatC(z, format = "f", digits = 3), "")
}
fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
cox_tab <- function(fit) {
  s <- summary(fit)
  data.frame(term = rownames(s$coefficients), HR = s$conf.int[, "exp(coef)"], lo = s$conf.int[, "lower .95"], hi = s$conf.int[, "upper .95"],
             p = s$coefficients[, ncol(s$coefficients)], row.names = NULL)
}
write_csv_utf8 <- function(df, path) { utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8", na = "") }
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
## collapse an Ensembl-keyed matrix to HGNC symbols (keep the gene with the highest mean per duplicated symbol)
to_symbol_matrix <- function(m, gmap) {
  s <- gmap$symbol[match(rownames(m), gmap$ensembl_full)]; ok <- !is.na(s); m <- m[ok, , drop = FALSE]; s <- s[ok]
  o <- order(rowMeans(m), decreasing = TRUE); m <- m[o, , drop = FALSE]; s <- s[o]; d <- !duplicated(s); m <- m[d, , drop = FALSE]; rownames(m) <- s[d]; m
}
## Mann-Whitney AUC (probability that x > y, ties = 0.5)
auc_mw <- function(x, y) { x <- x[!is.na(x)]; y <- y[!is.na(y)]; r <- rank(c(x, y)); (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y)) }
boot_ci <- function(fun, n, B = 2000, seed = SEED) { set.seed(seed); v <- replicate(B, fun(sample.int(n, n, replace = TRUE))); stats::quantile(v, c(0.025, 0.975), na.rm = TRUE) }
## accumulate key numbers for the results ledger
rm_write <- function(df, name) write_csv_utf8(df, file.path(RMD, paste0("rm_", name, ".csv")))
