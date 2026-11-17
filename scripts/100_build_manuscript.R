## 100_build_manuscript.R : assemble the revised manuscript and supplementary materials (docx) from the
## text sources in 01_Manuscript/src, with automatic Vancouver numbering of the verified references.
if (getRversion() < "4.5") { .libPaths(c(file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "_work", "Rlib"), .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
suppressPackageStartupMessages({ library(officer); library(flextable); library(magick) })
source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "100a_refdb.R"))

RB <- Sys.getenv("GLP2R_RB", "W:/Re_build"); SRC <- file.path(RB, "01_Manuscript/src"); FIGD <- file.path(RB, "02_Figures"); TABD <- file.path(RB, "03_Tables")
FONT <- "Times New Roman"; BODY <- 11; LEG <- 8.5
rd <- function(f) { con <- file(file.path(SRC, f), encoding = "UTF-8"); on.exit(close(con)); readLines(con, warn = FALSE) }
lines <- unlist(lapply(c("ms_01_front_intro.txt", "ms_02_methods.txt", "ms_03_results.txt", "ms_04_discussion.txt", "ms_05_legends.txt"), rd))
lines <- lines[nzchar(trimws(lines))]
sp <- which(lines == "!SUPPSTART")
main_l <- lines[seq_len(sp - 1)]; supp_l <- lines[-seq_len(sp)]

## ---------- block parser ----------
parse_blocks <- function(L) {
  out <- list()
  for (ln in L) {
    if (grepl("^!TITLE ", ln)) out[[length(out) + 1]] <- list(type = "title", text = sub("^!TITLE ", "", ln))
    else if (grepl("^!SHORTTITLE ", ln)) out[[length(out) + 1]] <- list(type = "short", text = sub("^!SHORTTITLE ", "", ln))
    else if (grepl("^!AUTHORS ", ln)) out[[length(out) + 1]] <- list(type = "authors", text = sub("^!AUTHORS ", "", ln))
    else if (grepl("^!KEYWORDS ", ln)) out[[length(out) + 1]] <- list(type = "keywords", text = sub("^!KEYWORDS ", "", ln))
    else if (grepl("^!TAB ", ln)) { p <- strsplit(sub("^!TAB ", "", ln), "|", fixed = TRUE)[[1]]
      out[[length(out) + 1]] <- list(type = "tab", csv = p[1], text = p[2], font = as.numeric(sub("font=", "", p[3]))) }
    else if (grepl("^!FIG\\|", ln)) { p <- strsplit(ln, "|", fixed = TRUE)[[1]]
      out[[length(out) + 1]] <- list(type = "fig", path = gsub("\\\\", "/", p[2]), label = p[3], text = paste(p[-(1:3)], collapse = "|")) }
    else if (grepl("^## ", ln)) out[[length(out) + 1]] <- list(type = "h2", text = sub("^## ", "", ln))
    else if (grepl("^# ", ln)) out[[length(out) + 1]] <- list(type = "h1", text = sub("^# ", "", ln))
    else if (grepl("^!", ln)) stop("unknown directive: ", ln)
    else out[[length(out) + 1]] <- list(type = "para", text = ln)
  }
  out
}
mb <- parse_blocks(main_l); sb <- parse_blocks(supp_l)

## ---------- citations ----------
cite_pat <- "\\[@[A-Za-z0-9;@ ]+\\]"
keys_of <- function(cc) trimws(strsplit(gsub("[\\[\\]@ ]", "", cc, perl = TRUE), ";")[[1]])
cited <- character(0)
register <- function(txt) { m <- regmatches(txt, gregexpr(cite_pat, txt))[[1]]; for (cc in m) for (k in keys_of(cc)) if (!(k %in% cited)) cited <<- c(cited, k) }
for (b in c(mb, sb)) register(b$text)
miss <- setdiff(cited, names(ref_keys)); if (length(miss)) stop("unknown citation keys: ", paste(miss, collapse = ", "))
dropped <- c("colaprico2015", "davis2007"); if (any(dropped %in% cited)) stop("a dropped reference is cited")
collapse_nums <- function(n) { n <- sort(unique(n)); o <- character(); i <- 1
  while (i <= length(n)) { j <- i; while (j < length(n) && n[j + 1] == n[j] + 1) j <- j + 1
    o <- c(o, if (j - i >= 2) paste0(n[i], "\u2013", n[j]) else as.character(n[i:j])); i <- j + 1 }
  paste(o, collapse = ",") }
fmt_cites <- function(txt) {
  m <- gregexpr(cite_pat, txt); mm <- regmatches(txt, m)
  if (length(mm[[1]])) regmatches(txt, m) <- list(vapply(mm[[1]], function(cc) paste0("[", collapse_nums(match(keys_of(cc), cited)), "]"), ""))
  txt }

## ---------- inline markup ----------
find_close <- function(ch, open_idx) { d <- 1; i <- open_idx + 1
  while (i <= length(ch)) { if (ch[i] == "{") d <- d + 1 else if (ch[i] == "}") { d <- d - 1; if (d == 0) return(i) }; i <- i + 1 }
  stop("unbalanced braces") }
parse_inline <- function(s, bold = FALSE, italic = FALSE, sup = FALSE, sub = FALSE) {
  ch <- strsplit(s, "")[[1]]; n <- length(ch); runs <- list(); buf <- character(); i <- 1
  flush <- function() { if (length(buf)) { runs[[length(runs) + 1]] <<- list(t = paste(buf, collapse = ""), b = bold, i = italic, sup = sup, sub = sub); buf <<- character() } }
  while (i <= n) {
    tag <- NULL
    if (ch[i] == "@" && i + 2 <= n && ch[i + 1] %in% c("b", "i") && ch[i + 2] == "{") tag <- list(k = ch[i + 1], o = i + 2)
    else if (ch[i] %in% c("^", "_") && i + 1 <= n && ch[i + 1] == "{") tag <- list(k = ch[i], o = i + 1)
    if (!is.null(tag)) {
      cl <- find_close(ch, tag$o); inner <- paste(ch[(tag$o + 1):(cl - 1)], collapse = ""); flush()
      sub_runs <- switch(tag$k,
        b = parse_inline(inner, TRUE, italic, sup, sub), i = parse_inline(inner, bold, TRUE, sup, sub),
        "^" = parse_inline(inner, bold, italic, TRUE, FALSE), "_" = parse_inline(inner, bold, italic, FALSE, TRUE))
      runs <- c(runs, sub_runs); i <- cl + 1
    } else { buf <- c(buf, ch[i]); i <- i + 1 }
  }
  flush(); runs
}
ital_pat <- "\\bP(?=\\s*[=<>\u2264\u2265]|\\s+values?\\b|\\s+for\\b)|\\bn(?=\\s?=)"
split_italics <- function(runs) {
  out <- list()
  for (r in runs) {
    if (r$i || r$sup || r$sub || !grepl(ital_pat, r$t, perl = TRUE)) { out[[length(out) + 1]] <- r; next }
    m <- gregexpr(ital_pat, r$t, perl = TRUE)[[1]]; pos <- 1
    for (k in seq_along(m)) {
      if (m[k] > pos) out[[length(out) + 1]] <- modifyList(r, list(t = substr(r$t, pos, m[k] - 1)))
      out[[length(out) + 1]] <- modifyList(r, list(t = substr(r$t, m[k], m[k] + attr(m, "match.length")[k] - 1), i = TRUE))
      pos <- m[k] + attr(m, "match.length")[k]
    }
    if (pos <= nchar(r$t)) out[[length(out) + 1]] <- modifyList(r, list(t = substr(r$t, pos, nchar(r$t))))
  }
  out
}
mk_runs <- function(txt, size, bold = FALSE, highlight = FALSE) {
  rr <- split_italics(parse_inline(fmt_cites(txt)))
  lapply(rr, function(r) ftext(r$t, fp_text(color = "black", font.size = size, bold = bold || r$b, italic = r$i, font.family = FONT,
    vertical.align = if (r$sup) "superscript" else if (r$sub) "subscript" else "baseline", shading.color = if (highlight) "yellow" else "transparent")))
}
plain_text <- function(txt) paste(vapply(parse_inline(fmt_cites(txt)), function(r) if (r$sup) paste0("^", r$t) else r$t, ""), collapse = "")
count_w <- function(t) length(strsplit(trimws(plain_text(t)), "\\s+")[[1]])

## ---------- document helpers ----------
add_par <- function(x, txt, size = BODY, bold = FALSE, align = "justify", before = 0, after = 6, spacing = 1.15, style = NULL, highlight = FALSE, keep_next = FALSE) {
  fp <- fp_par(text.align = align, line_spacing = spacing, padding.top = before, padding.bottom = after, keep_with_next = keep_next)
  p <- do.call(fpar, c(mk_runs(txt, size, bold, highlight), list(fp_p = fp)))
  if (is.null(style)) body_add_fpar(x, p) else body_add_fpar(x, p, style = style)
}
add_head <- function(x, txt, lvl) add_par(x, txt, size = if (lvl == 1) 14 else 12, bold = TRUE, align = "left", before = if (lvl == 1) 14 else 8, after = 6, style = paste("heading", lvl), keep_next = TRUE)
add_table <- function(x, b) {
  d <- read.csv(file.path(TABD, b$csv), check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
  x <- add_par(x, sub("^(Table [0-9S]+\\.)", "@b{\\1}", b$text), size = 9, align = "left", before = 6, after = 3, spacing = 1.0, keep_next = TRUE)
  ft <- flextable(d); ft <- theme_booktabs(ft); ft <- font(ft, fontname = FONT, part = "all"); ft <- fontsize(ft, size = b$font, part = "all")
  ft <- bold(ft, part = "header"); ft <- align(ft, j = 2:ncol(d), align = "center", part = "all"); ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- padding(ft, padding.top = 1.2, padding.bottom = 1.2, padding.left = 2.5, padding.right = 2.5, part = "all")
  if (ncol(d) == 8) { w <- c(1.55, 0.95, 0.95, 0.95, 0.95, 0.45, 0.65, 0.6); ft <- bold(ft, j = 1, part = "body") } else w <- c(1.9, 1.0, 0.55, 1.0, 0.55, 1.0, 0.55)
  w <- w / sum(w) * 6.69
  ft <- set_table_properties(ft, layout = "fixed"); for (j in seq_along(w)) ft <- width(ft, j = j, width = w[j])
  x <- body_add_flextable(x, ft, align = "center")
  add_par(x, "", size = 4, after = 4)
}
add_fig <- function(x, b, max_h = 7.1) {
  f <- file.path(FIGD, b$path); i <- image_info(image_read(f)); asp <- i$height / i$width
  w <- min(6.69, max_h / asp)
  x <- body_add_fpar(x, fpar(external_img(f, width = w, height = w * asp), fp_p = fp_par(text.align = "center", padding.bottom = 6, keep_with_next = TRUE)))
  add_par(x, b$text, size = LEG, align = "justify", spacing = 1.0, after = 0)
}
new_doc <- function() read_docx()

## ---------- plain text log (for number QA) ----------
plain <- character(0); addp <- function(t) plain <<- c(plain, plain_text(t))

## ---------- word counts (pre-pass) ----------
words_abs <- 0; words_text <- 0; sec <- ""; n_fig_main <- sum(vapply(mb, function(b) b$type == "fig", NA)); n_tab_main <- sum(vapply(mb, function(b) b$type == "tab", NA))
for (b in mb) { if (b$type == "h1") sec <- b$text
  if (b$type == "para") { if (sec == "Abstract") words_abs <- words_abs + count_w(b$text) else if (grepl("^[1-5]\\. ", sec)) words_text <- words_text + count_w(b$text) } }
rdb <- build_refdb()
rt <- paste0(seq_along(cited), ". ", rdb[cited])

## ---------- MAIN DOCUMENT ----------
x <- new_doc(); nm <- 0; refs_done <- FALSE
add_refs <- function(x) { x <- body_add_break(x); x <- add_head(x, "References", 1); for (r in rt) x <- add_par(x, r, size = 10, align = "left", after = 3, spacing = 1.05); x }
for (b in mb) {
  tp <- b$type
  if (tp == "title") {
    x <- add_par(x, b$text, size = 17, bold = TRUE, align = "center", before = 60, after = 14, spacing = 1.2); addp(b$text)
  } else if (tp == "short") {
    x <- add_par(x, paste0("@i{Running title:} ", b$text), size = 10, align = "center", after = 18, spacing = 1.0); addp(b$text)
  } else if (tp == "authors") {
    x <- add_par(x, b$text, size = BODY, align = "center", after = 18, highlight = TRUE); addp(b$text)
    x <- add_par(x, sprintf("Abstract: %d words. Main text (Introduction to Conclusions): %d words. Main figures: %d. Tables: %d. References: %d. Supplementary figures: 6; supplementary tables: 8.", words_abs, words_text, n_fig_main, n_tab_main, length(cited)), size = 10, align = "center", after = 6)
    x <- body_add_break(x)
  } else if (tp == "keywords") {
    x <- add_par(x, paste0("@b{Keywords:} ", b$text), size = 10.5, align = "left", before = 6, after = 6); addp(paste("Keywords:", b$text))
    x <- body_add_break(x)
  } else if (tp == "h1") {
    if (b$text == "Figure legends") { x <- add_refs(x); refs_done <- TRUE; x <- body_add_break(x); x <- add_head(x, "Figures and figure legends", 1); next }
    x <- add_head(x, b$text, 1); addp(paste0("\n", b$text))
  } else if (tp == "h2") {
    x <- add_head(x, b$text, 2); addp(b$text)
  } else if (tp == "para") {
    x <- add_par(x, b$text, highlight = grepl("^\\[Author", b$text), size = BODY); addp(b$text)
  } else if (tp == "tab") {
    x <- add_table(x, b); addp(b$text)
  } else if (tp == "fig") {
    if (nm > 0) x <- body_add_break(x)
    x <- add_fig(x, b); addp(b$text); nm <- nm + 1
  }
}
if (!refs_done) stop("no figure section found")
print(x, target = file.path(RB, "01_Manuscript/manuscript_revised_raw.docx"))

## ---------- reference outputs ----------
con <- file(file.path(RB, "06_References/references_final.txt"), "w", encoding = "UTF-8"); writeLines(rt, con); close(con)
ck <- data.frame(number = seq_along(cited), key = cited, orig_no = unname(ref_keys[cited]), stringsAsFactors = FALSE)
write.csv(ck, file.path(RB, "06_References/reference_key_map.csv"), row.names = FALSE, fileEncoding = "UTF-8")
unused <- setdiff(names(ref_keys), cited)
cat("citations:", length(cited), "| uncited keys in db:", paste(unused, collapse = ", "), "\n")
cat("abstract words:", words_abs, "| main text words (sections 1-5):", words_text, "| figures:", nm, "| tables:", n_tab_main, "\n")
con <- file(file.path(RB, "05_Analysis/results/manuscript_plaintext.txt"), "w", encoding = "UTF-8"); writeLines(plain, con); close(con)
writeLines(c(sprintf("abstract_words=%d", words_abs), sprintf("main_text_words=%d", words_text)), file.path(RB, "_work/wordcounts.txt"))

## ---------- SUPPLEMENTARY DOCUMENT ----------
y <- new_doc()
y <- add_par(y, "Supplementary materials", size = 17, bold = TRUE, align = "center", before = 10, after = 8, spacing = 1.2)
y <- add_par(y, "GLP2R expression loss in colorectal cancer: a multi-cohort reanalysis of its prognostic, immune and epigenetic correlates", size = 12, align = "center", after = 14)
first <- TRUE
for (b in sb) {
  if (b$type == "h1") { if (b$text == "Supplementary tables") y <- body_add_break(y); y <- add_head(y, b$text, 1) }
  else if (b$type == "fig") { if (!first) y <- body_add_break(y); first <- FALSE; y <- add_fig(y, b, max_h = 7.3) }
  else if (b$type == "para") y <- add_par(y, b$text, align = "left", spacing = 1.1)
}
print(y, target = file.path(RB, "01_Manuscript/supplementary_materials_raw.docx"))
cat("built raw docx files\n")
