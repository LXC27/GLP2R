## 110_build_response_docs.R : build the cover letter, the English response to reviewers and the Chinese revision log (raw docx)
## from the text sources in 04_Response/src. Word finishing (page set-up, footer, PDF) is done by 100b_word_finish.ps1.
if (getRversion() < "4.5") { .libPaths(c(file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "_work", "Rlib"), .libPaths())) } else { ul <- Sys.getenv("R_LIBS_USER"); if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths())) }
suppressPackageStartupMessages({ library(officer); library(flextable) })

RB <- Sys.getenv("GLP2R_RB", "W:/Re_build"); SRCD <- file.path(RB, "04_Response/src"); OUTD <- file.path(RB, "04_Response")
FONT <- "Times New Roman"; EA <- "Microsoft YaHei"; PAGE_W <- 6.69; CNQ <- FALSE
rd <- function(f) { con <- file(f, encoding = "UTF-8"); on.exit(close(con)); readLines(con, warn = FALSE) }

## ---------- inline markup (same grammar as 100_build_manuscript.R; no reference numbering) ----------
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
mk_runs <- function(txt, size, bold = FALSE, italic = FALSE, highlight = FALSE, color = "black") {
  rr <- split_italics(parse_inline(txt))
  lapply(rr, function(r) ftext(r$t, fp_text(color = color, font.size = size, bold = bold || r$b, italic = italic || r$i, font.family = FONT, eastasia.family = EA,
    vertical.align = if (r$sup) "superscript" else if (r$sub) "subscript" else "baseline", shading.color = if (highlight) "yellow" else "transparent")))
}
add_par <- function(x, txt, size = 11, bold = FALSE, italic = FALSE, align = "left", before = 0, after = 6, spacing = 1.15,
                    style = NULL, highlight = FALSE, keep_next = FALSE, color = "black", quote = FALSE, bullet = FALSE) {
  fp <- if (quote) fp_par(text.align = "left", line_spacing = 1.05, padding.top = before, padding.bottom = after, padding.left = 10,
                          border.left = fp_border(color = "grey55", width = 2.5), keep_with_next = keep_next)
        else fp_par(text.align = align, line_spacing = spacing, padding.top = before, padding.bottom = after, keep_with_next = keep_next)
  runs <- mk_runs(txt, size, bold, italic, highlight, color)
  if (bullet) runs <- c(list(ftext("\u2022", fp_text(font.size = size, font.family = FONT)), run_tab()), runs)
  p <- do.call(fpar, c(runs, list(fp_p = fp)))
  if (is.null(style)) body_add_fpar(x, p) else body_add_fpar(x, p, style = style)
}
add_head <- function(x, txt, lvl) add_par(x, txt, size = if (lvl == 1) 14 else 11.5, bold = TRUE, align = "left",
  before = if (lvl == 1) 14 else 9, after = 4, style = paste("heading", lvl), keep_next = TRUE)

## ---------- tables ----------
sup_cells <- function(ft, d, fs) {
  pat <- "10\\^[-\u2212]?[0-9]+"
  for (i in seq_len(nrow(d))) for (j in seq_len(ncol(d))) {
    v <- as.character(d[i, j]); if (is.na(v) || !grepl(pat, v)) next
    m <- gregexpr(pat, v)[[1]]; ml <- attr(m, "match.length"); pos <- 1; ch <- list()
    for (k in seq_along(m)) {
      if (m[k] > pos) ch[[length(ch) + 1]] <- as_chunk(substr(v, pos, m[k] - 1))
      tok <- substr(v, m[k], m[k] + ml[k] - 1); ch[[length(ch) + 1]] <- as_chunk("10"); ch[[length(ch) + 1]] <- as_sup(sub("^10\\^", "", tok))
      pos <- m[k] + ml[k]
    }
    if (pos <= nchar(v)) ch[[length(ch) + 1]] <- as_chunk(substr(v, pos, nchar(v)))
    ft <- compose(ft, i = i, j = j, value = do.call(as_paragraph, ch))
  }
  ft
}
add_table <- function(x, b) {
  csvf <- normalizePath(file.path(SRCD, gsub("\\\\", "/", b$csv)), winslash = "/", mustWork = TRUE)
  d <- read.csv(csvf, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character")
  if (CNQ) d[] <- lapply(d, function(v) chartr("\u201c\u201d", "\u300c\u300d", v))
  if (!is.null(b$hdr)) names(d) <- b$hdr
  x <- add_par(x, sub("^((Table|\u8868) ?[A-Z0-9]+[.\u3000 ])", "@b{\\1}", b$text), size = 9, align = "left", before = 6, after = 3, spacing = 1.0, keep_next = TRUE)
  ft <- flextable(d); ft <- theme_booktabs(ft)
  ft <- font(ft, fontname = FONT, eastasia.family = EA, part = "all"); ft <- fontsize(ft, size = b$font, part = "all")
  ft <- bold(ft, part = "header"); ft <- align(ft, align = "left", part = "all")
  ft <- padding(ft, padding.top = 1.5, padding.bottom = 1.5, padding.left = 2.5, padding.right = 2.5, part = "all")
  if (b$bold1) ft <- bold(ft, j = 1, part = "body")
  ft <- sup_cells(ft, d, b$font)
  w <- if (!is.null(b$widths)) b$widths else if (ncol(d) == 3) c(0.32, 0.46, 0.22) else rep(1, ncol(d))
  w <- w / sum(w) * PAGE_W
  ft <- set_table_properties(ft, layout = "fixed"); for (j in seq_along(w)) ft <- width(ft, j = j, width = w[j])
  x <- body_add_flextable(x, ft, align = "center")
  add_par(x, "", size = 4, after = 4)
}

## ---------- parser ----------
parse_blocks <- function(L) {
  out <- list()
  for (ln in L) {
    if (grepl("^!TITLE ", ln)) out[[length(out) + 1]] <- list(type = "title", text = sub("^!TITLE ", "", ln))
    else if (grepl("^!SUB ", ln)) out[[length(out) + 1]] <- list(type = "sub", text = sub("^!SUB ", "", ln))
    else if (grepl("^!TAB\\|", ln)) {
      p <- strsplit(ln, "|", fixed = TRUE)[[1]][-1]
      b <- list(type = "tab", csv = p[1], text = p[2], font = 9, hdr = NULL, widths = NULL, bold1 = FALSE)
      for (o in p[-(1:2)]) {
        if (grepl("^font=", o)) b$font <- as.numeric(sub("font=", "", o))
        else if (grepl("^hdr=", o)) b$hdr <- strsplit(sub("^hdr=", "", o), ";")[[1]]
        else if (grepl("^widths=", o)) b$widths <- as.numeric(strsplit(sub("^widths=", "", o), ";")[[1]])
        else if (o == "bold1") b$bold1 <- TRUE
      }
      out[[length(out) + 1]] <- b
    }
    else if (grepl("^## ", ln)) out[[length(out) + 1]] <- list(type = "h2", text = sub("^## ", "", ln))
    else if (grepl("^# ", ln)) out[[length(out) + 1]] <- list(type = "h1", text = sub("^# ", "", ln))
    else if (grepl("^> ", ln)) out[[length(out) + 1]] <- list(type = "quote", text = sub("^> ", "", ln))
    else if (grepl("^- ", ln)) out[[length(out) + 1]] <- list(type = "bullet", text = sub("^- ", "", ln))
    else if (grepl("^!", ln)) stop("unknown directive: ", ln)
    else out[[length(out) + 1]] <- list(type = "para", text = ln)
  }
  out
}

build_doc <- function(src, out, body_size = 11, spacing = 1.15, letter = FALSE, cn = FALSE) {
  CNQ <<- cn
  L <- rd(file.path(SRCD, src)); L <- L[nzchar(trimws(L))]
  if (cn) L <- chartr("\u201c\u201d", "\u300c\u300d", L)
  x <- read_docx(); nb <- 0
  for (b in parse_blocks(L)) {
    nb <- nb + 1
    x <- switch(b$type,
      title = add_par(x, b$text, size = 17, bold = TRUE, align = "center", before = 6, after = 8, spacing = 1.2),
      sub = add_par(x, b$text, size = 9.5, italic = TRUE, align = "center", after = 12, spacing = 1.05),
      h1 = add_head(x, b$text, 1),
      h2 = add_head(x, b$text, 2),
      quote = add_par(x, b$text, size = body_size - 1, italic = TRUE, color = "#404040", quote = TRUE, before = 2, after = 5),
      bullet = add_par(x, b$text, size = body_size, spacing = spacing, after = 3, bullet = TRUE),
      tab = add_table(x, b),
      para = {
        hl <- letter && grepl("^\\[", b$text)
        txt <- if (letter && grepl("^Re: ", b$text)) paste0("@b{", b$text, "}") else b$text
        add_par(x, txt, size = if (grepl("^@b\\{Changes", b$text)) body_size - 1 else body_size, spacing = spacing, highlight = hl,
                align = if (letter) "left" else "justify", after = if (grepl("^@b\\{Changes", b$text)) 8 else 6)
      })
  }
  print(x, target = file.path(OUTD, out)); cat("built", out, "(", nb, "blocks )\n")
}

build_doc("cover_letter.txt", "Cover_letter_revised_raw.docx", body_size = 11, spacing = 1.2, letter = TRUE)
build_doc("response_EN.txt", "Response_to_Reviewers_EN_raw.docx", body_size = 10.5, spacing = 1.12)
build_doc("revision_log_CN.txt", "revision_log_CN_raw.docx", body_size = 10.5, spacing = 1.2, cn = TRUE)
