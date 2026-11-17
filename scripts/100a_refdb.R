# Reference database: citation key -> Vancouver-style text, built from verified Crossref metadata.
# Sourced by 100_build_manuscript.R; can also be run directly to write 06_References/references_all_verified.txt

ref_keys <- c(
  sung2021 = "1", han2024 = "2", drucker2018 = "3", wang2024 = "5", holst2007 = "6", mayo2003 = "7",
  galon2006 = "11", le2015 = "12", andre2020 = "13", ganesh2019 = "14", baylin2016 = "15", wu2017 = "16",
  flavahan2017 = "17", bennett1998 = "19", popat2005 = "20", jeppesen2011 = "21", langfelder2008 = "23",
  hanzelmann2013 = "24", cristescu2018 = "26", love2014 = "27", simon2011 = "29", colaprico2015 = "30",
  davis2007 = "31", ayers2017 = "32", guinney2015 = "33",
  bushi2025 = "N1", zhong2025 = "N2", zhu2024 = "N3", koehler2015 = "N4", koehler2008 = "N5",
  iakoubov2009 = "N6", bengi2011 = "N7", korner2012 = "N8", yusta2000 = "N9", guan2006 = "N10",
  bjerknes2001 = "N11", drokhlyansky2020 = "N12", lee2020 = "N13", weisenberger2006 = "N14",
  marisa2013 = "N15", tcga2012 = "N16", zeng2018 = "N17", altman1994 = "N18", sterne2009 = "N19",
  vanbuuren2011 = "N20", yoshihara2013 = "N21", aran2015 = "N22", thorsson2018 = "N23", bindea2013 = "N24",
  charoentong2017 = "N25", collins2015 = "N26", mcshane2005 = "N27", gelman2006 = "N28",
  schemper1996 = "N29", benjamini1995 = "N30", grambsch1994 = "N31", blanche2013 = "N32",
  subramanian2005 = "N33", goldman2020 = "N34", mermel2011 = "N36", cortes2017 = "N37", becht2016 = "N39",
  liu2018 = "N40", cerami2012 = "N41", gao2013 = "N42", ellrott2018 = "N43", taylor2018 = "N44",
  kautto2016 = "N45", liberzon2015 = "N46", smith2010 = "N47", zhang2005 = "N48",
  meinshausen2010 = "N49", peduzzi1995 = "N50", grossman2016 = "N51", tibshirani1997 = "N52",
  korotkevich2016 = "N53", gattex2025 = "LABEL"
)

.fix_name <- function(a) {
  a <- trimws(gsub("[()]", "", a))
  p <- strsplit(a, " ")[[1]]
  if (length(p) >= 2) {
    s <- p[-length(p)]
    s <- vapply(s, function(w) if (nchar(w) > 1 && w == toupper(w)) paste0(substr(w, 1, 1), tolower(substring(w, 2))) else w, "")
    a <- paste(c(s, p[length(p)]), collapse = " ")
  }
  a
}

.fix_title <- function(t) {
  t <- gsub("&amp;", "&", gsub("<[^>]+>", "", t))
  t <- gsub("\\s+", " ", trimws(t))
  if (t == toupper(t)) {
    t <- tolower(t); t <- paste0(toupper(substr(t, 1, 1)), substring(t, 2))
    t <- gsub("\\bcox\\b", "Cox", t)
  }
  sub("\\.$", "", t)
}

build_refdb <- function(meta_file = file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "06_References", "crossref_metadata.csv")) {
  m <- read.csv(meta_file, stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = c("NA", ""))
  m$orig_no <- as.character(m$orig_no)
  out <- character(0)
  for (k in names(ref_keys)) {
    if (k == "gattex2025") {
      out[k] <- "Takeda Pharmaceuticals U.S.A., Inc. GATTEX (teduglutide) for injection, for subcutaneous use: prescribing information [Internet]. Silver Spring (MD): US Food and Drug Administration; revised 2025 Sep [cited 2026]. Available from: https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/203441s024lbl.pdf"
      next
    }
    r <- m[m$orig_no == ref_keys[[k]], ]
    if (nrow(r) != 1) stop("reference not found in metadata: ", k)
    au <- trimws(strsplit(r$authors, ";")[[1]])
    au <- vapply(au, .fix_name, "")
    au <- if (length(au) > 6) paste0(paste(au[1:6], collapse = ", "), ", et al") else paste(au, collapse = ", ")
    vol <- if (!is.na(r$volume)) as.character(r$volume) else ""
    iss <- if (!is.na(r$issue) && grepl("^[0-9]+$", r$issue)) paste0("(", r$issue, ")") else ""
    pg  <- if (!is.na(r$page)) r$page else ""
    pgp <- strsplit(pg, "-")[[1]]
    if (length(pgp) == 2 && pgp[1] == pgp[2]) pg <- pgp[1]
    loc <- paste0(r$year, if (nzchar(vol)) paste0(";", vol, iss) else "", if (nzchar(pg)) paste0(":", pg) else "")
    ti <- .fix_title(r$title)
    if (k == "vanbuuren2011") { au <- "van Buuren S, Groothuis-Oudshoorn K"; ti <- "mice: Multivariate imputation by chained equations in R" }
    if (k == "wang2024") loc <- "2024;10(2):256-258"
    if (k == "kautto2016") loc <- "2017;8(5):7452-7463"
    if (k == "gao2013") loc <- "2013;6(269):pl1"
    if (k == "simon2011") loc <- "2011;39(5):1-13"
    if (k == "vanbuuren2011") loc <- "2011;45(3):1-67"
    if (k == "zhang2005") loc <- "2005;4:Article17"
    if (k == "blanche2013") au <- "Blanche P, Dartigues JF, Jacqmin-Gadda H"
    jn <- gsub("&amp;", "&", r$journal)
    out[k] <- gsub("\u2010", "-", sprintf("%s. %s%s %s. %s. doi:%s", au, ti, if (grepl("[?]$", ti)) "" else ".", jn, loc, r$doi))
  }
  out
}

if (sys.nframe() == 0) {
  rd <- build_refdb()
  writeLines(paste0(seq_along(rd), ". ", rd), file.path(Sys.getenv("GLP2R_RB", "W:/Re_build"), "06_References", "references_all_verified.txt"), useBytes = TRUE)
  cat(length(rd), "references written\n")
}
