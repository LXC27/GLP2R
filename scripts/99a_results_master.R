source(file.path(Sys.getenv("GLP2R_SCRIPTS", "W:/Re_build/05_Analysis/scripts"), "00_common.R"))
fs <- list.files(RES, pattern = "^rm_.*\\.csv$", full.names = TRUE)
rm_all <- do.call(rbind, lapply(fs, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(module = sub("^rm_", "", tools::file_path_sans_ext(basename(f))), key = d$key, value = d$value,
             description = if ("description" %in% names(d)) d$description else NA_character_, stringsAsFactors = FALSE)
}))
stopifnot(!anyDuplicated(paste(rm_all$module, rm_all$key)))
write_csv_utf8(rm_all, file.path(RB, "05_Analysis", "results_master.csv"))
cat("rows:", nrow(rm_all), "\n")
for (i in seq_len(nrow(rm_all))) cat(sprintf("%s|%s|%s|%s\n", rm_all$module[i], rm_all$key[i], rm_all$value[i], substr(rm_all$description[i], 1, 70)))
