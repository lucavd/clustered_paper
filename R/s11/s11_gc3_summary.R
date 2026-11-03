# s11_gc3_summary.R -- S11 gate GC3, summary: (a) chain tables, rows of V4_perf_all outside the marginalized model that
# differ between variants nc and cv, by set; (b) every top-level CSV of tables/ archived before the step
# (results_v4/S11/tables_pre_S11.tgz, extracted in /tmp/pre11) against the current file: rows that differ (1e-12) and,
# among them, rows with no mention of the marginalized model in any text column. Output: tables/V8_S11_gc3_summary.csv,
# tables/V8_S11_gc3_perf_all_sets.csv.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setwd("/home/user/MMM")
sm <- function(a, b) {   # TRUE where equal: both NA, identical (also equal infinities), or numeric within 1e-12 relative
  if (is.numeric(a) && is.numeric(b)) { r <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & (a == b | abs(a - b) <= 1e-12 * pmax(1, abs(a)))) }
  else { a <- as.character(a); b <- as.character(b); r <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b) }
  r[is.na(r)] <- FALSE; r }
x <- fread("tables/variant_nc/V4_perf_all.csv.gz"); y <- fread("tables/variant_cv/V4_perf_all.csv.gz")
d <- Reduce(`|`, lapply(names(x), function(cl) !sm(x[[cl]], y[[cl]])))
P <- x[d & !grepl("mmm", method), .N, by = .(set, usab)]; P[, total_rows_set := x[!grepl("mmm", method), .N, by = .(set, usab)][P, on = .(set, usab), N]]
fwrite(P, "tables/V8_S11_gc3_perf_all_sets.csv"); print(P)
pre <- list.files("/tmp/pre11/tables", pattern = "\\.csv$", full.names = TRUE)
OUT <- list()
for (f in pre) {
  g <- file.path("tables", basename(f)); if (!file.exists(g)) { OUT[[f]] <- data.table(table = basename(f), status = "removed"); next }
  a <- tryCatch(fread(f), error = function(e) NULL); b <- tryCatch(fread(g), error = function(e) NULL)
  if (is.null(a) || is.null(b)) { OUT[[f]] <- data.table(table = basename(f), status = "unreadable"); next }
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) { OUT[[f]] <- data.table(table = basename(f), status = "shape differs", rows_pre = nrow(a), rows_now = nrow(b)); next }
  if (ncol(a) == 0 || nrow(a) == 0) { OUT[[f]] <- data.table(table = basename(f), status = "empty"); next }
  dd <- Reduce(`|`, lapply(names(a), function(cl) !sm(a[[cl]], b[[cl]]))); dd[is.na(dd)] <- TRUE   # NA: counted as differing
  ch <- names(a)[sapply(a, is.character)]
  mm <- if (length(ch)) Reduce(`|`, lapply(ch, function(cl) grepl("mmm|MMM|marginali", a[[cl]]))) else rep(FALSE, nrow(a))
  same_bytes <- unname(tools::md5sum(f) == tools::md5sum(g))
  OUT[[f]] <- data.table(table = basename(f), status = if (any(dd)) "differs" else "equal", same_bytes = same_bytes, rows_pre = nrow(a), rows_now = nrow(b),
                         differ = sum(dd), differ_no_mmm = sum(dd & !mm))
}
S <- rbindlist(OUT, fill = TRUE); fwrite(S, "tables/V8_S11_gc3_summary.csv")
print(S[, .N, by = status]); print(S[status != "equal"][order(-differ_no_mmm)], nrows = 200)
