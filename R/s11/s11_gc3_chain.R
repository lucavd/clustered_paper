# s11_gc3_chain.R -- S11 gate GC3, chain part: every table of the variant cv (tables/variant_cv) against the same table of
# the variant nc (tables/variant_nc): rows that differ (tolerance 1e-12 on numbers) must belong to the marginalized model
# (method/strategy column containing "mmm") of scenarios M3-M5 or be summaries over scenarios that include them.
# Output: tables/V8_S11_gc3_chain.csv (one row per table: rows, rows differing, rows differing outside mmm).
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages(library(data.table)); setDTthreads(2L); setwd("/home/user/MMM")
A <- "tables/variant_nc"; B <- "tables/variant_cv"
fs <- intersect(list.files(A, "\\.csv(\\.gz)?$"), list.files(B, "\\.csv(\\.gz)?$"))
fs <- fs[!sapply(fs, function(f) Sys.readlink(file.path(B, f)) != "" && Sys.readlink(file.path(A, f)) != "")]   # both symlinks: shared study table
OUT <- list()
for (f in fs) {
  x <- tryCatch(fread(file.path(A, f)), error = function(e) NULL); y <- tryCatch(fread(file.path(B, f)), error = function(e) NULL)
  if (is.null(x) || is.null(y)) { OUT[[f]] <- data.table(table = f, note = "unreadable"); next }
  if (!identical(names(x), names(y)) || nrow(x) != nrow(y)) { OUT[[f]] <- data.table(table = f, rows = nrow(x), rows_cv = nrow(y), note = "shape differs"); next }
  same <- function(a, b) if (is.numeric(a) && is.numeric(b)) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(a))) else (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
  d <- Reduce(`|`, lapply(names(x), function(cl) !same(x[[cl]], y[[cl]])))
  mc <- intersect(c("method", "strategy", "m", "name", "id", "claim", "variant"), names(x))
  lab <- if (length(mc)) do.call(paste, c(lapply(mc, function(cl) as.character(x[[cl]])), sep = "|")) else rep("", nrow(x))
  sc <- if ("scen" %in% names(x)) as.character(x$scen) else rep(NA_character_, nrow(x))
  in_mmm <- grepl("mmm", lab, ignore.case = TRUE)
  in_scen_ok <- is.na(sc) | sc %in% c("M3", "M4", "M5", "all", "")
  OUT[[f]] <- data.table(table = f, rows = nrow(x), differ = sum(d), differ_outside_mmm = sum(d & !in_mmm), differ_mmm_outside_M3M5 = sum(d & in_mmm & !in_scen_ok),
                         example_outside = if (any(d & !in_mmm)) lab[which(d & !in_mmm)[1]] else "", note = "")
}
R <- rbindlist(OUT, fill = TRUE); fwrite(R, "tables/V8_S11_gc3_chain.csv"); print(R[differ > 0 | note != ""])
