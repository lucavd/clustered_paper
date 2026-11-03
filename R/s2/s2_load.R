# =============================================================================
# s2_load.R -- S2: read every chunk of results_full/ once (read-only) and keep
# the rows S2 needs: term == "time:trt" plus the failure rows (term NA).
# Output: results_v4/S2/rows_timetrt.rds (not tracked; md5 in the S2 report).
# =============================================================================
LIB <- path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(LIB)
suppressMessages(library(data.table))
RES <- path.expand("~/MMM/results_full")
OUT <- path.expand("~/MMM/results_v4/S2"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fs <- list.files(RES, pattern = "[.]rds$", full.names = TRUE)
stopifnot(length(fs) == 28800L)
keyre <- "^(GenBinaryY|GLMMnative|GLMMAR1)_m([1-5])_N([0-9]+)_T([0-9]+)_(mmlb|manual)_c([0-9]+)[.]rds$"
rd <- function(f) {
  x <- as.data.table(readRDS(f))
  b <- basename(f)
  x <- x[is.na(term) | term == "time:trt",
         .(method, term, estimate, se, stat, p, converged, status, sim_id = as.integer(sim_id),
           D = if ("D" %in% names(x)) D else NA_real_)]
  x[, `:=`(file_dgp = sub(keyre, "\\1", b), model = as.integer(sub(keyre, "\\2", b)),
           N = as.integer(sub(keyre, "\\3", b)), ntime = as.integer(sub(keyre, "\\4", b)),
           label = sub(keyre, "\\5", b), chunk = as.integer(sub(keyre, "\\6", b)))]
  x
}
t0 <- Sys.time()
L <- parallel::mclapply(fs, rd, mc.cores = 40, mc.preschedule = TRUE)
bad <- which(vapply(L, function(z) !is.data.table(z), logical(1)))
if (length(bad)) stop("unreadable chunks: ", paste(basename(fs[bad]), collapse = ", "))
A <- rbindlist(L)
A[, dgp := c(GenBinaryY = "GenBinaryY", GLMMnative = "GLMM-native", GLMMAR1 = "GLMM-AR1")[file_dgp]]
A[, file_dgp := NULL]
saveRDS(A, file.path(OUT, "rows_timetrt.rds"), compress = FALSE)
cat(sprintf("rows %d | files %d | %.1f s\n", nrow(A), length(fs),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(A[, .N, by = .(method)][order(method)])
print(A[is.na(term), .N, by = .(method)][order(method)])
