# s10t_extract2.R -- S10 "Testo": read-only extraction of S6b ns cc records, S7d2 interval rows, S7g quadrature rows.
.libPaths("/home/user/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressMessages({library(data.table); library(qs2)}); setDTthreads(2L)
setwd("/home/user/MMM"); OUT <- "results_v4/S10/testo"
t0 <- Sys.time(); tm <- function(s) cat(sprintf("[%.1f min] %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), s))
f <- list.files("results_v4/S6b/ns/rows", pattern = "\\.rds$", full.names = TRUE)
CC <- rbindlist(parallel::mclapply(f, function(x) { y <- readRDS(x)$cc; y[, c("A_alpha", "A_alpha_names") := NULL]; y }, mc.cores = 3L), fill = TRUE)
qs_save(CC, file.path(OUT, "x_nscc.qs2"), nthreads = 2L); tm(sprintf("ns cc %d rows", nrow(CC)))
for (d in c("results_v4/S7d2/rows", "results_v4/S7g/rows")) {
  f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  y <- readRDS(f[1]); cat("==", d, length(f), "files; first file structure:\n"); str(y, max.level = 2, list.len = 60)
  L <- parallel::mclapply(f, function(x) { y <- readRDS(x); if (is.data.frame(y)) list(main = as.data.table(y)) else lapply(y, function(z) if (is.data.frame(z)) as.data.table(z) else NULL) }, mc.cores = 3L)
  nm <- unique(unlist(lapply(L, names)))
  X <- setNames(lapply(nm, function(k) rbindlist(lapply(L, `[[`, k), fill = TRUE)), nm)
  for (k in nm) cat(d, k, nrow(X[[k]]), "rows; cols:", paste(names(X[[k]]), collapse = ","), "\n")
  qs_save(X, file.path(OUT, paste0("x_", basename(dirname(d)), ".qs2")), nthreads = 2L); rm(L, X); invisible(gc()); tm(paste("saved", d))
}
