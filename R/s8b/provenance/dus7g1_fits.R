# dus7g1_fits.R -- S8 review, track Dubbi, doubt 1 received from S7g: number of model fits per run and method,
# counted from the per-replicate chunk files. A "fit" = one method applied to one dataset (one sim_id of one chunk);
# small-sample corrections (gee_md, gee_kc, gee_fg, gee_ss rows of the old format; se_* columns of the new one) are
# computed from an existing GEE fit and are counted separately as "corrections", not as fits. A likelihood-ratio test
# is counted as one fit row (it refits the null model: 2 model fits). Datasets = rows of recs (or distinct sim_id).
# Read-only. Single R process. Output: dus7g1_fits.csv
.libPaths(path.expand("~/MMM/renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))
suppressMessages(library(data.table)); setDTthreads(1L)
R <- "/home/user/MMM"; t0 <- Sys.time(); OUT <- list()
add <- function(run, comp, X) { X[, `:=`(run = run, component = comp)]; OUT[[length(OUT) + 1L]] <<- X }
files <- function(d) list.files(file.path(R, d), pattern = "\\.rds$", full.names = TRUE)
CORR <- c("gee_md", "gee_kc", "gee_fg", "gee_ss")
old_rows <- function(y) { y <- as.data.table(y); y[, .(n = uniqueN(sim_id)), by = method] }
tally <- function(L) rbindlist(L)[, .(n = sum(n)), by = method]
## results_full and S3 blocks (old tibble format: one row per term per method per sim_id)
for (d in c("results_full", "results_v4/S3/L", "results_v4/S3/C", "results_v4/S3/F1", "results_v4/S3/F2")) {
  f <- files(d); L <- vector("list", length(f)); nd <- 0L
  for (i in seq_along(f)) { y <- as.data.table(readRDS(f[i])); L[[i]] <- old_rows(y); nd <- nd + uniqueN(y$sim_id) }
  X <- tally(L); X <- rbind(X, data.table(method = "datasets (distinct sim_id per chunk)", n = nd)); add(sub("results_v4/", "", d), "rows", X)
  cat(d, length(f), "chunks", sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
## S6cc (MMM convergence check: fit A = refit with time trend, fit B = refit without, when at the boundary)
f <- files("results_v4/S6cc/rows"); L <- vector("list", length(f))
for (i in seq_along(f)) { y <- readRDS(f[i]); L[[i]] <- data.table(method = c("mmm fit A", "mmm fit B"), n = c(sum(!is.na(y$A_ok)), sum(!is.na(y$B_ok)))) }
add("S6cc", "cc", tally(L))
## S6b ns (rows old format + cc + recs)
f <- files("results_v4/S6b/ns/rows"); L <- list(); M <- list(); nd <- 0L
for (i in seq_along(f)) { y <- readRDS(f[i]); L[[i]] <- old_rows(y$rows)
  M[[i]] <- data.table(method = c("mmm fit A (cc)", "mmm fit B (cc)"), n = c(sum(!is.na(y$cc$A_ok)), sum(!is.na(y$cc$B_ok)))); nd <- nd + nrow(y$recs) }
add("S6b_ns", "rows", rbind(tally(L), tally(M), data.table(method = "datasets (recs)", n = nd)))
## S7c pilot (tests) and posthoc
f <- files("results_v4/S7c/pilot"); L <- list(); nd <- 0L
for (i in seq_along(f)) { y <- readRDS(f[i]); L[[i]] <- y$tests[, .(n = .N), by = .(method = paste("test", test))]; nd <- nd + nrow(y$recs) }
add("S7c_pilot", "tests", rbind(tally(L), data.table(method = "datasets (recs)", n = nd)))
f <- files("results_v4/S7c/posthoc"); nd <- 0L; for (x in f) nd <- nd + nrow(readRDS(x))
add("S7c_posthoc", "rows", data.table(method = "datasets (rows; tests recomputed)", n = nd))
## S7d (GEE three working correlations, corrections as columns; sign-flip in recs)
f <- files("results_v4/S7d/rows"); L <- list(); nd <- 0L; nf <- 0L
for (i in seq_along(f)) { y <- readRDS(f[i]); L[[i]] <- y$rows[, .(n = .N), by = .(method = paste0("gee_", cor))]; nd <- nd + nrow(y$recs); nf <- nf + sum(!is.na(y$recs$flip_n_ref) & y$recs$flip_n_ref > 0) }
add("S7d", "rows", rbind(tally(L), data.table(method = c("datasets (recs)", "sign-flip tests with reference distribution"), n = c(nd, nf))))
cat("S7d", sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
## S7d2 (inverted sign-flip interval and GEE AR(1) FG per dataset)
f <- files("results_v4/S7d2/rows"); nd <- 0L; ng <- 0L; ne <- 0
for (x in f) { y <- readRDS(x); nd <- nd + nrow(y); ng <- ng + sum(!is.na(y$ar1_gee_ok)); ne <- ne + sum(y$n_eval, na.rm = TRUE) }
add("S7d2", "rows", data.table(method = c("datasets (rows)", "gee_ar1", "sign-flip evaluations (n_eval)"), n = c(nd, ng, ne)))
## S7e (copula: study methods, convergence check, GEE three correlations)
f <- files("results_v4/S7e/rows"); L <- list(); M <- list(); G <- list(); nd <- 0L; sec <- 0
for (i in seq_along(f)) { y <- readRDS(f[i]); st <- as.data.table(y$study)
  L[[i]] <- st[, .(n = uniqueN(sim_id)), by = method]
  M[[i]] <- data.table(method = c("mmm fit A (cc)", "mmm fit B (cc)"), n = c(sum(!is.na(y$cc$A_ok)), sum(!is.na(y$cc$B_ok))))
  G[[i]] <- y$gee[, .(n = .N), by = .(method = paste0("gee_", cor))]; nd <- nd + nrow(y$recs); sec <- sec + sum(y$secs) }
add("S7e", "rows", rbind(tally(L), tally(M), tally(G), data.table(method = c("datasets (recs)", "worker seconds (sum of secs)"), n = c(nd, sec))))
## S7g (GLMM refits by nAGQ)
f <- files("results_v4/S7g/rows"); L <- list(); nd <- 0L; sec <- 0
for (i in seq_along(f)) { y <- readRDS(f[i]); L[[i]] <- y$rows[, .(n = .N), by = .(method = paste0("glmer nAGQ=", nAGQ))]; nd <- nd + nrow(y$recs); sec <- sec + sum(y$rows$secs, na.rm = TRUE) }
add("S7g", "rows", rbind(tally(L), data.table(method = c("datasets (recs)", "fit seconds (sum of secs)"), n = c(nd, sec))))
X <- rbindlist(OUT, use.names = TRUE); X[, kind := fifelse(method %in% CORR, "correction", fifelse(grepl("^datasets|seconds|evaluations|sign-flip tests", method), "count", "fit"))]
setcolorder(X, c("run", "component", "method", "kind", "n")); fwrite(X, "dus7g1_fits.csv"); print(X, nrows = 200)
print(X[kind == "fit", .(fits = sum(n)), by = run]); cat("total fits", X[kind == "fit", sum(n)], "\n")
cat(sprintf("done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
